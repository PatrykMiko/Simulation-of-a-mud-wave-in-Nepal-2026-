"""Reproducible Nepal v2 preparation, Optuna calibration and resumable jobs.

python -m ravaflow_tools.nepal --help
Only completed simulations are checkpoints; an interrupted solver run restarts.
"""
from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import shutil
import subprocess
import sys
import time

import numpy as np
from scipy import ndimage
import yaml

from . import metrics as M
from . import runner

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'configs/nepal_v2.yaml'
LAYERS = ('elevation', 'release_height', 'entrainment_height', 'reference_extent')


def write_json(path, value):
    def clean(x):
        if isinstance(x, dict):
            return {str(k): clean(v) for k, v in x.items()}
        if isinstance(x, (list, tuple)):
            return [clean(v) for v in x]
        if isinstance(x, (np.integer, np.floating)):
            x = x.item()
        if isinstance(x, float) and not math.isfinite(x):
            return None
        return x
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(clean(value), indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def load(path):
    cfg = yaml.safe_load(Path(path).read_text())
    for key in ('data_dir', 'study_dir'):
        cfg[key] = str((ROOT / cfg[key]).resolve())
    if cfg['dem'] not in ('copernicus', 'alos', 'srtm'):
        raise ValueError('Nieznany DEM')
    for key in ('coarse_factor', 'refine_factor'):
        if not isinstance(cfg[key], int) or cfg[key] < 1:
            raise ValueError(f'{key}: wymagany dodatni mnożnik siatki')
    validate_initial_trials(cfg)
    return cfg


def validate_initial_trials(cfg):
    seeds = cfg.get('initial_trials', [])
    if seeds and len(seeds) > cfg['n_trials']:
        raise ValueError('Budżet prób jest mniejszy niż liczba przygotowanych hipotez')
    ids = set()
    for seed in seeds:
        if not seed.get('id') or seed['id'] in ids:
            raise ValueError('Każda hipoteza musi mieć unikalne id')
        ids.add(seed['id'])
        if set(seed['params']) != set(cfg['search']):
            raise ValueError(f"Niekompletne parametry hipotezy {seed['id']}")
        for key, value in seed['params'].items():
            spec = cfg['search'][key]
            if not math.isfinite(value) or not spec['low'] <= value <= spec['high']:
                raise ValueError(f"Hipoteza {seed['id']}: {key} poza zakresem")


def enqueue_initial_trials(study, cfg):
    """Seed the study once, also recovering interruption halfway through initialization."""
    seeds = cfg.get('initial_trials', [])
    existing = {t.user_attrs.get('seed_id') for t in study.trials}
    for seed in seeds:
        if seed['id'] not in existing:
            study.enqueue_trial(seed['params'], user_attrs={
                'seed_id': seed['id'], 'design_rationale': seed['rationale'],
                'origin': 'designed_hypothesis'})
    if not seeds and not study.trials:
        study.enqueue_trial(cfg['baseline'])


def source_hashes(cfg):
    files = [Path(cfg['data_dir']) / cfg['dem'] / f'{x}.tif' for x in LAYERS]
    return {str(p.relative_to(cfg['data_dir'])): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}


def fingerprint(cfg):
    operational = {'n_trials', 'max_hours', 'timeout_seconds', 'final_timeout_seconds', 'keep_trial_rasters'}
    settings = {k: v for k, v in cfg.items() if k not in operational}
    # Changes to preparation, scoring or the wrapper must not silently reuse an old study.
    code = {n: hashlib.sha256((Path(__file__).parent / n).read_bytes()).hexdigest()
            for n in ('nepal.py', 'runner.py', 'metrics.py')}
    return hashlib.sha256(json.dumps([settings, source_hashes(cfg), code], sort_keys=True).encode()).hexdigest()


def check_study(cfg):
    directory = Path(cfg['study_dir'])
    directory.mkdir(parents=True, exist_ok=True)
    manifest = directory / 'study_manifest.json'
    digest = fingerprint(cfg)
    if manifest.exists():
        if json.loads(manifest.read_text())['fingerprint'] != digest:
            raise ValueError('Zmieniono dane, kod lub ustawienia fizyczne. Użyj nowego study_dir; nie mieszaj badań.')
    else:
        write_json(manifest, {'fingerprint': digest, 'config': cfg, 'sources': source_hashes(cfg)})


def block_mean(arr, factor):
    """Zero padding and block averages conserve the integral of a thickness raster."""
    if factor < 1:
        raise ValueError('factor >= 1')
    rows, cols = arr.shape
    padded = np.pad(arr, ((0, (-rows) % factor), (0, (-cols) % factor)))
    return padded.reshape(padded.shape[0] // factor, factor,
                          padded.shape[1] // factor, factor).mean(axis=(1, 3))


def coarse_meta(meta, factor):
    top, left = grid_padding(meta, factor)
    nr, nc = math.ceil((meta['nrows'] + top) / factor), math.ceil((meta['ncols'] + left) / factor)
    cs = meta['cellsize'] * factor
    north = meta['yll'] + (meta['nrows'] + top) * meta['cellsize']
    return {**meta, 'nrows': nr, 'ncols': nc, 'cellsize': cs,
            'xll': meta['xll'] - left * meta['cellsize'], 'yll': north - nr * cs, 'nodata': -9999.0}


def grid_padding(meta, factor):
    """Align to global GRASS cell multiples before aggregation, not by later resampling."""
    cs = meta['cellsize'] * factor
    north = meta['yll'] + meta['nrows'] * meta['cellsize']
    left = (meta['xll'] - math.floor(meta['xll'] / cs) * cs) / meta['cellsize']
    top = (math.ceil(north / cs) * cs - north) / meta['cellsize']
    if not np.allclose([left, top], np.round([left, top]), atol=1e-8, rtol=0):
        raise ValueError('Źródłowa siatka nie jest wyrównana do wielokrotności 30 m')
    return int(round(top)), int(round(left))


def prepare(cfg, factor):
    """Preserve source data. Extrapolate terrain outside domain without zero cliffs.

    Nearest valid elevation is a numerical extension, NOT measured terrain. Its
    affected footprint is reported for every run; results there are not scored.
    """
    import rasterio
    directory = Path(cfg['study_dir']) / 'grids' / f"{cfg['dem']}_{30 * factor}m"
    manifest = directory / 'manifest.json'
    hashes = source_hashes(cfg)
    if manifest.exists():
        saved = json.loads(manifest.read_text())
        if saved['sources'] == hashes and all((directory / f'{n}.tif').exists() for n in LAYERS):
            return directory
        raise ValueError(f'Niezgodny cache danych: {directory}')
    directory.mkdir(parents=True, exist_ok=True)
    arrays = {}
    meta = None
    for name in LAYERS:
        p = Path(cfg['data_dir']) / cfg['dem'] / f'{name}.tif'
        arr, m = M.read_raster(p)
        with rasterio.open(p) as ds:
            if ds.crs is None or ds.crs.to_epsg() != 32645:
                raise ValueError(f'Wymagany EPSG:32645: {p}')
        if meta is not None and not M.same_grid(meta, m):
            raise ValueError(f'Niezgodna siatka: {p}')
        if not np.isfinite(arr).all():
            raise ValueError(f'NaN/inf: {p}')
        arrays[name] = arr
        meta = m
    if meta['cellsize'] != 30:
        raise ValueError('Oczekiwana źródłowa siatka 30 m')
    ref = arrays['reference_extent']
    if not np.isin(ref, [-1, 0, 1]).all():
        raise ValueError('Maska referencyjna musi zawierać -1, 0, 1')
    valid = ref >= 0
    rel = arrays['release_height']
    ent = arrays['entrainment_height']
    if ((rel < 0).any() or (ent < 0).any() or (rel[~valid] != 0).any()
            or (ent[~valid] != 0).any() or not valid.any()):
        raise ValueError('Niepoprawne grubości lub źródło poza domeną')
    volume = float(rel.sum() * 900)
    if not np.isclose(volume, cfg['volume_m3'], rtol=1e-3):
        raise ValueError(f'Objętość {volume} != {cfg["volume_m3"]}')
    # Average only measured elevations in mixed edge blocks, never with zeros.
    top, left = grid_padding(meta, factor)
    def aggregate(arr):
        return block_mean(np.pad(arr, ((top, 0), (left, 0))), factor)
    coverage = aggregate(valid.astype(float))
    terrain = aggregate(np.where(valid, arrays['elevation'], 0))
    terrain = np.divide(terrain, coverage, out=np.zeros_like(terrain), where=coverage > 0)
    outside = coverage == 0
    nearest = ndimage.distance_transform_edt(outside, return_distances=False, return_indices=True)
    terrain[outside] = terrain[tuple(nearest[:, outside])]
    cm = coarse_meta(meta, factor)
    output = {'elevation': terrain, 'release_height': aggregate(rel),
              'entrainment_height': aggregate(ent),
              'reference_extent': np.where(coverage > 0,
                  aggregate((ref > 0).astype(float)) > 0, -1).astype(float)}
    # Coarse reference is for inspection only. Objective always uses native 30 m data.
    for name, arr in output.items():
        M.write_raster(directory / f'{name}.tif', arr, cm, epsg=32645)
    out_rel, _ = M.read_raster(directory / 'release_height.tif')
    after = float(out_rel.sum() * cm['cellsize']**2)
    if not np.isclose(after, volume, rtol=1e-6):
        raise ValueError('Agregacja zmieniła objętość')
    write_json(manifest, {'sources': hashes, 'factor': factor, 'shape': list(terrain.shape),
                         'volume_before_m3': volume, 'volume_after_m3': after,
                         'boundary': 'nearest measured terrain; extrapolated cells excluded from evaluation',
                         'reference': 'any-overlap preview only; objective uses native reference', 'meta': cm})
    return directory


def ensure_mapset(workdir):
    location = Path(workdir) / 'grass'
    tempdir = Path(workdir) / '.create_tmp'
    tempdir.mkdir(exist_ok=True)
    subprocess.run(['grass', '-c', 'EPSG:32645', '-e', str(location)],
                   env={**runner.grass_environment(), 'TMPDIR': str(tempdir)}, check=True, stdin=subprocess.DEVNULL,
                   stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return location / 'PERMANENT'


def model_config(cfg, params, directory, grid, mapset, duration):
    p = params
    return {'workdir': str(directory), 'grass_mapset': str(mapset), 'flags': '-e',
            'region_from': 'elev', 'timeout_seconds': cfg['timeout_seconds'],
            'inputs': {'elev': str(grid / 'elevation.tif'), 'rel': str(grid / 'release_height.tif'),
                       'ent': str(grid / 'entrainment_height.tif')},
            'params': {'prefix': 'nepal', 'phases': 3, 'elevation': 'elev',
                       'hrelease': 'rel', 'rhrelease1': p['solid_fraction'],
                       'hentrmax1': 'ent', 'rhentrmax1': cfg['entrained_solid_fraction'],
                       'centrainment': 1, 'entrainment': [p['entrainment_log10'], 0],
                       'density': [2700, 1800, 1000],
                       'friction': [p['phi'], 0, 0, p['delta'], 0, 0, 0, 0, p['manning']],
                       'deformation': [1, 1, 1], 'cmelt': 0, 'ctopo': 0, 'csurface': 0,
                       'cfl': [0.4, 0.001], 'thresholds': [0.1, 10000, 10000, 0.1, '0.000001'],
                       **({'ctrlpoints': [v for p in cfg['control_points'] for v in (p['x'], p['y'])]}
                          if cfg.get('control_points') else {}),
                       'time': [min(cfg['output_interval_s'], duration), duration]}}


def evaluate(cfg, hmax, meta):
    ref, rm = M.read_raster(Path(cfg['data_dir']) / cfg['dem'] / 'reference_extent.tif')
    rel, _ = M.read_raster(Path(cfg['data_dir']) / cfg['dem'] / 'release_height.tif')
    native = M.sample_to_grid(hmax, meta, rm) if not M.same_grid(meta, rm) else hmax
    observed, valid = M.reference_masks(ref, rm)
    protocol_valid = valid.copy()
    # Exclude the prescribed source from both fitting and held-out evaluation.
    valid &= rel == 0
    x = rm['xll'] + (np.arange(rm['ncols']) + 0.5) * rm['cellsize']
    split = float(np.quantile(np.broadcast_to(x, ref.shape)[observed & valid], cfg['holdout_fraction']))
    train = valid & (x[None, :] >= split)
    holdout = valid & (x[None, :] < split)
    tables = {}
    for threshold in (0.1, 0.5, 1.0):
        model = M.model_extent(native, {**rm, 'nodata': meta['nodata']}, threshold)
        tables[str(threshold)] = {key: M.score(model, observed, mask, rm['cellsize'])
                                 for key, mask in [('train', train), ('holdout', holdout), ('all', valid),
                                                   ('protocol_including_source', protocol_valid)]}
    domain = ref >= 0
    edge = domain & ~ndimage.binary_erosion(domain, iterations=2)
    model = native >= cfg['threshold_m']
    return {'objective': tables[str(cfg['threshold_m'])]['train']['jaccard'],
            'thresholds': tables, 'holdout_west_of_x': split,
            'boundary_affected_cells_30m': int((edge & model).sum()),
            'outside_domain_affected_cells_30m': int((~domain & model).sum()),
            'source_excluded': True}


def unique_directory(parent):
    parent.mkdir(parents=True, exist_ok=True)
    for i in range(10000):
        d = parent / f'attempt_{i:03d}'
        try:
            d.mkdir()
            return d
        except FileExistsError:
            pass
    raise RuntimeError('Za dużo prób')


def simulate(cfg, params, factor, label, duration=None, final=False):
    started = time.monotonic()
    duration = duration or cfg['duration_s']
    parent = Path(cfg['study_dir']) / label
    done = parent / 'completed.json'
    if done.exists():
        return json.loads(done.read_text())
    grid = prepare(cfg, factor)
    directory = unique_directory(parent)
    mapset = ensure_mapset(directory)
    rcfg = model_config(cfg, params, directory, grid, mapset, duration)
    if final:
        rcfg['timeout_seconds'] = cfg.get('final_timeout_seconds')
        if shutil.disk_usage(directory).free < cfg.get('final_min_free_gib', 10) * 1024**3:
            raise RuntimeError('Za mało miejsca na finał 30 m; zwolnij miejsce i wznów start/final-start')
    (directory / 'run.yaml').write_text(yaml.safe_dump(rcfg, sort_keys=False))
    write_json(Path(cfg['study_dir']) / 'current.json', {'label': label, 'directory': str(directory),
               'factor': factor, 'simulated_seconds': duration, 'started': time.time()})
    result = runner.run(rcfg)
    parameter_file = Path(result['results_dir']) / 'nepal_files/nepal_paramcomm.txt'
    applied = dict(line.split('\t', 1) for line in parameter_file.read_text().splitlines() if '\t' in line)
    if applied.get('Maximum height of entrainment MIXTURE or PHASE 1 map') != 'ent':
        raise RuntimeError('Solver nie odczytał rastra limitującego erozję')
    ascii_dir = Path(result['ascii_dir'])
    hmax, meta = runner.read_output_ascii(ascii_dir / 'nepal_hflow_max.asc')
    _, input_meta = M.read_raster(grid / 'elevation.tif')
    if not M.same_grid(meta, input_meta):
        raise RuntimeError(f'Solver zmienił siatkę wejściową: {input_meta} -> {meta}')
    if not np.isfinite(hmax).all() or np.any((hmax != meta['nodata']) & (hmax < 0)):
        raise RuntimeError('Niepoprawny wynik hflow_max')
    # A numerical-success flag alone does not prove the requested horizon was reached.
    summary = Path(result['results_dir']) / 'nepal_files/nepal_summary.txt'
    if not summary.exists():
        raise RuntimeError('Brak podsumowania kroków solvera')
    lines = [x.split() for x in summary.read_text().splitlines() if x.strip()]
    header = next((row for row in lines if 'tsum' in row), None)
    if header is None:
        raise RuntimeError('Brak czasu tsum w podsumowaniu')
    rows = [row for row in lines if len(row) == len(header) and row[0].isdigit()]
    last = dict(zip(header, map(float, rows[-1]))) if rows else {}
    if not all(math.isfinite(value) for value in last.values()):
        raise RuntimeError('NaN/inf w podsumowaniu solvera')
    if last.get('tsum', 0) < duration - 1:
        raise RuntimeError(f'Solver nie osiągnął czasu {duration}: {last}')
    metrics = evaluate(cfg, hmax, meta)
    volume_file = Path(result['results_dir']) / 'nepal_files/nepal_volumes.txt'
    volume_rows = [line.split() for line in volume_file.read_text().splitlines() if line.strip()]
    volume_last = dict(zip(volume_rows[0], map(float, volume_rows[-1])))
    volume_initial = dict(zip(volume_rows[0], map(float, volume_rows[1])))
    initial_total = sum(value for key, value in volume_initial.items() if '_vflow' in key)
    if not np.isclose(initial_total, cfg['volume_m3'], rtol=1e-3):
        raise RuntimeError(f'Solver odczytał błędną objętość początkową: {initial_total}')
    erosion = -sum(v for k, v in volume_last.items() if 'vchange' in k)
    available, am = M.read_raster(grid / 'entrainment_height.tif')
    cap = float(available.sum() * am['cellsize']**2)
    if erosion > cap * 1.01:
        raise RuntimeError(f'Erozja {erosion:g} przekracza zasób {cap:g}')
    metrics['net_entrainment_m3'] = erosion
    metrics['initial_volume_m3'] = initial_total
    metrics['available_entrainment_m3'] = cap
    metrics['volume_balance_residual_m3'] = (sum(last.get(f'volume{i}', 0) for i in (1, 2, 3))
                                              - cfg['volume_m3'] - erosion)
    # A residual can include material leaving the rectangular domain, plus numerical error.
    reach_path = ascii_dir / 'nepal_treach.asc'
    if reach_path.exists():
        reach, reach_meta = runner.read_output_ascii(reach_path)
        arrived = (reach >= 0) & (reach <= duration) & (reach != reach_meta['nodata'])
        metrics['latest_arrival_s'] = float(reach[arrived].max()) if arrived.any() else None
        points = []
        for point in cfg.get('control_points', []):
            col = math.floor((point['x'] - reach_meta['xll']) / reach_meta['cellsize'])
            row = reach_meta['nrows'] - 1 - math.floor((point['y'] - reach_meta['yll']) / reach_meta['cellsize'])
            inside = 0 <= row < reach.shape[0] and 0 <= col < reach.shape[1]
            value = float(reach[row, col]) if inside and arrived[row, col] else None
            points.append({**point, 'arrival_s_from_model_release': value})
        metrics['control_points'] = points
    report = {'params': params, 'factor': factor, 'duration_s': duration,
              'wall_seconds': time.monotonic() - started, 'solver': result, 'last_step': last,
              'metrics': metrics, 'run_config': str(directory / 'run.yaml')}
    np.savez_compressed(directory / 'hmax.npz', hmax=hmax.astype('float32'), meta=json.dumps(meta))
    write_json(directory / 'metrics.json', metrics)
    if final:
        from . import viz
        elevation, em = M.read_raster(grid / 'elevation.tif')
        if not M.same_grid(em, meta):
            elevation = M.sample_to_grid(elevation, em, meta)
        viz.plot_hmax(elevation, hmax, meta, directory / 'hmax.png', threshold=0.1)
        M.write_raster(directory / 'hmax.tif', hmax, meta, epsg=32645)
        ref, rm = M.read_raster(Path(cfg['data_dir']) / cfg['dem'] / 'reference_extent.tif')
        if not M.same_grid(rm, meta):
            ref = M.sample_to_grid(ref, rm, meta)
        obs, valid = M.reference_masks(ref, {**meta, 'nodata': rm['nodata']})
        model = M.model_extent(hmax, meta, cfg['threshold_m'])
        viz.plot_extent_comparison(elevation, model, obs, valid, meta,
                                   directory / 'extent_comparison.png',
                                   metrics=metrics['thresholds'][str(cfg['threshold_m'])]['protocol_including_source'])
    elif not cfg.get('keep_trial_rasters', False):
        # Exact, job-owned raster output only; retain configs, logs, summaries and compressed hmax.
        shutil.rmtree(ascii_dir)
        shutil.rmtree(Path(result['results_dir']) / 'nepal_vr', ignore_errors=True)
    # Each attempt owns its mapset; completed mapsets are regenerable and may be large.
    shutil.rmtree(mapset.parent)
    write_json(done, report)
    return report


def storage(cfg):
    import optuna
    return optuna.create_study(study_name='nepal_v2', direction='maximize', load_if_exists=True,
        storage='sqlite:///' + str(Path(cfg['study_dir']) / 'study.db'),
        sampler=optuna.samplers.TPESampler(seed=cfg['seed'], n_startup_trials=cfg['startup_trials']))


def sample(trial, cfg):
    return {key: trial.suggest_float(key, spec['low'], spec['high'], log=spec.get('log', False))
            for key, spec in cfg['search'].items()}


def calibrate(cfg):
    import optuna
    study = storage(cfg)
    State = optuna.trial.TrialState
    # Only one coordinator can enter this method (flock). Recover interrupted jobs.
    for trial in study.get_trials(states=(State.RUNNING,)):
        checkpoint = Path(cfg['study_dir']) / f'trials/t{trial.number:04d}/completed.json'
        if checkpoint.exists():
            r = json.loads(checkpoint.read_text())
            study.tell(trial.number, r['metrics']['objective'])
        else:
            study.tell(trial.number, state=State.FAIL)
            retry_params = trial.system_attrs.get('fixed_params', trial.params)
            if retry_params:
                study.enqueue_trial(retry_params, user_attrs={**trial.user_attrs, 'retry_of': trial.number})
    enqueue_initial_trials(study, cfg)
    deadline = time.monotonic() + cfg['max_hours'] * 3600 if cfg.get('max_hours') else math.inf
    failures = 0
    while len(study.get_trials(states=(State.COMPLETE,))) < cfg['n_trials']:
        if time.monotonic() >= deadline:
            print('Osiągnięto budżet czasu. Wznów poleceniem start.', flush=True)
            return False
        # Per-trial seed + stored history: restart does not repeat the startup RNG sequence.
        study.sampler = optuna.samplers.TPESampler(seed=cfg['seed'] + len(study.trials),
                                                 n_startup_trials=cfg['startup_trials'])
        trial = study.ask()
        params = sample(trial, cfg)
        try:
            r = simulate(cfg, params, cfg['coarse_factor'], f'trials/t{trial.number:04d}')
            trial.set_user_attr('wall_seconds', r['wall_seconds'])
            study.tell(trial, r['metrics']['objective'])
            failures = 0
        except (KeyboardInterrupt, SystemExit):
            # Preserve RUNNING for deterministic retry on the next invocation.
            raise
        except Exception as exc:
            trial.set_user_attr('error', str(exc))
            study.tell(trial, state=State.FAIL)
            failures += 1
            print(f'Próba {trial.number} nieudana: {exc}', flush=True)
            if failures >= 3:
                raise RuntimeError('Trzy kolejne błędy. Sprawdź logi przed wznowieniem.') from exc
        export_best(cfg, study)
    refine(cfg, study)
    if cfg.get('compare_dems', False):
        sensitivity(cfg)
    if cfg.get('auto_final', False):
        final_run(cfg)
    return True


def export_best(cfg, study):
    import optuna
    trials = study.get_trials(states=(optuna.trial.TrialState.COMPLETE,))
    csv_path = Path(cfg['study_dir']) / 'trials.csv'
    tmp = csv_path.with_suffix('.tmp')
    with tmp.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['trial', 'state', 'objective', 'wall_seconds', 'seed_id', *cfg['search']])
        writer.writeheader()
        for trial in study.trials:
            writer.writerow({'trial': trial.number, 'state': trial.state.name, 'objective': trial.value,
                             'wall_seconds': trial.user_attrs.get('wall_seconds'),
                             'seed_id': trial.user_attrs.get('seed_id'), **trial.params})
    tmp.replace(csv_path)
    if trials:
        best = max(trials, key=lambda t: t.value)
        write_json(Path(cfg['study_dir']) / 'best_coarse.json', {'trial': best.number,
                   'objective': best.value, 'params': best.params, 'completed': len(trials)})
        durations = [t.user_attrs['wall_seconds'] for t in trials if 'wall_seconds' in t.user_attrs]
        if durations:
            write_json(Path(cfg['study_dir']) / 'eta.json', {
                'median_trial_seconds': float(np.median(durations)),
                'remaining_coarse_seconds': max(0, cfg['n_trials'] - len(trials)) * float(np.median(durations)),
                'note': 'Sequential runs; refinement and final simulation excluded.'})


def refine(cfg, study):
    import optuna
    trials = sorted(study.get_trials(states=(optuna.trial.TrialState.COMPLETE,)),
                    key=lambda t: t.value, reverse=True)[:cfg['refine_top_k']]
    results = []
    for trial in trials:
        r = simulate(cfg, trial.params, cfg['refine_factor'], f'refine/t{trial.number:04d}')
        results.append({'trial': trial.number, **r})
    best = max(results, key=lambda r: r['metrics']['objective'])
    write_json(Path(cfg['study_dir']) / 'best_refined.json', best)
    print('Kalibracja i kontrola na 60 m zakończone.', flush=True)


def sensitivity(cfg):
    best = json.loads((Path(cfg['study_dir']) / 'best_refined.json').read_text())
    reports = {}
    for dem in ('alos', 'srtm'):
        if dem != cfg['dem']:
            alternate = {**cfg, 'dem': dem}
            reports[dem] = simulate(alternate, best['params'], cfg['refine_factor'],
                                    f"sensitivity/{dem}_t{best['trial']:04d}")
    write_json(Path(cfg['study_dir']) / 'dem_sensitivity.json', reports)


def final_run(cfg):
    best_path = Path(cfg['study_dir']) / 'best_refined.json'
    if not best_path.exists():
        raise RuntimeError('Najpierw zakończ kalibrację i etap refine (start).')
    best = json.loads(best_path.read_text())
    simulate(cfg, best['params'], 1, f"final/t{best['trial']:04d}", final=True)


def benchmark(cfg):
    reports = []
    for factor, seconds in cfg['benchmark']:
        reports.append(simulate(cfg, cfg['baseline'], factor, f'benchmark/f{factor}_{seconds}s', seconds))
        write_json(Path(cfg['study_dir']) / 'benchmark.json', reports)
    return reports


@contextmanager
def job_lock(cfg):
    directory = Path(cfg['study_dir'])
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'job.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('To badanie już pracuje. Użyj status/stop.') from None
        yield


def proc_identity(pid):
    try:
        raw = Path(f'/proc/{pid}/stat').read_text()
        fields = raw[raw.rindex(')') + 2:].split()
        if fields[0] == 'Z':
            return None
        return fields[19]  # starttime: field 22, fields starts at field 3
    except (FileNotFoundError, PermissionError):
        return None


def running(cfg):
    p = Path(cfg['study_dir']) / 'job.json'
    if not p.exists():
        return None
    job = json.loads(p.read_text())
    return job if proc_identity(job['pid']) == job['identity'] else None


def status(cfg):
    directory = Path(cfg['study_dir'])
    job = running(cfg)
    print('DZIAŁA' if job else 'ZATRZYMANE / ZAKOŃCZONE', json.dumps(job or {}, ensure_ascii=False))
    for name in ('job_result.json', 'best_coarse.json', 'eta.json'):
        if (directory / name).exists():
            print(name, (directory / name).read_text())
    if job and (directory / 'current.json').exists():
        current = json.loads((directory / 'current.json').read_text())
        print('Bieżący przebieg:', current)
        log = Path(current['directory']) / 'nepal.run.log'
        if log.exists():
            with log.open('rb') as f:
                f.seek(max(0, log.stat().st_size - 2000))
                print(f.read().decode(errors='replace')[-1200:])


def start(cfg, config_path, command):
    directory = Path(cfg['study_dir'])
    directory.mkdir(parents=True, exist_ok=True)
    if running(cfg):
        raise RuntimeError('Proces już działa')
    with (directory / 'job.log').open('a') as log:
        proc = subprocess.Popen([sys.executable, '-B', '-u', '-m', 'ravaflow_tools.nepal',
                                 command, '--config', str(Path(config_path).resolve())],
                                cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                start_new_session=True, close_fds=True)
    for _ in range(50):
        if proc.poll() is not None:
            raise RuntimeError(f'Nie udało się uruchomić. Zobacz {directory / "job.log"}')
        if running(cfg):
            print(f'Uruchomiono PID {proc.pid}. Log: {directory / "job.log"}')
            return
        time.sleep(0.1)
    raise RuntimeError('Brak potwierdzenia startu; sprawdź status/job.log')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['prepare', 'benchmark', 'benchmark-start', 'run', 'start',
                    'refine', 'sensitivity', 'final', 'final-start', 'status', 'stop'])
    ap.add_argument('--config', default=str(DEFAULT))
    args = ap.parse_args()
    cfg = load(args.config)
    os.environ.setdefault('MPLCONFIGDIR', str(Path(cfg['study_dir']) / '.matplotlib'))
    if args.command == 'status':
        status(cfg)
        return
    if args.command == 'stop':
        job = running(cfg)
        if job:
            os.kill(job['pid'], signal.SIGTERM)
            print('Wysłano stop. Ukończone próby zostają; bieżąca będzie powtórzona po wznowieniu.')
        else:
            print('Brak działającego procesu.')
        return
    if args.command in ('start', 'final-start', 'benchmark-start'):
        start(cfg, args.config, {'start': 'run', 'final-start': 'final', 'benchmark-start': 'benchmark'}[args.command])
        return
    def interrupted(signum, frame):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, interrupted)
    signal.signal(signal.SIGINT, interrupted)
    with job_lock(cfg):
        check_study(cfg)
        directory = Path(cfg['study_dir'])
        write_json(directory / 'job.json', {'pid': os.getpid(), 'identity': proc_identity(os.getpid()),
                   'command': args.command, 'started': time.time()})
        write_json(directory / 'job_result.json', {'state': 'running', 'command': args.command})
        complete = True
        try:
            if args.command == 'prepare':
                for factor in sorted({1, cfg['coarse_factor'], cfg['refine_factor']}):
                    print(prepare(cfg, factor))
            elif args.command == 'benchmark':
                benchmark(cfg)
            elif args.command == 'run':
                complete = calibrate(cfg)
            elif args.command == 'refine':
                refine(cfg, storage(cfg))
            elif args.command == 'final':
                final_run(cfg)
            elif args.command == 'sensitivity':
                sensitivity(cfg)
        except KeyboardInterrupt:
            write_json(directory / 'job_result.json', {'state': 'stopped', 'command': args.command})
            print('Zatrzymano. Zapisane próby można wznowić.', flush=True)
        except Exception as exc:
            write_json(directory / 'job_result.json', {'state': 'error', 'error': str(exc), 'command': args.command})
            raise
        else:
            write_json(directory / 'job_result.json', {'state': 'finished' if complete else 'budget_exhausted',
                                                       'command': args.command})
        finally:
            (directory / 'job.json').unlink(missing_ok=True)


if __name__ == '__main__':
    main()

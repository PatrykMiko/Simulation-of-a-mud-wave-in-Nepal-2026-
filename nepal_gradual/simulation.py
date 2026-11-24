"""Gradual-input simulation; isolated from checkpointed legacy packages.

Output/metric pipeline adapted from ravaflow_tools.nepal.simulate. Only the
release forcing and its initial-volume audit change; solver remains original.
"""
import json
import math
from pathlib import Path
import shutil
import time
import numpy as np
import yaml
from ravaflow_tools import nepal as N, metrics as M, runner
from ravaflow_tools.nepal import prepare, unique_directory, ensure_mapset, write_json, evaluate


def model_config(cfg, params, directory, grid, mapset, duration):
    start = float(cfg['release_start_s'])
    span = float(params['release_duration_s'])
    if not (math.isfinite(start) and math.isfinite(span) and start == 0 and 0 < span < duration):
        raise ValueError('Wymagany start 0 i 0 < czas uwalniania < horyzont')
    height, meta = M.read_raster(grid / 'release_height.tif')
    mask = height > 0
    # Native extrusion with trelease=0 initializes one second-equivalent of rate as height.
    # Account for this initial seed explicitly, so seed + integral stays at V.
    seed_seconds = 1.0
    rate = height / (span + seed_seconds)
    directory.mkdir(parents=True, exist_ok=True)
    paths = {}
    for name, values in [('release_rate', rate),
                         ('release_start', np.where(mask, start, 0.0)),
                         ('release_stop', np.where(mask, start + span, 0.0))]:
        paths[name] = str(directory / (name + '.tif'))
        M.write_raster(paths[name], values, meta, epsg=32645)
    actual, _ = M.read_raster(paths['release_rate'])
    volume = float(np.sum(actual, dtype=np.float64) * meta['cellsize']**2 * (span + seed_seconds))
    if not np.isclose(volume, cfg['volume_m3'], rtol=1e-5):
        raise RuntimeError('Całka tempa uwalniania nie zachowuje objętości źródła')
    write_json(directory / 'release_schedule.json', {
        'mode': 'extrusion', 'rate_units': 'm/s', 'start_s': start,
        'stop_s': start + span, 'duration_s': span,
        'integrated_volume_m3': volume, 'mean_discharge_m3_s': volume/(span+seed_seconds),
        'initial_seed_volume_m3': volume*seed_seconds/(span+seed_seconds),
        'initial_seed_fraction': seed_seconds/(span+seed_seconds),
        'phase1_fraction': params['solid_fraction'], 'phase2_fraction': 0,
        'phase3_fraction': 1-params['solid_fraction'],
        'note': 'Prescribed distributed source, not a measured breach hydrograph; solver integration has timestep error.'})
    result = N.model_config(cfg, params, directory, grid, mapset, duration)
    result['inputs'].pop('rel')
    result['inputs'].update(paths)
    result['params'].update(hrelease='release_rate', trelease='release_start', trelstop='release_stop')
    return result


def simulate(cfg, params, factor, label, duration=None, final=False):
    started = time.time()
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
    release = json.loads((directory / 'release_schedule.json').read_text())
    if not np.isclose(initial_total, release['initial_seed_volume_m3'], rtol=1e-3):
        raise RuntimeError(f'Błędna mała masa początkowa: {initial_total}')
    for key, expected in [('Release time map', 'release_start'), ('Release stop time map', 'release_stop')]:
        if applied.get(key) != expected:
            raise RuntimeError(f'Solver nie zastosował mapy czasu: {key}')
    metrics['release_schedule'] = release
    metrics['prescribed_source_volume_m3'] = release['integrated_volume_m3']
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
              'wall_seconds': time.time() - started, 'solver': result, 'last_step': last,
              'metrics': metrics, 'run_config': str(directory / 'run.yaml')}
    np.savez_compressed(directory / 'hmax.npz', hmax=hmax.astype('float32'), meta=json.dumps(meta))
    write_json(directory / 'metrics.json', metrics)
    if final:
        from ravaflow_tools import viz
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

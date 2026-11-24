"""One independent GRASS mapset and output tree per candidate."""
import csv
import json
import math
import shlex
import shutil
import time
from pathlib import Path

import numpy as np
from ravaflow_tools import nepal as N, runner, metrics as M
from nepal_gradual.simulation import model_config as gradual_config
from .build import PREFIX, ROOT
from . import objective as O

OUT = ROOT/'wyniki'
CAP = 43200
INTERVAL = 120
# Centre of the observed 30 m pixel nearest the downstream river-axis endpoint.
OUTLET = (300225.,3076875.)


def config():
    cfg=N.load(N.DEFAULT)
    cfg.update(study_dir=str(OUT/'work'), timeout_seconds=None, release_start_s=0,
               output_interval_s=INTERVAL, dynamics_objective=O.load())
    return cfg


def simulate(cfg, params, directory, grid, cap=CAP, outlet=OUTLET):
    frames=directory/'frames'
    frames.mkdir(parents=True)
    mapset=N.ensure_mapset(directory)
    builder=gradual_config if params['mode']=='gradual' else N.model_config
    rc=builder(cfg,params,directory,grid,mapset,cap)
    observed, rm = M.read_raster(grid/'reference_extent.tif')
    source, _ = M.read_raster(grid/'release_height.tif')
    mask = directory/'valid.u8'
    ((observed >= 0) & (source == 0)).astype('uint8').tofile(mask)
    del observed, source
    rc['pre_commands']=['export PATH='+shlex.quote(str(PREFIX))+':"$PATH"',
                        'export NEPAL_VALID='+shlex.quote(str(mask)),
                        'export NEPAL_FRAME_DIR='+shlex.quote(str(frames)),
                        f'export NEPAL_OUTLET_X={outlet[0]} NEPAL_OUTLET_Y={outlet[1]}']
    spec=cfg.get('dynamics_objective')
    if spec:
        samples=O.prepare(spec,grid,directory)
        rc['pre_commands'] += ['export NEPAL_DYNAMICS_SPEC='+shlex.quote(str(samples)),
                               'export NEPAL_DYNAMICS_CSV='+shlex.quote(str(directory/'dynamics.csv'))]
    N.write_json(directory/'run.json',rc)
    began = time.time()
    try:
        result=runner.run(rc)
    except RuntimeError:
        with (directory/'nepal.run.log').open('rb') as log:
            log.seek(max(0,log.seek(0,2)-4500))
            print(log.read().decode(errors='replace'),flush=True)
        raise
    fields=('frame','time','cells','outlet_cells','arrival','reason')
    with (frames/'timeline.csv').open() as f:
        timeline=list(csv.DictReader(f,fieldnames=fields))
    last=timeline[-1]
    if last['reason'] not in ('flow_stopped','outlet_stable','stalled','safety_cap'):
        raise RuntimeError(f'Brak potwierdzenia zakończenia: {last}')
    raster=Path(result['ascii_dir'])/'nepal_hflow_max.asc'
    h,meta=runner.read_output_ascii(raster)
    _,rm=M.read_raster(grid/'elevation.tif')
    if not M.same_grid(meta,rm) or not np.isfinite(h).all() or ((h<0)&(h!=meta['nodata'])).any():
        raise ValueError('Błędny raster końcowy')
    mask.unlink()
    metrics=N.evaluate(cfg,h,meta)['thresholds']['0.1']
    iou=metrics['all']['jaccard']
    if not math.isfinite(iou):
        raise ValueError('Nieokreślony IoU')
    objective=O.evaluate(directory/'dynamics.csv',iou,spec) if spec else {}
    # Export portable final fields, then remove redundant GRASS/ASCII.
    import rasterio
    from rasterio.transform import from_origin
    for field in ('hflow_max','treach','vfront','vhmax','basechange_fin'):
        src=raster.parent/f'nepal_{field}.asc'
        arr,m=runner.read_output_ascii(src)
        if not M.same_grid(meta,m):raise ValueError(field)
        cs=m['cellsize']
        with rasterio.open(directory/f'{field}.tif','w',driver='GTiff',width=m['ncols'],
                height=m['nrows'],count=1,dtype='float32',crs='EPSG:32645',nodata=m['nodata'],
                transform=from_origin(m['xll'],m['yll']+cs*m['nrows'],cs,cs),
                compress='deflate',predictor=3) as dst:
            dst.write(arr.astype('float32'),1)
    for name in ('summary','volumes','paramcomm','ctrlpoints'):
        src=Path(result['results_dir'])/'nepal_files'/f'nepal_{name}.txt'
        if src.exists():shutil.copy2(src,directory/src.name)
    table = (directory/'nepal_volumes.txt').read_text().splitlines()
    headings = table[0].split()
    initial = dict(zip(headings,map(float,table[1].split()),strict=True))
    final = dict(zip(headings,map(float,table[-1].split()),strict=True))
    volume = sum(v for k,v in initial.items() if '_vflow' in k)
    expected = cfg['volume_m3']
    applied = dict(line.split('\t',1) for line in (directory/'nepal_paramcomm.txt').read_text().splitlines() if '\t' in line)
    if applied.get('Maximum height of entrainment MIXTURE or PHASE 1 map') != 'ent':
        raise ValueError('Solver nie potwierdził mapy erozji.')
    if params['mode']=='gradual':
        schedule=json.loads((directory/'release_schedule.json').read_text())
        expected=schedule['initial_seed_volume_m3']
        if applied.get('Release time map')!='release_start' or applied.get('Release stop time map')!='release_stop':
            raise ValueError('Solver nie potwierdził map czasu uwalniania.')
    if not np.isclose(volume, expected,rtol=1e-3) or not all(map(math.isfinite,final.values())):
        raise ValueError('Nieprawidłowa objętość początkowa lub bilans solvera.')
    audit=dict(prescribed_m3=cfg['volume_m3'],initial_m3=volume,
               entrainment_net_m3=-sum(v for k,v in final.items() if 'vchange' in k))
    # Optional outputs not needed after validated portable exports.
    shutil.rmtree(directory/'grass')
    shutil.rmtree(result['results_dir'])
    return dict(**objective,params=params,iou_all=iou,metrics=metrics,time_s=float(last['time']),
                reached=float(last['arrival'])>=0,arrival_s=float(last['arrival']),
                end_reason=last['reason'],volume=audit,elapsed_s_including_pauses=time.time()-began)


def worker(number):
    directory=OUT/'work'/f't{number:02}'
    cfg=config()
    params=json.loads((OUT/'plans'/f'{number:02}.json').read_text())
    grid=Path(cfg['study_dir'])/'grids/copernicus_30m'
    result=simulate(cfg,params,directory,grid)
    N.write_json(directory/'report.json',dict(number=number,**result))

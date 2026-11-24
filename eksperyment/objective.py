"""One versioned objective for Optuna, retention and the published ranking.

The literature targets are explicit proxies with experimental tolerances.
They must not be described as measured confidence intervals or validation.
"""
import csv
import json
import math
from pathlib import Path

import numpy as np

from .build import ROOT

PATH = ROOT/'configs/nepal30_objective.json'


def load(path=PATH):
    spec = json.loads(Path(path).read_text())
    validate(spec)
    return spec


def validate(spec):
    if spec['crs'] != 'EPSG:32645' or not spec['version'] or not spec['points']:
        raise ValueError('Niekompletny protokół dynamiki.')
    for key in ('sample_interval_s', 'wet_depth_m', 'speed_window_s'):
        if not math.isfinite(spec[key]) or spec[key] <= 0:
            raise ValueError(key)
    for key in ('min_wet_cells', 'min_speed_samples'):
        if type(spec[key]) is not int or spec[key] < 1:
            raise ValueError(key)
    if spec['speed_window_s'] < spec['sample_interval_s']*(spec['min_speed_samples']-1):
        raise ValueError('Okno prędkości jest za krótkie.')
    if len(spec['points']) > 32:
        raise ValueError('Zbyt wiele punktów.')
    ids = set()
    for p in spec['points']:
        if not p['id'] or p['id'] in ids:
            raise ValueError('Powtórzony identyfikator punktu.')
        ids.add(p['id'])
        for key in ('x', 'y', 'radius_m', 'speed_target_m_s', 'arrival_target_s'):
            if not math.isfinite(p[key]) or (key not in ('x','y') and p[key] <= 0):
                raise ValueError(key)
        for key in ('speed_tolerance_factor', 'arrival_tolerance_factor'):
            if not math.isfinite(p[key]) or p[key] <= 1:
                raise ValueError(key)
        for key in ('location_basis', 'speed_basis', 'arrival_basis', 'uncertainty_basis'):
            if not p.get(key):
                raise ValueError(f'Brak uzasadnienia {key}')


def prepare(spec, grid, directory):
    """Freeze valid native-grid samples, excluding the prescribed source."""
    from ravaflow_tools import metrics as M
    validate(spec)
    ref, meta = M.read_raster(grid/'reference_extent.tif')
    source, sm = M.read_raster(grid/'release_height.tif')
    if not M.same_grid(meta,sm):
        raise ValueError('Niezgodna siatka punktów.')
    valid = np.isfinite(ref) & (ref >= 0) & (source == 0)
    cs = meta['cellsize']
    x = meta['xll']+(np.arange(meta['ncols'])+.5)*cs
    y = meta['yll']+(meta['nrows']-np.arange(meta['nrows'])-.5)*cs
    lines = [f"{len(spec['points'])} {spec['sample_interval_s']} {spec['wet_depth_m']}"]
    for p in spec['points']:
        mask = ((x[None,:]-p['x'])**2+(y[:,None]-p['y'])**2 <= p['radius_m']**2) & valid
        cells = np.flatnonzero(mask)
        if not spec['min_wet_cells'] <= len(cells) <= 4096:
            raise ValueError(f"Punkt {p['id']}: brak poprawnej lokalizacji w domenie lub za duży promień.")
        lines.append(str(len(cells))+' '+' '.join(map(str,cells)))
    path = directory/'dynamics-spec.txt'
    path.write_text('\n'.join(lines)+'\n')
    return path


def quality(value, target, factor):
    if value is None or value <= 0 or not math.isfinite(value):
        return 0., False
    distance = abs(math.log(value/target))/math.log(factor)
    return math.exp(-.5*distance*distance), target/factor <= value <= target*factor


def constraint(value, target, factor):
    """Nonpositive means feasible; a missing observation is explicitly violated."""
    if value is None or value <= 0 or not math.isfinite(value):
        return 2.
    error=abs(math.log(value/target))/math.log(factor)-1
    # Match the inclusive bounds used by the published ranking, including roundoff.
    return min(0.,error) if target/factor <= value <= target*factor else max(1e-12,error)


def score(iou, summaries, spec):
    """Feasible trials always win; among them rank strictly by full-domain IoU."""
    validate(spec)
    if not math.isfinite(iou) or not 0 <= iou <= 1 or len(summaries)!=len(spec['points']):
        raise ValueError('Niepoprawne dane funkcji celu.')
    details=[]
    constraints=[]
    for p, raw in zip(spec['points'],summaries,strict=True):
        if raw['id']!=p['id']:
            raise ValueError('Niezgodne punkty dynamiki.')
        aq, ap = quality(raw['arrival_s'],p['arrival_target_s'],p['arrival_tolerance_factor'])
        vq, vp = quality(raw['speed_m_s'],p['speed_target_m_s'],p['speed_tolerance_factor'])
        enough=raw['speed_samples'] >= spec['min_speed_samples'] and raw['window_complete']
        if not enough:
            vq, vp = 0., False
        constraints.extend([constraint(raw['arrival_s'],p['arrival_target_s'],p['arrival_tolerance_factor']),
                            constraint(raw['speed_m_s'],p['speed_target_m_s'],p['speed_tolerance_factor']) if enough else 2.])
        details.append({**raw,'arrival_quality':aq,'speed_quality':vq,
                        'pass':bool(ap and vp and enough)})
    dynamics_quality=sum((p['arrival_quality']+p['speed_quality'])/2 for p in details)/len(details)
    passed=all(p['pass'] for p in details)
    # No-arrival cases still teach TPE about spatial fit. A good final footprint
    # at an implausible time can never outrank any screening-qualified trial.
    value=.5+.5*iou if passed else .49*iou*(.8+.2*dynamics_quality)
    return dict(objective_version=spec['version'],objective_score=value,
                dynamics_constraints=constraints,
                dynamics_pass=passed,dynamics_quality=dynamics_quality,
                model_status='screening_pass' if passed else 'no_acceptable_dynamics',
                dynamics=details)


def evaluate(path, iou, spec):
    samples=[[] for _ in spec['points']]
    last=[-math.inf for _ in samples]
    with Path(path).open() as stream:
        for row in csv.reader(stream):
            if len(row)!=5:
                raise ValueError('Niekompletny zapis dynamiki.')
            t, index, wet, depth_sum, numerator = map(float,row)
            if not all(map(math.isfinite,(t,index,wet,depth_sum,numerator))):
                raise ValueError('Niefinitywna dynamika.')
            if index != int(index) or not 0 <= index < len(samples) or wet != int(wet):
                raise ValueError('Błędny punkt lub liczba komórek.')
            i=int(index)
            if t <= last[i] or min(t,wet,depth_sum,numerator)<0:
                raise ValueError('Niepoprawna kolejność lub wartości dynamiki.')
            last[i]=t
            samples[i].append((t,int(wet),depth_sum,numerator))
    summaries=[]
    for p, rows in zip(spec['points'],samples,strict=True):
        if not rows:
            raise ValueError(f"Brak szeregu dla {p['id']}")
        wet=[r for r in rows if r[1]>=spec['min_wet_cells'] and r[2]>0]
        arrival=wet[0][0] if wet else None
        window=[r for r in wet if r[0] <= arrival+spec['speed_window_s']] if wet else []
        # Mean of spatially depth-weighted phase speed, using only the first wave.
        speed=sum(r[3]/r[2] for r in window)/len(window) if window else None
        summaries.append(dict(id=p['id'],arrival_s=arrival,speed_m_s=speed,
            speed_samples=len(window),window_complete=bool(wet and rows[-1][0]>=arrival+spec['speed_window_s'])))
    return score(iou,summaries,spec)


def ranking_key(report):
    if report.get('objective_version') != load()['version']:
        raise ValueError('Stary ranking. Nowy protokół wymaga oddzielnej serii, bez starych checkpointów.')
    return (-report['objective_score'],-report['iou_all'],report['number'])

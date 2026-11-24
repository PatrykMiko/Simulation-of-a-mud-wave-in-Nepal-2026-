"""Four persistent batches. Only the parent publishes results and removes outputs."""
import csv
import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

from ravaflow_tools import checkpoint as C, nepal as N
from .design import batch
from .run import OUT, config
from .render import render
from . import objective as O


def records(root):
    return [C.read(p) for p in sorted((root/'records').glob('*.json'))]


def refresh(root):
    rows = sorted(records(root), key=O.ranking_key)
    if not rows:
        return
    best = rows[0]
    winner = root/'artifacts'/f"t{best['number']:02}"
    if not (winner/'report.json').is_file():
        raise RuntimeError(f'Brakuje plików lidera {winner}')
    C.atomic(root/'best.json', {**best, 'directory': str(winner.relative_to(root))})
    keys = ['number', 'objective_score', 'dynamics_pass', 'model_status', 'iou_all',
            'dynamics_quality', 'reached', 'arrival_s', 'time_s', 'end_reason']
    tmp = root/'ranking.tmp'
    with tmp.open('w') as f:
        writer = csv.DictWriter(f, keys, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(root/'ranking.csv')
    for row in rows[1:]:
        loser = root/'artifacts'/f"t{row['number']:02}"
        if loser.exists():
            shutil.rmtree(loser)


def accept(root, number):
    name = f't{number:02}'
    staged = root/'artifacts'/name
    if not staged.exists():
        (root/'work'/name).rename(staged)
    report = C.read(staged/'report.json')
    if report['number'] != number or report['iou_all'] != report['metrics']['all']['jaccard']:
        raise ValueError(f'Niespójny raport {number}')
    recomputed=O.score(report['iou_all'],report['dynamics'],O.load())
    if any(report.get(k)!=recomputed[k] for k in recomputed):
        raise ValueError('Niespójna funkcja celu lub inny protokół dynamiki.')
    C.atomic(root/'records'/f'{number:02}.json', report)
    refresh(root)


def launch(number):
    subprocess.run([sys.executable,'-B','-u','-m','eksperyment.cli','worker',str(number)],check=True)


def finish(root, cfg):
    best = C.read(root/'best.json')
    folder = root/best['directory']
    render(folder,Path(cfg['data_dir'])/cfg['dem'])
    tmp = root/'paczka.tmp'
    files = [root/'ranking.csv', root/'best.json', root/'protocol.json',
             *sorted((root/'records').glob('*.json')), *sorted(folder.rglob('*'))]
    with ZipFile(tmp,'w',ZIP_DEFLATED,compresslevel=3,allowZip64=True) as archive:
        for doc in ('SYMULACJA.md','OPIS.md'):
            archive.write(C.ROOT/doc,'protokol/'+doc)
        for p in files:
            if p.is_file():
                archive.write(p,p.relative_to(root))
    with ZipFile(tmp) as archive:
        if bad := archive.testzip():
            raise IOError(f'Błąd CRC w ZIP: {bad}')
    with tmp.open('rb') as stream:
        os.fsync(stream.fileno())
    tmp.replace(root/'paczka.zip')
    C.atomic(root/'finished.json',dict(trials=20,winner=best['number'],iou_all=best['iou_all'],
        objective_score=best['objective_score'],dynamics_pass=best['dynamics_pass'],model_status=best['model_status']))
    if not best['dynamics_pass']:
        print('BRAK ZAAKCEPTOWANEJ DYNAMIKI. Zachowano najlepszą próbę diagnostyczną, nie zatwierdzony model.',flush=True)


def study(root=OUT):
    cfg = config()
    workers = int(os.environ.get('NEPAL_WORKERS','5'))
    if not 1 <= workers <= 5:
        raise ValueError('NEPAL_WORKERS musi być w zakresie 1–5.')
    for name in ('work','artifacts','plans','records'):
        (root/name).mkdir(parents=True,exist_ok=True)
    grid=N.prepare(cfg,1)
    from .cli import signature
    from .run import CAP, INTERVAL, OUTLET
    spec=O.load()
    protocol=root/'protocol.json'
    if protocol.exists() and C.read(protocol).get('objective')!=spec:
        raise ValueError('Zmieniono protokół. Zachowaj starą serię i użyj pustego katalogu wyniki.')
    C.atomic(protocol,dict(objective=spec,threshold_m=.1,
        source_excluded=True,concurrent_workers=workers,cap_s=CAP,interval_s=INTERVAL,outlet_epsg32645=OUTLET,
        grid=C.read(grid/'manifest.json'),environment=signature()))
    for first in (1,6,11,16):
        prior = [r for r in records(root) if r['number'] < first]
        if len(prior) != first-1:
            raise RuntimeError('Niekompletna poprzednia seria.')
        proposals = batch(first, prior)
        pending = []
        for number, params in enumerate(proposals,first):
            plan = root/'plans'/f'{number:02}.json'
            if not plan.exists():
                C.atomic(plan,params)
            if (root/'records'/plan.name).exists():
                continue
            name=f't{number:02}'
            if any((root/base/name/'report.json').is_file() for base in ('artifacts','work')):
                accept(root,number)
            else:
                # Only the exclusive study manager may retry its own unfinished trial.
                scratch=root/'work'/name
                if scratch.exists():
                    shutil.rmtree(scratch)
                pending.append(number)
        errors=[]
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures={pool.submit(launch,n):n for n in pending}
            for f in as_completed(futures):
                n=futures[f]
                try:
                    f.result()
                except subprocess.CalledProcessError:
                    errors.append(n)
                    continue
                accept(root,n)
                row=C.read(root/'records'/f'{n:02}.json')
                print(f'Próba {n}: IoU={row["iou_all"]:.6f}, cel={row["objective_score"]:.6f}, dynamika={row["dynamics_pass"]}',flush=True)
        if errors:
            raise RuntimeError(f'Błędy prób {errors}. Szczegóły w logach solvera.')
    refresh(root)
    finish(root,cfg)

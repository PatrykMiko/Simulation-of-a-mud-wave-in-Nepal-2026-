"""Exercise a planned save and restore on the target host with a synthetic solver."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
import tempfile
import time

from ravaflow_tools import checkpoint as C, nepal as N, synthetic
from .build import PREFIX, ROOT
from .design import BASE
from .run import config,simulate


def sample(root, name):
    cfg=config()
    cfg.update(data_dir=str(root/'data'),study_dir=str(root/'work'),volume_m3=200000,
               output_interval_s=5,release_start_s=0,control_points=[])
    p={**BASE,'mode':'gradual','release_duration_s':4}
    result=simulate(cfg,p,root/name,root/'work/grids/copernicus_30m',cap=90,outlet=(340750,3111350))
    C.atomic(root/name/'done.json', {k:result[k] for k in ('iou_all','metrics','volume','time_s','end_reason')})


def wait(session, condition):
    deadline=time.monotonic()+600
    while not condition():
        if time.monotonic()>deadline or not C.alive(session.state().get('coordinator')):
            raise RuntimeError(f'Brak wyniku testu. Sprawdź {session.directory}/job.log')
        time.sleep(.05)


def verify():
    # Keep failed diagnostics. Remove only a fully completed verification tree.
    import shutil
    root=Path(tempfile.mkdtemp(prefix='test-restore-',dir=PREFIX))
    cfg=config()
    cfg.update(data_dir=str(root/'data'),study_dir=str(root/'work'),volume_m3=200000)
    synthetic.generate(root/'data/copernicus',volume=200000,nrows=96,ncols=96,xll=339990,yll=3110010)
    N.prepare(cfg,1)
    sample(root,'continuous')
    s=C.Session(root/'checkpoint')
    env=dict(test=1,path=str(root))
    s.start(['/usr/bin/env','NEPAL_RELATIVE_SLEEP=1',
             f'PYTHONPATH={ROOT}/nepal_gradual/bootstrap:{ROOT}',sys.executable,'-Bu','-m',
             'eksperyment.check_restore',str(root),'group'],env)
    try:
        names = ('resumed1','resumed2')
        wait(s,lambda:all((root/name/'frames/timeline.csv').exists() for name in names))
        s.save_stop()
        s.start([],env)
        wait(s,lambda:all((root/name/'done.json').exists() for name in names))
        deadline=time.monotonic()+60
        while C.alive(s.state().get('coordinator')) and time.monotonic()<deadline:
            time.sleep(.1)
        if C.alive(s.state().get('coordinator')):
            raise RuntimeError('Proces testowy nadal działa.')
        def hashes(name):
            d=root/name
            return {str(p.relative_to(d)):C.digest(p) for p in [*d.glob('*.tif'),*(d/'frames').glob('*'),d/'done.json']}
        expected = hashes('continuous')
        if any(hashes(name)!=expected for name in names):
            raise ValueError(f'Niezgodne wyniki po wznowieniu. Zachowano {root}')
    except BaseException:
        if C.alive(s.state().get('coordinator')):
            s.command(s.state(),'--kill')
        raise
    shutil.rmtree(root)
    print('Zgodne rastry, klatki, metryki i bilans dwóch równoległych prób po wznowieniu.')


if __name__=='__main__':
    if sys.argv[2]=='group':
        def run_peer(name):
            subprocess.run([sys.executable,'-Bu','-m','eksperyment.check_restore',sys.argv[1],name],check=True)
        with ThreadPoolExecutor(max_workers=2) as pool:
            list(pool.map(run_peer,('resumed1','resumed2')))
    else:
        sample(Path(sys.argv[1]),sys.argv[2])

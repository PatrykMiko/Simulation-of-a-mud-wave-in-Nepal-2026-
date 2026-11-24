"""Background controller for the 30 m Nepal experiment on Linux."""
import argparse
import os
import shutil
import subprocess
import sys
import time

import optuna
from nepal_final import disk
from ravaflow_tools import checkpoint as C, nepal as N
from .build import PREFIX, ROOT
from .run import OUT, config
from .resources import resources, describe, worker_count

SESSION = C.Session(OUT/'checkpoint')
GIB = disk.GIB


def signature():
    paths = [*sorted((ROOT/'eksperyment').glob('*.py')), ROOT/'eksperyment/monitor.c.inc',
             ROOT/'nepal_final/disk.py', ROOT/'nepal_gradual/simulation.py',
             ROOT/'nepal_gradual/bootstrap/sitecustomize.py', ROOT/'configs/nepal_v2.yaml',
             ROOT/'configs/nepal_v2_final30.yaml', PREFIX/'r.avaflow.main', PREFIX/'r.avaflow.40G']
    return dict(environment=C.environment_signature(), optuna=optuna.__version__,
                inputs=N.source_hashes(config()), files={str(p.relative_to(ROOT)):C.digest(p) for p in paths})


def status():
    if (OUT/'finished.json').exists():
        print('Zakończono.',C.read(OUT/'finished.json'),OUT/'paczka.zip')
        return
    SESSION.status()
    if (OUT/'guard.json').is_file():
        print('Nadzór:',C.read(OUT/'guard.json'))
    if (OUT/'frozen.json').exists():
        print('Pauza w RAM. Zwolnij miejsce i użyj save-stop przed wyłączeniem komputera.')
    print('Zakończone:',len(list((OUT/'records').glob('*.json'))),'/ 20')
    for path in sorted((OUT/'work').glob('t*/frames/timeline.csv')):
        with path.open() as stream:
            lines=stream.read().splitlines()
        if lines:
            print(path.parents[1].name, 'ostatnia klatka (s):',lines[-1].split(',')[1])
    if (OUT/'best.json').exists():
        best=C.read(OUT/'best.json')
        print('Najlepsza próba:',best['number'],'IoU całości:',best['iou_all'])


def save(reason):
    with C.control_lock(OUT/'operations'):
        state=SESSION.state()
        if state['state']=='saved':
            print('Zapis potwierdzony. Możesz wyłączyć komputer.')
            return
        if state['state'] not in ('running','saving'):
            raise RuntimeError('Sprawdź status i poczekaj na zakończenie startu.')
        if state['state']=='running':
            peers=disk.members(OUT)
            needed=disk.estimate(OUT,peers,2*GIB)
            disk.release_reserve(OUT)
            if shutil.disk_usage(OUT).free < needed:
                disk.freeze(OUT,state,peers,reason)
                raise RuntimeError('Pauza w RAM. Zwolnij miejsce i ponów save-stop.')
            disk.thaw(OUT)
        SESSION.save_stop()
        C.atomic(OUT/'guard.json',dict(state='saved',reason=reason))


def supervise():
    with C.control_lock(OUT/'supervision'):
        state=SESSION.state()
        try:
            if state['state'] in ('new','saved'):
                if shutil.disk_usage(OUT).free < 4*GIB:
                    raise RuntimeError('Zwolnij miejsce przed startem (poniżej 4 GiB).')
                disk.allocate_reserve(OUT,2*GIB)
                bootstrap=ROOT/'nepal_gradual/bootstrap'
                argv=['/usr/bin/env','NEPAL_RELATIVE_SLEEP=1',f'PYTHONPATH={bootstrap}:{ROOT}',
                      sys.executable,'-B','-u','-m','eksperyment.cli','study']
                SESSION.start(argv,signature(),startup_timeout=3600)
                disk.prune_consumed(OUT,state)
            while True:
                state=SESSION.state()
                if state['state']=='saved':
                    return
                if state['state']=='saving':
                    try:
                        save('checkpoint')
                    except RuntimeError as exc:
                        if 'Inne polecenie sterujące' not in str(exc):
                            raise
                        time.sleep(2)
                        continue
                    return
                if not C.alive(state.get('coordinator')):
                    if (OUT/'finished.json').exists():
                        disk.release_reserve(OUT)
                        C.atomic(OUT/'guard.json',dict(state='finished'))
                        return
                    raise RuntimeError('Zadanie przerwane. Zobacz wyniki/checkpoint/job.log.')
                if (OUT/'frozen.json').exists():
                    return
                peers=disk.members(OUT)
                required=disk.estimate(OUT,peers,2*GIB)
                free=shutil.disk_usage(OUT).free
                C.atomic(OUT/'guard.json',dict(state='running',pid=os.getpid(),
                    free_gib=round(free/GIB,2),checkpoint_gib=round(required/GIB,2)))
                if free+disk.reserve_bytes(OUT) <= required+3*GIB:
                    try:
                        save('low-disk')
                    except RuntimeError as exc:
                        if 'Inne polecenie sterujące' not in str(exc):
                            raise
                        time.sleep(2)
                        continue
                    return
                time.sleep(2)
        except Exception as exc:
            state=SESSION.state()
            if C.alive(state.get('coordinator')) and state['state']=='running':
                if not (OUT/'frozen.json').exists():
                    disk.freeze(OUT,state,disk.members(OUT),str(exc))
            C.atomic(OUT/'guard.json',dict(state='error',message=str(exc)))
            raise


def start(retry=False):
    OUT.mkdir(exist_ok=True)
    if (OUT/'finished.json').exists():
        status()
        return
    C.tools_available()
    C.restart_prefix()
    env=signature()
    verification=PREFIX/'verified.json'
    if not verification.exists() or C.read(verification)!=env:
        raise RuntimeError('Wykonaj ./scripts/eksperyment.sh verify po instalacji.')
    with C.control_lock(OUT/'launch'):
        state=SESSION.state()
        supervisor=OUT/'supervisor.json'
        if supervisor.exists() and C.alive(C.read(supervisor)):
            print('Nadzór już działa. Użyj status.')
            return
        if (OUT/'frozen.json').exists():
            raise RuntimeError('Najpierw zwolnij miejsce i wykonaj save-stop.')
        if state['state'] not in ('new','saved') and not C.alive(state.get('coordinator')):
            if not retry or disk.members(OUT):
                raise RuntimeError('Sprawdź log. Polecenie retry powtórzy niedokończone próby po awarii.')
            if state.get('signature') != env:
                raise RuntimeError('Kod lub dane zmieniły się od poprzedniego badania.')
            SESSION.directory.rename(OUT/f'checkpoint-failed-{time.time_ns()}')
        launch_env=os.environ.copy()
        if SESSION.state()['state']=='new':
            measured=resources()
            describe(measured)
            count=worker_count(measured)
            launch_env['NEPAL_WORKERS']=str(count)
            C.atomic(OUT/'resources.json',{**measured,'workers':count})
            print(f'Uruchamiam do {count} prób jednocześnie.')
        with (OUT/'supervisor.log').open('a') as f:
            proc=subprocess.Popen([sys.executable,'-Bu','-m','eksperyment.cli','supervise'],
                stdin=subprocess.DEVNULL,stdout=f,stderr=f,start_new_session=True,cwd=ROOT,env=launch_env)
        C.atomic(supervisor,C.identity(proc.pid))
    print('Uruchomiono nadzór w tle. Sprawdź ./scripts/eksperyment.sh status oraz wyniki/supervisor.log.')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('cmd',choices=('zasoby','start','status','save-stop','retry','worker','study','supervise','verify'))
    ap.add_argument('number',nargs='?',type=int)
    a=ap.parse_args()
    if a.cmd=='zasoby':
        describe(resources())
    elif a.cmd=='worker':
        from .run import worker
        worker(a.number)
    elif a.cmd=='study':
        from .study import study
        study()
    elif a.cmd=='verify':
        from .smoke import verify
        C.tools_available()
        C.restart_prefix()
        verify()
        from .check_restore import verify as verify_restore
        verify_restore()
        C.atomic(PREFIX/'verified.json',signature())
    elif a.cmd=='status':
        status()
    elif a.cmd=='save-stop':
        save('manual')
    elif a.cmd=='supervise':
        supervise()
    else:
        start(a.cmd=='retry')


if __name__=='__main__':
    main()

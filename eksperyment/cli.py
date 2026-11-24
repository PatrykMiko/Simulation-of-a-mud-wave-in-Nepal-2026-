"""Background controller for the 30 m Nepal experiment on Linux."""
import argparse
from contextlib import contextmanager
import fcntl
import os
import shutil
import subprocess
import sys
import time

import optuna
from nepal_final import disk
from ravaflow_tools import checkpoint as C, nepal as N
from . import autosave
from .build import PREFIX, ROOT
from .run import OUT, config
from .resources import resources, describe, worker_count

SESSION = C.Session(OUT/'checkpoint')
GIB = disk.GIB


def signature():
    paths = [*sorted((ROOT/'eksperyment').glob('*.py')), ROOT/'eksperyment/monitor.c.inc',
             ROOT/'eksperyment/dynamics.c.inc', ROOT/'configs/nepal30_objective.json',
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
    backup = autosave.latest(OUT)
    if backup:
        print('Ostatni automatyczny zapis:', autosave.age(backup), 'ukończone próby:', backup['completed'])
    if (autosave.store(OUT)/'restoring.json').exists():
        print('Odtwarzanie plików przerwane. Ponów polecenie recover.')
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
        print('Cel:',best.get('objective_score'),'Dynamika:',best.get('model_status'))
        for p in best.get('dynamics',[]):
            print(p['id'],'dojście [s]:',p['arrival_s'],'prędkość [m/s]:',p['speed_m_s'])


def save_locked(reason):
    state=SESSION.state()
    if state['state']=='saved':
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
    # During autosave only the supervisor log sees DMTCP's stop messages.
    SESSION.save_stop()
    C.atomic(OUT/'guard.json',dict(state='saved',reason=reason))


@contextmanager
def operation_lock():
    directory = OUT/'operations'
    directory.mkdir(parents=True,exist_ok=True)
    with (directory/'control.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        yield


def save(reason):
    # An already scheduled supervisor must not launch a consumed image after
    # a simultaneous manual request has confirmed that it is safe to power off.
    if reason=='manual':
        with C.control_lock(OUT/'launch'):
            C.atomic(OUT/'pause-request.json',{'requested':time.time()})
    with operation_lock():
        if (OUT/'finished.json').is_file():
            print('Seria ukończona. Wyniki zapisane.')
            return
        save_locked(reason)
        print('ZAPIS UKOŃCZONY. Możesz wyłączyć komputer.',flush=True)


def periodic_save(required):
    with C.control_lock(OUT/'operations'):
        if SESSION.state()['state']!='running' or (OUT/'finished.json').is_file():
            return False
        save_locked('automatic')
        old=SESSION.state()
        C.atomic(OUT/'guard.json',dict(state='autosave-copying'))
        try:
            info=autosave.publish(OUT)
        except (OSError, subprocess.CalledProcessError) as exc:
            C.atomic(OUT/'guard.json',dict(state='saved',reason='backup-failed',message=str(exc)))
            print('Zachowano bieżący checkpoint. Sprawdź wolne miejsce przed start.',flush=True)
            return False
        print('Automatyczna kopia potwierdzona:',autosave.age(info),flush=True)
        if (OUT/'pause-request.json').exists():
            C.atomic(OUT/'guard.json',dict(state='saved',reason='manual'))
            return False
        if shutil.disk_usage(OUT).free < required+3*GIB:
            C.atomic(OUT/'guard.json',dict(state='saved',reason='low-disk-after-backup'))
            return False
        C.atomic(OUT/'guard.json',dict(state='autosave-resuming'))
        disk.allocate_reserve(OUT,2*GIB)
        SESSION.start([],signature(),startup_timeout=3600)
        disk.prune_consumed(OUT,old)
        return True


def supervise():
    with C.control_lock(OUT/'supervision'):
        state=SESSION.state()
        try:
            if state['state'] in ('new','saved'):
                with operation_lock():
                    if (OUT/'pause-request.json').exists():
                        return
                    if shutil.disk_usage(OUT).free < 4*GIB:
                        raise RuntimeError('Zwolnij miejsce przed startem (poniżej 4 GiB).')
                    disk.allocate_reserve(OUT,2*GIB)
                    bootstrap=ROOT/'nepal_gradual/bootstrap'
                    argv=['/usr/bin/env','NEPAL_RELATIVE_SLEEP=1',f'PYTHONPATH={bootstrap}:{ROOT}',
                          sys.executable,'-B','-u','-m','eksperyment.cli','study']
                    SESSION.start(argv,signature(),startup_timeout=3600)
                    disk.prune_consumed(OUT,state)
            next_save = time.monotonic()+autosave.PERIOD
            saved_batch = len(list((OUT/'records').glob('*.json')))//5
            while True:
                state=SESSION.state()
                if state['state']=='saved':
                    return
                if state['state']=='saving':
                    save('checkpoint')
                    return
                if not C.alive(state.get('coordinator')):
                    if (OUT/'finished.json').exists():
                        disk.release_reserve(OUT)
                        autosave.remove(autosave.store(OUT))
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
                    save('low-disk')
                    return
                batch = len(list((OUT/'records').glob('*.json')))//5
                if time.monotonic() >= next_save or saved_batch < batch < 4:
                    try:
                        resumed = periodic_save(required)
                    except RuntimeError as exc:
                        if 'Inne polecenie sterujące' not in str(exc):
                            raise
                        continue
                    if not resumed:
                        return
                    saved_batch = batch
                    next_save = time.monotonic()+autosave.PERIOD
                time.sleep(2)
        except Exception as exc:
            state=SESSION.state()
            if C.alive(state.get('coordinator')) and state['state']=='running':
                if not (OUT/'frozen.json').exists():
                    disk.freeze(OUT,state,disk.members(OUT),str(exc))
            C.atomic(OUT/'guard.json',dict(state='error',message=str(exc)))
            raise


def start(retry=False, recover=False):
    OUT.mkdir(exist_ok=True)
    if (OUT/'finished.json').exists() and not (autosave.store(OUT)/'restoring.json').exists():
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
        restoring=(autosave.store(OUT)/'restoring.json').exists()
        if recover:
            if C.alive(state.get('coordinator')) or disk.members(OUT):
                raise RuntimeError('Aktywne procesy. Użyj save-stop.')
            if state['state']=='saved' and not restoring and not autosave.latest(OUT):
                print('Dostępny jest bieżący checkpoint. Wznawiam go.')
            else:
                info=autosave.restore(OUT,env)
                print('Przywrócono automatyczny zapis:',autosave.age(info))
                state=SESSION.state()
        elif restoring:
            raise RuntimeError('Niedokończone odtworzenie. Użyj recover.')
        if (OUT/'frozen.json').exists():
            raise RuntimeError('Najpierw zwolnij miejsce i wykonaj save-stop.')
        if state['state'] not in ('new','saved') and not C.alive(state.get('coordinator')):
            if not retry and autosave.latest(OUT):
                raise RuntimeError('Znaleziono automatyczny zapis. Użyj recover, aby wrócić do tego punktu.')
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
        (OUT/'pause-request.json').unlink(missing_ok=True)
        with (OUT/'supervisor.log').open('a') as f:
            proc=subprocess.Popen([sys.executable,'-Bu','-m','eksperyment.cli','supervise'],
                stdin=subprocess.DEVNULL,stdout=f,stderr=f,start_new_session=True,cwd=ROOT,env=launch_env)
        C.atomic(supervisor,C.identity(proc.pid))
    print('Uruchomiono nadzór w tle. Sprawdź ./scripts/eksperyment.sh status oraz wyniki/supervisor.log.')


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('cmd',choices=('zasoby','start','status','save-stop','retry','recover','worker','study','supervise','verify'))
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
        start(retry=a.cmd=='retry', recover=a.cmd=='recover')


if __name__=='__main__':
    main()

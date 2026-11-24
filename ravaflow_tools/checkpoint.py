"""Planned shutdown/restart of the whole Nepal computation using DMTCP 4.2.0.

Only save-and-exit checkpoints are supported: the filesystem must not advance
beyond a saved memory image. An image is consumed when restarted. This is NOT
automatic recovery from a power failure during a running computation.
"""
from __future__ import annotations

import argparse
from contextlib import contextmanager
import fcntl
import hashlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import subprocess
import sys
import time

from . import nepal

ROOT = nepal.ROOT
BIN = ROOT / '.local-tools/dmtcp-4.2.0/bin'
FILE_PLUGIN = ROOT / '.local-tools/libnepal_checkpoint_files.so'
FILE_POLICY = ROOT / 'scripts/dmtcp-file-policy.c'
DEFAULT = ROOT / 'configs/nepal_v2_checkpoint.yaml'


def read(path):
    return json.loads(Path(path).read_text())


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4 * 1024**2), b''):
            h.update(chunk)
    return h.hexdigest()


def atomic(path, value):
    nepal.write_json(path, value)
    with Path(path).open('rb') as f:
        os.fsync(f.fileno())
    fd = os.open(Path(path).parent, os.O_RDONLY | os.O_DIRECTORY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def identity(pid):
    start = nepal.proc_identity(pid)
    return {'pid': pid, 'start': start,
            'boot': Path('/proc/sys/kernel/random/boot_id').read_text().strip()}


def alive(process):
    return bool(process and process.get('start') and identity(process['pid']) == process)


def tools_available():
    for name in ('dmtcp_coordinator', 'dmtcp_launch', 'dmtcp_restart', 'mtcp_restart', 'dmtcp_command'):
        if not (BIN / name).is_file():
            raise RuntimeError('Brak DMTCP 4.2.0. Wykonaj: ./scripts/install-dmtcp.sh')
    if not FILE_PLUGIN.is_file():
        raise RuntimeError('Brak ochrony plików środowiska. Wykonaj: ./scripts/build-checkpoint-plugin.sh')


def restart_prefix():
    """Keep mtcp_restart's initial brk below the saved non-PIE Python heap.

    With ASLR, the loader can start above the saved brk. Linux then refuses
    to shrink below start_brk, leaving libc's cached break inconsistent with
    the kernel. A later numpy.loadtxt realloc can segfault. Limit the change
    to the restart process tree; never change the system-wide ASLR setting.
    Probe before consuming a checkpoint, also before launching a new study.
    """
    executable = shutil.which('setarch')
    if not executable:
        raise RuntimeError('Brak setarch (pakiet util-linux), wymaganego do odtwarzania pamięci.')
    prefix = [executable, '-R']
    try:
        result = subprocess.run([*prefix, '/bin/true'], capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise RuntimeError('Nie można sprawdzić setarch -R; nie uruchomiono odtwarzania.') from exc
    if result.returncode:
        raise RuntimeError('System nie pozwala na setarch -R, wymagane do odtworzenia sterty. '
                           f'Nie uruchomiono zadania: {result.stderr.strip()}')
    return prefix


def environment_signature():
    files = sorted((ROOT / 'ravaflow_tools').glob('*.py'))
    files += sorted(BIN.glob('dmtcp_*'))
    files.append(BIN / 'mtcp_restart')
    if setarch := shutil.which('setarch'):
        files.append(Path(setarch).resolve())
    files += sorted((BIN.parent / 'lib').rglob('*.so'))
    files += [Path(sys.executable).resolve(), Path('/usr/bin/python3').resolve()]
    files += [FILE_PLUGIN, FILE_POLICY]
    for base in (Path.home() / '.grass8/addons', Path.home() / '.local/share/grass8/addons'):
        files += sorted((base / 'bin').glob('r.avaflow*'))
        files += sorted((base / 'scripts').glob('r.avaflow*'))
    grass = shutil.which('grass')
    if grass:
        files.append(Path(grass).resolve())
    return {'files': {str(p): digest(p) for p in files if p.is_file()},
            'file_modes': {str(p): stat.S_IMODE(p.stat().st_mode) for p in files if p.is_file()},
            'python': sys.version, 'kernel': platform.release(), 'machine': platform.machine(),
            'libc': list(platform.libc_ver()), 'root': str(ROOT)}


@contextmanager
def control_lock(directory):
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / 'control.lock').open('a') as f:
        try:
            fcntl.flock(f, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            raise RuntimeError('Inne polecenie sterujące już działa; poczekaj na jego koniec.') from None
        yield


class Session:
    """External controller; stays outside the checkpointed process tree."""

    def __init__(self, directory):
        self.directory = Path(directory).resolve()
        self.path = self.directory / 'state.json'

    def state(self):
        return read(self.path) if self.path.exists() else {'state': 'new'}

    def command(self, state, command):
        if not alive(state.get('coordinator')):
            raise RuntimeError('Koordynator nie działa; nie wysłano sygnału do starego portu/PID.')
        result = subprocess.run([str(BIN / 'dmtcp_command'), '--coord-host', '127.0.0.1',
            '--coord-port', str(state['port']), command], capture_output=True, text=True,
            timeout=15)
        # v4.2.0 can return 2 after Kc: its second request sees the ongoing checkpoint.
        # Completion is checked independently using the coordinator and saved files.
        if result.returncode and not (command == '--kcheckpoint' and result.returncode == 2
                and 'checkpointed and killed' in result.stdout):
            raise RuntimeError(f'DMTCP: {result.stdout}\n{result.stderr}')
        return result.stdout

    def start(self, argv, signature, *, cwd=ROOT, startup_timeout=900):
        if startup_timeout <= 0:
            raise ValueError('startup_timeout musi być dodatni')
        tools_available()
        with control_lock(self.directory):
            old = self.state()
            if alive(old.get('coordinator')):
                raise RuntimeError('To zadanie już działa. Użyj status albo save-stop.')
            if old['state'] not in ('new', 'saved'):
                raise RuntimeError('Brak bezpiecznego punktu wznowienia. Sprawdź log; nie odtwarzam zużytego zapisu.')
            restarting = old['state'] == 'saved'
            if restarting:
                if old['signature'] != signature:
                    raise RuntimeError('Zmieniono środowisko, kod lub dane. Nie można bezpiecznie odtworzyć zapisu.')
                for item in old['files']:
                    p = self.directory / item['path']
                    if not p.is_file() or p.stat().st_size != item['size'] or digest(p) != item['sha256']:
                        raise RuntimeError(f'Niekompletny lub zmieniony checkpoint: {p}')
            restart = restart_prefix()
            generation = self.directory / f"generation_{old.get('generation', -1) + 1:04d}"
            generation.mkdir()
            port_file = generation / 'port'
            temp = self.directory / 'tmp'
            temp.mkdir(exist_ok=True)
            env = {k: v for k, v in os.environ.items() if not k.startswith('DMTCP_')}
            env.update(DMTCP_CHECKPOINT_DIR=str(generation), DMTCP_TMPDIR=str(temp),
                       TMPDIR=str(temp), DMTCP_UNIQUE_CKPT_PLUGIN='0', DMTCP_GZIP='1',
                       PYTHONFAULTHANDLER='1',
                       NEPAL_DMTCP_FILE_ROOT=str(self.directory.parent))
            with (generation / 'coordinator.log').open('a') as log:
                coord = subprocess.Popen([str(BIN / 'dmtcp_coordinator'), '--coord-port', '0',
                    '--port-file', str(port_file), '--ckptdir', str(generation), '--tmpdir', str(temp),
                    '--exit-on-last', '--interval', '0', '--stale-timeout', str(int(startup_timeout) + 60)],
                    cwd=cwd, env=env, stdin=subprocess.DEVNULL, stdout=log, stderr=log,
                    start_new_session=True)
            state = {'state': 'starting', 'generation': old.get('generation', -1) + 1,
                     'signature': signature, 'coordinator': identity(coord.pid)}
            # Consume the previous image BEFORE any restored process can run.
            atomic(self.path, state)
            for _ in range(100):
                if port_file.exists():
                    break
                if coord.poll() is not None:
                    raise RuntimeError(f'Nie wystartował koordynator: {generation / "coordinator.log"}')
                time.sleep(.1)
            if not port_file.exists():
                coord.terminate()
                raise RuntimeError('Koordynator nie zapisał portu.')
            state['port'] = int(port_file.read_text())
            common = ['--join-coordinator', '--coord-host', '127.0.0.1', '--coord-port',
                      str(state['port']), '--ckptdir', str(generation), '--tmpdir', str(temp)]
            if restarting:
                images = [str(self.directory / x['path']) for x in old['files'] if x['path'].endswith('.dmtcp')]
                command = [*restart, str(BIN / 'dmtcp_restart'), *common, *images]
            else:
                command = [str(BIN / 'dmtcp_launch'), *common,
                           '--with-plugin', str(FILE_PLUGIN), *argv]
            with (self.directory / 'job.log').open('a') as log:
                worker = subprocess.Popen(command, cwd=cwd, env=env, stdin=subprocess.DEVNULL,
                    stdout=log, stderr=log, start_new_session=True)
            state['worker'] = identity(worker.pid)
            atomic(self.path, state)
            self.wait_ready(state, worker, restarting=restarting, timeout=startup_timeout)

    def wait_ready(self, state, worker, *, restarting, timeout):
        """Large compressed images can take minutes to restore before RUNNING=yes."""
        started = time.monotonic()
        deadline = started + timeout
        next_notice = started
        latest = 'brak statusu'
        while time.monotonic() < deadline:
            returncode = worker.poll()
            if returncode is not None:
                raise RuntimeError(f'Proces zakończył się przy starcie (kod {returncode}). '
                                   f'Zobacz {self.directory / "job.log"}')
            latest = self.command(state, '--status')
            if 'RUNNING=yes' in latest and re.search(r'NUM_PEERS=[1-9]', latest):
                state['state'] = 'running'
                atomic(self.path, state)
                print('WZNOWIONO Z PAMIĘCI' if restarting else 'URUCHOMIONO Z OBSŁUGĄ CHECKPOINT', flush=True)
                print(f'Log: {self.directory / "job.log"}', flush=True)
                return
            now = time.monotonic()
            if now >= next_notice:
                peers = re.search(r'NUM_PEERS=(\d+)', latest)
                stage = 'Odtwarzanie pamięci' if restarting else 'Uruchamianie procesów'
                print(f'{stage}: {now - started:.0f} s / limit {timeout:.0f} s; '
                      f'połączone procesy: {peers.group(1) if peers else "?"}. Czekam...', flush=True)
                next_notice = now + 15
            time.sleep(min(1, max(0, deadline - time.monotonic())))
        raise RuntimeError(f'Brak potwierdzenia działania po {timeout:.0f} s. '
                           f'Ostatni status: {latest.strip()}. '
                           f'Sprawdź {self.directory / "job.log"} oraz log koordynatora. '
                           'Nie uruchamiaj drugiego zadania bez sprawdzenia statusu.')

    def save_stop(self):
        with control_lock(self.directory):
            state = self.state()
            if state['state'] == 'saved':
                print('Zapis już ukończony; zadanie jest zatrzymane.')
                return
            if state['state'] != 'saving':
                status = self.command(state, '--status')
                if 'RUNNING=yes' not in status:
                    raise RuntimeError('Procesy nie są gotowe do zapisu; ponów po zakończeniu startu.')
                # Disk estimate is conservative, but not a guarantee: file sizes can still grow.
                reserve = checkpoint_reserve_bytes(state['worker']['pid'])
                free = shutil.disk_usage(self.directory).free
                if free < reserve:
                    raise RuntimeError(f'Za mało miejsca na zapis: wolne {free / 1024**3:.1f} GiB, '
                                       f'szacunek z zapasem {reserve / 1024**3:.1f} GiB. Proces nadal działa.')
                state['state'] = 'saving'
                atomic(self.path, state)
                self.command(state, '--kcheckpoint')
            print('Zapisywanie pamięci i otwartych plików. Nie wyłączaj komputera...', flush=True)
            while alive(state['coordinator']):
                time.sleep(1)
            generation = self.directory / f"generation_{state['generation']:04d}"
            images = sorted(generation.rglob('*.dmtcp'))
            scripts = list(generation.rglob('dmtcp_restart_script*.sh'))
            if not images or not scripts or any(p.stat().st_size == 0 for p in images):
                raise RuntimeError('Brak kompletnego zapisu. Nie potwierdzono bezpiecznego wznowienia.')
            # Hash all saved files (including GRASS/open-file copies), not only memory images.
            files = [p for p in sorted(generation.rglob('*')) if p.is_file() and not p.is_symlink()
                     and p.name not in ('port', 'coordinator.log')]
            state['files'] = [{'path': str(p.relative_to(self.directory)), 'size': p.stat().st_size,
                              'sha256': digest(p)} for p in files]
            os.sync()
            state.update(state='saved', saved_at=time.time())
            atomic(self.path, state)
            print('ZAPIS UKOŃCZONY. Procesy zakończone. Możesz wyłączyć komputer.', flush=True)

    def status(self):
        state = self.state()
        print('Stan zapisu:', state['state'])
        if alive(state.get('coordinator')):
            print(self.command(state, '--status'))
        elif state['state'] == 'saved':
            print('Zatrzymane z punktem wznowienia. Polecenie start odtworzy pamięć.')
        elif state['state'] != 'new':
            print('Procesy nie działają. Sprawdź wynik zadania i job.log; brak nowego potwierdzonego checkpointu.')


def checkpoint_reserve_bytes(root_pid):
    """Sum RSS and sizes of distinct regular open files across a process tree."""
    rows = {}
    for p in Path('/proc').glob('[0-9]*'):
        try:
            raw = (p / 'stat').read_text().rsplit(')', 1)[1].split()
            rows[int(p.name)] = (int(raw[1]), int(raw[21]) * os.sysconf('SC_PAGE_SIZE'))
        except (OSError, ValueError, IndexError):
            continue
    children = {root_pid}
    while True:
        expanded = children | {pid for pid, (parent, _) in rows.items() if parent in children}
        if expanded == children:
            break
        children = expanded
    size = sum(rows.get(pid, (0, 0))[1] for pid in children)
    files = set()
    import stat
    for pid in children:
        for fd in Path(f'/proc/{pid}/fd').glob('*'):
            try:
                st = fd.stat()
                key = st.st_dev, st.st_ino
                if stat.S_ISREG(st.st_mode) and key not in files:
                    size += st.st_size
                    files.add(key)
            except OSError:
                pass
    return int(size * 1.3) + 2 * 1024**3


def discard_old(cfg):
    """Explicitly requested removal of one legacy study, only after it has stopped."""
    path = Path(cfg['study_dir'])
    if path.parent != ROOT / 'calibration_runs' or path.is_symlink():
        raise RuntimeError('Usuwanie dozwolone tylko dla jednego katalogu bezpośrednio w calibration_runs.')
    if (path / 'checkpoint').exists():
        raise RuntimeError('To badanie ma checkpointy; discard-old jest tylko dla starego zwykłego zadania.')
    if not path.exists():
        print('Katalog już nie istnieje.')
        return
    for _ in range(30):
        if not nepal.running(cfg):
            break
        time.sleep(1)
    if nepal.running(cfg):
        raise RuntimeError('Zadanie nadal działa. Najpierw wykonaj nepal.sh stop z tym samym config.')
    with nepal.job_lock(cfg):
        shutil.rmtree(path)
    print(f'Usunięto tylko: {path}')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('command', choices=['start', 'save-stop', 'status', 'discard-old', 'worker'])
    ap.add_argument('--config', default=str(DEFAULT))
    ap.add_argument('--yes', action='store_true')
    args = ap.parse_args()
    cfg = nepal.load(args.config)
    session = Session(Path(cfg['study_dir']) / 'checkpoint')
    if args.command == 'discard-old':
        if not args.yes:
            raise RuntimeError('Dodaj --yes, aby usunąć wyniki wskazanego starego badania.')
        discard_old(cfg)
    elif args.command == 'worker':
        # Use epoch time for durations across reboot; includes time spent powered off.
        # No wall-time deadline is permitted for a checkpointed calculation.
        from types import SimpleNamespace
        clock = SimpleNamespace(monotonic=time.time, perf_counter=time.time, time=time.time, sleep=time.sleep)
        nepal.time = clock
        nepal.runner.time = clock
        sys.argv = [sys.argv[0], 'run', '--config', str(Path(args.config).resolve())]
        nepal.main()
    elif args.command == 'status':
        session.status()
        result = Path(cfg['study_dir']) / 'job_result.json'
        if result.exists():
            print('Ostatni raport aplikacji (przy pauzie może pozostać running):', read(result))
        current = Path(cfg['study_dir']) / 'current.json'
        if current.exists():
            print('Przebieg:', read(current))
    elif args.command == 'save-stop':
        session.save_stop()
    else:
        if any(cfg.get(k) is not None for k in ('timeout_seconds', 'final_timeout_seconds', 'max_hours')):
            raise RuntimeError('Checkpoint wymaga timeout_seconds, final_timeout_seconds i max_hours: null.')
        tools_available()
        result = Path(cfg['study_dir']) / 'job_result.json'
        if result.exists() and read(result).get('state') == 'finished':
            print('Badanie już zakończone. Wyniki:', cfg['study_dir'])
            return
        # A live legacy job must not be duplicated by the new launcher.
        with nepal.job_lock(cfg):
            pass
        signature = {'environment': environment_signature(), 'study': nepal.fingerprint(cfg)}
        validation = ROOT / '.local-tools/checkpoint_validation.json'
        if not validation.exists() or read(validation).get('environment') != signature['environment']:
            raise RuntimeError('Najpierw sprawdź zgodność z solverem: ./scripts/nepal-checkpoint-verify.sh')
        session.start([sys.executable, '-B', '-u', '-m', 'ravaflow_tools.checkpoint', 'worker',
                       '--config', str(Path(args.config).resolve())], signature)


if __name__ == '__main__':
    main()

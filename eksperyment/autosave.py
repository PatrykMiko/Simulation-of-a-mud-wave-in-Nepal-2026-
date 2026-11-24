"""Keep one consistent offline copy of the process images AND their mutable files."""
from pathlib import Path
import os
import errno
import shutil
import subprocess
import time

from ravaflow_tools import checkpoint as C

PERIOD = 3600  # elapsed running seconds between automatic saves
# Only experiment-owned entries. Controller locks, PID files and the backup
# store itself must stay outside the restored tree.
DATA = ('work', 'plans', 'artifacts', 'records', 'protocol.json', 'resources.json',
        'best.json', 'ranking.csv', 'paczka.tmp', 'paczka.zip')


def store(root):
    return Path(root)/'.autosave'


def latest(root):
    path = store(root)/'latest.json'
    return C.read(path) if path.exists() else None


def remove(path):
    if path.is_symlink() or path.is_file():
        path.unlink()
    elif path.is_dir():
        shutil.rmtree(path)


def copy(source, dest):
    # GNU cp uses CoW where supported, otherwise independent regular files.
    # Archive mode preserves GRASS permissions, symlinks and special files.
    subprocess.run(['cp','-a','--reflink=auto','--', source, dest],check=True)


def publish(root):
    root = Path(root)
    state = C.Session(root/'checkpoint').state()
    if state['state'] != 'saved' or C.alive(state.get('coordinator')):
        raise RuntimeError('Kopia wymaga zatrzymanego i potwierdzonego checkpointu.')
    base = store(root)
    base.mkdir(exist_ok=True)
    previous = latest(root)
    slot = 1 - previous['slot'] if previous else 0
    dest = base/f'copy-{slot}'
    remove(dest)
    (dest/'checkpoint').mkdir(parents=True)
    try:
        entries = [Path(name) for name in DATA]
        entries.extend(Path('checkpoint')/name for name in
                       ('state.json', 'tmp', 'job.log', f"generation_{state['generation']:04d}"))
        entries = [path for path in entries if (root/path).exists()]
        # Budget a non-CoW copy while both the old backup and the current
        # stopped job are still intact. Leave room for the OS and metadata.
        def allocated(path):
            files = path.rglob('*') if path.is_dir() else (path,)
            return sum(p.lstat().st_blocks * 512 for p in files)
        needed = sum(allocated(root/path) for path in entries)
        if shutil.disk_usage(base).free <= needed + 2*1024**3:
            raise OSError(errno.ENOSPC,'Za mało miejsca na niezależną kopię punktu zapisu.')
        for path in entries:
            copy(root/path,dest/path)
        # Publish after all files (including un-opened maps and frames) are durable.
        os.sync()
        info = dict(slot=slot, saved_at=state['saved_at'], generation=state['generation'],
                    completed=len(list((dest/'records').glob('*.json'))))
        C.atomic(base/'latest.json',info)
    except BaseException:
        # Never remove a slot if its publication may already have succeeded.
        if (latest(root) or {}).get('slot') != slot:
            remove(dest)
        raise
    if previous:
        remove(base/f"copy-{previous['slot']}")
    return info


def restore(root, expected_signature):
    """Called with no live peers and the launch lock held. Retryable after power loss."""
    root = Path(root)
    base = store(root)
    info = latest(root)
    if not info:
        raise RuntimeError('Brak zatwierdzonego automatycznego zapisu.')
    src = base/f"copy-{info['slot']}"
    state = C.Session(src/'checkpoint').state()
    if state['state'] != 'saved' or state['signature'] != expected_signature:
        raise RuntimeError('Środowisko zapisanego przebiegu nie odpowiada bieżącemu.')
    for entry in state['files']:
        p = src/'checkpoint'/entry['path']
        if not p.is_file() or p.stat().st_size != entry['size'] or C.digest(p) != entry['sha256']:
            raise RuntimeError(f'Uszkodzony automatyczny checkpoint: {p}')
    # Keep the published copy intact until a successor has been published. If
    # interrupted half-way, start/retry cannot use the half-restored workspace.
    marker = base/'restoring.json'
    C.atomic(marker,info)
    for name in (*DATA, 'checkpoint'):
        remove(root/name)
        if (src/name).exists():
            copy(src/name, root/name)
    for name in ('guard.json','frozen.json'):
        (root/name).unlink(missing_ok=True)
    os.sync()
    marker.unlink()
    os.sync()
    return info


def age(info):
    return time.strftime('%Y-%m-%d %H:%M:%S %Z',time.localtime(info['saved_at']))

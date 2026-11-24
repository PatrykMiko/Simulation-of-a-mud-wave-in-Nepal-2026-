"""External-only disk guard; never runs inside the checkpointed solver tree."""
import os
import re
from pathlib import Path
import shutil
import signal
import stat

from ravaflow_tools import checkpoint as C

GIB = 1024**3


def listed_pids(text):
    """DMTCP reports virtual and native PIDs; /proc needs the latter."""
    return [int(pid) for pid in re.findall(r'\[\d+:(\d+)\]@', text)]


def members(root):
    """Native coordinator peer list survives restart's rewritten /proc/environ."""
    session = C.Session(Path(root) / 'checkpoint')
    state = session.state()
    coordinator = state.get('coordinator')
    pids = set()
    if C.alive(coordinator):
        try:
            pids.update(listed_pids(session.command(state, '--list')))
        except RuntimeError:
            # Concurrent save may retire the coordinator while sampling.
            if session.state().get('state') not in ('saving', 'saved') and C.alive(coordinator):
                raise
        pids.add(coordinator['pid'])
    # Include newly spawned descendants not yet registered with DMTCP.
    if C.alive(state.get('worker')):
        pids.add(state['worker']['pid'])
    marker = os.fsencode('NEPAL_DMTCP_FILE_ROOT=' + str(Path(root).resolve()))
    parents = {}
    for p in Path('/proc').glob('[0-9]*'):
        try:
            pid = int(p.name)
            raw = (p / 'stat').read_text().rsplit(')', 1)[1].split()
            parents[pid] = int(raw[1])
            if marker in (p / 'environ').read_bytes().split(b'\0'):
                pids.add(pid)
        except (OSError, ValueError, IndexError):
            continue
    while True:
        expanded = pids | {pid for pid, parent in parents.items() if parent in pids}
        if expanded == pids:
            break
        pids = expanded
    result = []
    for pid in sorted(pids):
        try:
            ident = C.identity(pid)
            if C.alive(ident):
                result.append(ident)
        except (OSError, ValueError):
            continue
    return result


def estimate(root, peers, minimum):
    """Conservative RSS + distinct open regular files; NOT a hard upper bound."""
    size = 0
    files = set()
    for peer in peers:
        if not C.alive(peer):
            continue
        pid = peer['pid']
        try:
            rss = int(Path(f'/proc/{pid}/statm').read_text().split()[1])
            size += rss * os.sysconf('SC_PAGE_SIZE')
        except OSError:
            pass
        for fd in Path(f'/proc/{pid}/fd').glob('*'):
            try:
                st = fd.stat()
                key = (st.st_dev, st.st_ino)
                if stat.S_ISREG(st.st_mode) and key not in files:
                    files.add(key)
                    size += st.st_size
            except OSError:
                continue
    return max(minimum, int(size * 1.5) + 2 * GIB)


def reserve_path(root):
    return Path(root) / '.checkpoint-reserve'


def reserve_bytes(root):
    p = reserve_path(root)
    return p.stat().st_blocks * 512 if p.exists() else 0


def allocate_reserve(root, size):
    """Real allocated blocks, not a sparse placeholder. Only called before running."""
    p = reserve_path(root)
    if p.is_symlink():
        raise RuntimeError('Plik rezerwy nie może być dowiązaniem.')
    if reserve_bytes(root) >= size:
        return
    fd = os.open(p, os.O_CREAT | os.O_RDWR, 0o600)
    try:
        os.posix_fallocate(fd, 0, size)
        os.fsync(fd)
    except BaseException:
        os.close(fd)
        p.unlink(missing_ok=True)
        raise
    else:
        os.close(fd)
    if reserve_bytes(root) < size:
        raise RuntimeError('System plików nie przydzielił wymaganej fizycznej rezerwy.')


def release_reserve(root):
    p = reserve_path(root)
    if p.is_symlink():
        raise RuntimeError('Nieprawidłowa ścieżka rezerwy.')
    p.unlink(missing_ok=True)


def tree_bytes(path):
    return sum(p.stat().st_blocks * 512 for p in Path(path).rglob('*')
               if p.is_file() and not p.is_symlink())


def pressure(free, required, margin):
    return free <= required + margin


def freeze(root, state, peers, reason):
    """Last resort: retain RAM, stop file growth. This is NOT a durable checkpoint."""
    ids = [p for p in peers if p != state.get('coordinator') and C.alive(p)]
    C.atomic(Path(root) / 'frozen.json', {'processes': ids, 'reason': reason,
        'warning': 'Tylko pauza w RAM. Nie wyłączaj komputera. Zwolnij miejsce i użyj save-stop.'})
    for peer in ids:
        if C.alive(peer):
            try:
                os.kill(peer['pid'], signal.SIGSTOP)
            except ProcessLookupError:
                pass


def thaw(root):
    p = Path(root) / 'frozen.json'
    if not p.exists():
        return
    ids = C.read(p)['processes']
    if not all(C.alive(i) for i in ids):
        raise RuntimeError('Proces z pauzy awaryjnej zniknął; nie ma potwierdzonego checkpointu.')
    for peer in ids:
        os.kill(peer['pid'], signal.SIGCONT)
    p.unlink()


def prune_consumed(root, old):
    """Only after confirmed restart; never delete the current or saved generation."""
    session = C.Session(Path(root) / 'checkpoint')
    current = session.state()
    if old.get('state') != 'saved' or current.get('state') != 'running':
        return 0
    if current['generation'] <= old['generation']:
        raise RuntimeError('Nie potwierdzono nowej generacji checkpointu.')
    directory = session.directory / f"generation_{old['generation']:04d}"
    if directory.is_symlink():
        raise RuntimeError('Nieprawidłowy katalog checkpointu.')
    # Refuse to unlink any image/file still held by the restored group.
    for peer in members(root):
        for fd in Path(f"/proc/{peer['pid']}/fd").glob('*'):
            try:
                if os.readlink(fd).startswith(str(directory) + '/'):
                    raise RuntimeError('Odtworzony proces nadal używa starego obrazu pamięci.')
            except FileNotFoundError:
                pass
    count = tree_bytes(directory) if directory.exists() else 0
    if directory.exists():
        shutil.rmtree(directory)
    return count

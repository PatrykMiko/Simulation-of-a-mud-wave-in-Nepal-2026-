"""Conservative concurrency estimate before loading a full 30 m grid."""
import os
import shutil
from pathlib import Path

from .build import ROOT

GIB = 1024**3


def resources():
    with Path('/proc/meminfo').open() as f:
        available = next(int(line.split()[1])*1024 for line in f if line.startswith('MemAvailable:'))
    disk = shutil.disk_usage(ROOT).free
    cpus = len(os.sched_getaffinity(0))
    # Allow for the manager, the incumbent, compression and fluctuations in
    # solver working sets. Values are estimates for these particular rasters.
    estimate = min(5, cpus, max(0, (available//GIB-3)//10), max(0, (disk//GIB-10)//32))
    return dict(cpu=cpus, ram_gib=round(available/GIB,1), disk_gib=round(disk/GIB,1),
                recommended=int(estimate))


def describe(measurement):
    print(f"Dostępne CPU: {measurement['cpu']}, RAM: {measurement['ram_gib']} GiB, "
          f"wolny dysk: {measurement['disk_gib']} GiB.")
    print('Orientacyjna liczba jednoczesnych prób (maks. 5):', measurement['recommended'])
    if not measurement['recommended']:
        print('Szacunek nie daje zapasu dla jednej próby. Zwolnij zasoby lub wybierz większy komputer.')


def worker_count(measurement):
    explicit = os.environ.get('NEPAL_WORKERS')
    if explicit is None:
        if measurement['recommended']==0:
            raise RuntimeError('Za mały zapas według oszacowania. Sprawdź komendę zasoby i SYMULACJA.md.')
        return measurement['recommended']
    number = int(explicit)
    if not 1 <= number <= 5:
        raise ValueError('NEPAL_WORKERS musi być od 1 do 5.')
    if number > measurement['recommended']:
        print('Wybrana ręcznie liczba przekracza oszacowanie zasobów. Obserwuj RAM i dysk.')
    return number

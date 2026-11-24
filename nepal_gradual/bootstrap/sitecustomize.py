"""Only loaded in the NEW job. Relative waits survive a monotonic-clock reset.

CPython's time.sleep uses an absolute CLOCK_MONOTONIC deadline on Linux.
A saved deadline from the previous boot can be hours ahead after restarting.
Use libc nanosleep's relative interval, retaining its remaining interval on EINTR.
No solver calculations, Python installation or system clocks are changed.
"""
import ctypes
import errno
import math
import os
import time


class Timespec(ctypes.Structure):
    _fields_ = [('tv_sec', ctypes.c_long), ('tv_nsec', ctypes.c_long)]


_libc = ctypes.CDLL(None, use_errno=True)
_nanosleep = _libc.nanosleep
_nanosleep.argtypes = [ctypes.POINTER(Timespec), ctypes.POINTER(Timespec)]
_nanosleep.restype = ctypes.c_int


def relative_sleep(seconds):
    seconds = float(seconds)
    if not math.isfinite(seconds) or seconds < 0:
        raise ValueError('sleep length must be finite and non-negative')
    whole = math.floor(seconds)
    req = Timespec(whole, min(999999999, math.ceil((seconds-whole)*1e9)))
    while True:
        remaining = Timespec()
        if _nanosleep(ctypes.byref(req), ctypes.byref(remaining)) == 0:
            return
        code = ctypes.get_errno()
        if code != errno.EINTR:
            raise OSError(code, os.strerror(code))
        req = remaining


if os.environ.get('NEPAL_RELATIVE_SLEEP') == '1':
    time.sleep = relative_sleep

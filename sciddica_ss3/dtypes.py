"""Dtype helpers for SCIDDICA-SS3."""

from __future__ import annotations

import os

import numpy as np

EPS64 = 1.0e-9
EPS32 = np.float32(1e-12)


def compute_dtype(use_float32: bool) -> np.dtype:
    return np.dtype(np.float32) if use_float32 else np.dtype(np.float64)


def compute_eps(dtype) -> float:
    dt = np.dtype(dtype)
    if dt == np.float32:
        return float(EPS32)
    return float(EPS64)


def as_compute_array(arr, dtype, copy: bool = False) -> np.ndarray:
    arr = np.asarray(arr)
    if arr.dtype == dtype and not copy:
        return arr
    return np.asarray(arr, dtype=dtype).copy() if copy else np.asarray(arr, dtype=dtype)


def resolve_numba_flag(value, env=None) -> bool:
    if value is not None:
        return bool(value)
    if env is None:
        env = os.environ
    raw = env.get("SCIDDICA_USE_NUMBA")
    if raw is None:
        return False
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def _cpp_available() -> bool:
    try:
        import sciddica_ss3_cpp  # noqa: F401
        return True
    except ImportError:
        return False


def resolve_cpp_flag(value, env=None) -> bool:
    if value is not None:
        return bool(value) and _cpp_available()
    if env is None:
        env = os.environ
    raw = env.get("SCIDDICA_USE_CPP")
    if raw is not None:
        requested = str(raw).strip().lower() in {"1", "true", "yes", "on"}
        return requested and _cpp_available()
    return False

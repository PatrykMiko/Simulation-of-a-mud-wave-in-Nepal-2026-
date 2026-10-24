#!/usr/bin/env python3
"""Grid geometry and neighbourhood operations."""

from __future__ import annotations

import numpy as np

from .constants import NEIGHBOUR_OFFSETS, SQRT3


def shifted(
    array: np.ndarray,
    dq: int,
    dr: int,
    fill_value: float = 0.0,
) -> np.ndarray:
    """Fetch neighbour values without wrapping at the grid boundary."""
    if array.ndim != 2:
        raise ValueError("shifted() expects a two-dimensional array")

    height, width = array.shape
    result = np.full_like(array, fill_value)

    r0, r1 = max(0, -dr), min(height, height - dr)
    q0, q1 = max(0, -dq), min(width, width - dq)

    if r0 < r1 and q0 < q1:
        result[r0:r1, q0:q1] = array[
            r0 + dr : r1 + dr,
            q0 + dq : q1 + dq,
        ]

    return result


def neighbour_stack(
    array: np.ndarray,
    fill_value: float = 0.0,
) -> np.ndarray:
    """Return the six neighbouring values with shape (6, H, W)."""
    return np.stack(
        [
            shifted(array, dq, dr, fill_value)
            for dq, dr in NEIGHBOUR_OFFSETS
        ],
        axis=0,
    )


def cell_centres(
    height: int,
    width: int,
    side_length: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Cartesian cell centres for axial hex coordinates, in metres."""
    if height < 1 or width < 1:
        raise ValueError("Grid dimensions must be positive")

    r, q = np.mgrid[0:height, 0:width].astype(float)

    x = side_length * SQRT3 * (q + 0.5 * r)
    y = side_length * 1.5 * r

    return x, y


def hex_vertices(
    x: np.ndarray,
    y: np.ndarray,
    side_length: float,
) -> np.ndarray:
    """Cartesian vertices of each hexagon, shape (H*W, 6, 2)."""
    angles = np.deg2rad(30.0 + 60.0 * np.arange(6, dtype=float))

    vx = x.ravel()[:, None] + side_length * np.cos(angles)[None, :]
    vy = y.ravel()[:, None] + side_length * np.sin(angles)[None, :]

    return np.stack((vx, vy), axis=-1)

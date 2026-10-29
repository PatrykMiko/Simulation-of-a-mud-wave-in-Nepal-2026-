"""Resample north-up rasters onto the hexagonal cell centres of the SS3 model.

Hexagon (r, q) is stored at array[r, q] and has its centre at
    x = pc * sqrt(3) * (q + r / 2) - xoff,    y = 1.5 * pc * r
measured from the raster's top-left corner (x to the east, y downwards).
"""
from __future__ import annotations

import numpy as np
from scipy.ndimage import map_coordinates

SQRT3 = float(np.sqrt(3.0))


def hex_shape(ncols: int, nrows: int, cellsize: float, pc: float) -> tuple[int, int]:
    """(H, W) of an axial array that covers the whole raster."""
    H = int(nrows * cellsize // (1.5 * pc)) + 1
    W = int(ncols * cellsize // (SQRT3 * pc)) + H // 2 + 2
    return H, W


def to_hex(arr: np.ndarray, meta: dict, pc: float, order: int = 1, chunk_rows: int = 256):
    """Sample `arr` at the centre of every hexagon.

    order=1 is bilinear (use for elevation, thickness); order=0 is nearest
    (use for masks). Returns (values, valid, xoff), where `valid` marks hexagons
    whose centre lies inside the raster, and `xoff` is the x shift [m] that makes
    the leftmost hexagons start left of the raster.
    """
    arr = np.asarray(arr, dtype=np.float64)
    nrows, ncols = arr.shape
    cs = float(meta["cellsize"])
    H, W = hex_shape(ncols, nrows, cs, pc)
    xoff = pc * SQRT3 * H / 2.0

    values = np.empty((H, W), dtype=np.float64)
    valid = np.empty((H, W), dtype=bool)
    q = np.arange(W, dtype=np.float64)[None, :]

    for r0 in range(0, H, chunk_rows):
        r1 = min(r0 + chunk_rows, H)
        r = np.arange(r0, r1, dtype=np.float64)[:, None]

        x = pc * SQRT3 * (q + 0.5 * r) - xoff
        y = 1.5 * pc * r

        col = x / cs - 0.5               # pixel centres sit at (j + 0.5) * cs
        row = np.broadcast_to(y / cs - 0.5, col.shape)

        valid[r0:r1] = (
            (col >= -0.5) & (col <= ncols - 0.5)
            & (row >= -0.5) & (row <= nrows - 0.5)
        )
        values[r0:r1] = map_coordinates(arr, [row, col], order=order, mode="nearest")

    return values, valid, xoff
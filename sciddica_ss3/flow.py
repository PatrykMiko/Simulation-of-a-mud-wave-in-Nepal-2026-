#!/usr/bin/env python3
"""Minimization Algorithm of the Differences and momentum corrections."""

from __future__ import annotations

import numpy as np

from .constants import EPS
from .state import ModelParameters


def minimization_algorithm(
    movable_mass: np.ndarray,
    central_unmovable_height: np.ndarray,
    neighbour_heights: np.ndarray,
    allowed: np.ndarray,
) -> np.ndarray:
    """Vectorised Minimization Algorithm of the Differences."""
    mass = np.asarray(movable_mass, dtype=float)
    centre = np.asarray(central_unmovable_height, dtype=float)
    neighbours = np.asarray(neighbour_heights, dtype=float)
    admissible = np.asarray(allowed, dtype=bool)

    if mass.ndim != 2 or centre.shape != mass.shape:
        raise ValueError(
            "movable_mass and central_unmovable_height must share shape (H, W)"
        )

    if neighbours.shape != (6, *mass.shape) or admissible.shape != neighbours.shape:
        raise ValueError("neighbour_heights and allowed must have shape (6, H, W)")

    if np.any(mass < -EPS):
        raise ValueError("movable_mass cannot be negative")

    if not np.all(np.isfinite(mass)) or not np.all(np.isfinite(centre)):
        raise ValueError("MA inputs must be finite")

    mass = np.maximum(mass, 0.0)

    active = admissible & (mass[None, :, :] > EPS)

    if not np.any(active):
        return np.zeros_like(neighbours)

    levels = np.concatenate(
        (
            centre[None, :, :],
            np.where(active, neighbours, np.inf),
        ),
        axis=0,
    )

    sorted_levels = np.sort(levels, axis=0)
    finite_levels = np.isfinite(sorted_levels)

    cumulative = np.cumsum(
        np.where(finite_levels, sorted_levels, 0.0),
        axis=0,
    )

    counts = np.arange(1, 8, dtype=float)[:, None, None]
    averages = (mass[None, :, :] + cumulative) / counts

    filled = finite_levels & (averages > sorted_levels)
    k_filled = filled.sum(axis=0)

    gather_index = np.maximum(k_filled - 1, 0)[None, :, :]

    fill_level = np.take_along_axis(
        averages,
        gather_index,
        axis=0,
    )[0]

    outflow = np.maximum(
        fill_level[None, :, :] - levels[1:],
        0.0,
    )

    outflow = np.where(active, outflow, 0.0)

    total = outflow.sum(axis=0)

    scale = np.minimum(
        1.0,
        np.divide(
            mass,
            total,
            out=np.ones_like(mass),
            where=total > EPS,
        ),
    )

    return np.maximum(outflow * scale[None, :, :], 0.0)


def momentum_elevation_correction(
    thickness: np.ndarray,
    momentum_x: np.ndarray,
    momentum_y: np.ndarray,
    params: ModelParameters,
) -> np.ndarray:
    """Compute SS3 fictitious elevation corrections for six neighbours."""
    th = np.asarray(thickness, dtype=float)
    px = np.asarray(momentum_x, dtype=float)
    py = np.asarray(momentum_y, dtype=float)

    if px.shape != th.shape or py.shape != th.shape:
        raise ValueError("thickness and momentum arrays must share shape")

    magnitude = np.hypot(px, py) / (params.rho * params.cell_area)

    moving = (magnitude > EPS) & (th > EPS)

    angle_sector = (np.arctan2(py, px) % (2.0 * np.pi)) / (np.pi / 3.0)

    lower = np.floor(angle_sector).astype(int) % 6
    fraction = angle_sector - np.floor(angle_sector)

    shares = np.zeros((6, *th.shape), dtype=float)

    for direction in range(6):
        shares[direction] = magnitude * (
            np.where(lower == direction, 1.0 - fraction, 0.0)
            + np.where(
                (lower + 1) % 6 == direction,
                fraction,
                0.0,
            )
        )

    shares = np.where(moving[None, :, :], shares, 0.0)

    return params.alpha * (np.roll(shares, 3, axis=0) - shares)

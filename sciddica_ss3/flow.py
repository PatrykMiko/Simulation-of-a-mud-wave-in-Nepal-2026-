#!/usr/bin/env python3
"""Minimization Algorithm of the Differences and momentum corrections."""

from __future__ import annotations

import numpy as np

from .constants import EPS, GRAVITY, DIRECTION_X, DIRECTION_Y
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
    kinetic_head: np.ndarray,
    barycentre_x: np.ndarray,
    barycentre_y: np.ndarray,
    params: ModelParameters,
) -> np.ndarray:
    """Fictitious-slope correction w(i) for the six neighbours (SS3, phase 1).

    - Magnitude from kinetic head: p/(rho*A) = th * sqrt(2 g kh)   (paper's formula),
      with speed capped at one cell spacing per step (the CA cannot move farther).
    - Direction from the momentum vector, split between the two faces whose
      directions (seen from the barycentre) bracket it.
    - Result is bounded: |w| <= 2 * th * alpha_scale-independent cap (see below).
    """
    th = np.asarray(thickness, dtype=float)
    px = np.asarray(momentum_x, dtype=float)
    py = np.asarray(momentum_y, dtype=float)
    kh = np.maximum(np.asarray(kinetic_head, dtype=float), 0.0)
    bx = np.asarray(barycentre_x, dtype=float)
    by = np.asarray(barycentre_y, dtype=float)

    if not (px.shape == py.shape == kh.shape == bx.shape == by.shape == th.shape):
        raise ValueError("all (H, W) inputs must share shape")

    spacing = params.centre_spacing
    apothem = spacing / 2.0

    # 1. magnitude from energy, bounded by one cell per step
    speed = np.minimum(np.sqrt(2.0 * GRAVITY * kh), spacing / params.pt)
    magnitude = th * speed                       # = p / (rho * A)

    moving = (magnitude > EPS) & (th > EPS) & (np.hypot(px, py) > EPS)
    theta = np.arctan2(py, px)

    # 2. keep the barycentre well inside the polygon of face midpoints
    #    (its inradius is 0.866*apothem) so the six face directions surround it
    #    in index order and every angular gap stays < 180 degrees.
    r_max = 0.6 * apothem
    r = np.hypot(bx, by)
    scale = np.minimum(1.0, np.divide(r_max, r, out=np.ones_like(r), where=r > EPS))
    bx = bx * scale
    by = by * scale

    # 3. angles of the six face midpoints seen from the barycentre.
    #    Relative to face 0 they increase monotonically with the face index.
    face_x = apothem * DIRECTION_X[:, None, None] - bx[None]
    face_y = apothem * DIRECTION_Y[:, None, None] - by[None]
    face_angles = np.arctan2(face_y, face_x)
    rel = np.mod(face_angles - face_angles[0], 2.0 * np.pi)       # rel[0] == 0
    t_rel = np.mod(theta - face_angles[0], 2.0 * np.pi)

    lower = np.clip(np.sum(t_rel[None] >= rel, axis=0) - 1, 0, 5)
    upper = (lower + 1) % 6

    rel_lower = np.take_along_axis(rel, lower[None], axis=0)[0]
    rel_upper = np.take_along_axis(rel, upper[None], axis=0)[0]
    rel_upper = np.where(lower == 5, 2.0 * np.pi, rel_upper)

    fraction = np.clip(
        (t_rel - rel_lower) / np.maximum(rel_upper - rel_lower, EPS), 0.0, 1.0
    )

    # 4. shares on two ADJACENT faces (never opposite, so no cancellation)
    k = np.arange(6)[:, None, None]
    shares = magnitude[None] * (
        np.where(k == lower[None], 1.0 - fraction[None], 0.0)
        + np.where(k == upper[None], fraction[None], 0.0)
    )
    shares = np.where(moving[None], shares, 0.0)

    # 5. lower the privileged faces, raise their opposites, and bound the result:
    #    a drop of 2*th already moves the entire mass one cell.
    w = params.alpha * (np.roll(shares, 3, axis=0) - shares)
    cap = 2.0 * th[None]
    return np.clip(w, -cap, cap)
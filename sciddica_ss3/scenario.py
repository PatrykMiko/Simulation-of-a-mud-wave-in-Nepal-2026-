#!/usr/bin/env python3
"""Synthetic scenario and validation metrics."""

from __future__ import annotations

import numpy as np

from .automaton import SciddicaSS3
from .geometry import cell_centres
from .state import ModelParameters
from .constants import SQRT3


def make_albano_like(
    n: int = 80,
    params: ModelParameters | None = None,
) -> tuple[SciddicaSS3, np.ndarray, np.ndarray]:
    """Create a synthetic channel and lake."""
    if n < 3:
        raise ValueError("n must be at least 3")

    p = params or ModelParameters()

    x, y = cell_centres(n, n, p.pc)
    max_y = float(y.max())

    if max_y <= 0:
        raise ValueError("Scenario requires at least two rows")

    base = 300.0 - 250.0 * (y / max_y)

    row_mid = p.pc * SQRT3 * (n / 2.0 + np.arange(n)[:, None] / 2.0)

    channel_centre = (
        row_mid
        + 0.08 * n * p.pc * np.sin(2.0 * np.pi * 1.5 * y / max_y)
    )

    channel = np.exp(-((x - channel_centre) / (4.0 * p.pc)) ** 2)

    base = np.maximum(base, 60.0)

    altitude = base - 8.0 * channel
    soil = 1.5 * channel

    distance_squared = (x - channel_centre) ** 2 + (y - 0.08 * max_y) ** 2

    debris = 4.0 * np.exp(-distance_squared / (5.0 * p.pc) ** 2)
    debris[debris < 0.05] = 0.0

    altitude = altitude - debris

    water_level = 300.0 - 250.0 * 0.80

    model = SciddicaSS3(
        altitude,
        debris,
        soil,
        p,
        water_level=water_level,
    )

    return model, x, y


def jaccard_fitness(
    real_mask: np.ndarray,
    simulated_mask: np.ndarray,
) -> float:
    """SS3 overlap fitness: |R ∩ S| / |R ∪ S|."""
    real = np.asarray(real_mask, dtype=bool)
    simulated = np.asarray(simulated_mask, dtype=bool)

    if real.shape != simulated.shape:
        raise ValueError("real_mask and simulated_mask must have the same shape")

    union = np.logical_or(real, simulated).sum()

    if union == 0:
        return 1.0

    return float(np.logical_and(real, simulated).sum() / union)


def ss2_overlap_fitness(
    real_mask: np.ndarray,
    simulated_mask: np.ndarray,
) -> float:
    """SS2 paper's overlap fitness: sqrt(|R ∩ S| / |R ∪ S|)."""
    return float(np.sqrt(jaccard_fitness(real_mask, simulated_mask)))

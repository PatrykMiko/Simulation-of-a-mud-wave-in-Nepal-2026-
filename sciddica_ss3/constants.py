#!/usr/bin/env python3
"""Core constants for SCIDDICA-SS3."""

import numpy as np

GRAVITY = 9.81
EPS = 1.0e-9
MIN_THICKNESS = 1.0e-6
SQRT3 = float(np.sqrt(3.0))

# Axial coordinates (q, r) are stored as array[r, q].
# The six neighbours are separated by 60 degrees.
NEIGHBOUR_OFFSETS: tuple[tuple[int, int], ...] = (
    (1, 0),
    (0, 1),
    (-1, 1),
    (-1, 0),
    (0, -1),
    (1, -1),
)

DIRECTION_ANGLES = np.deg2rad(60.0 * np.arange(6, dtype=float))
DIRECTION_X = np.cos(DIRECTION_ANGLES)
DIRECTION_Y = np.sin(DIRECTION_ANGLES)

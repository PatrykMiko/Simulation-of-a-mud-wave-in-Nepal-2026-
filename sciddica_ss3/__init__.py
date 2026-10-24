#!/usr/bin/env python3
"""SCIDDICA-SS3 package."""

from .constants import (
    GRAVITY,
    EPS,
    MIN_THICKNESS,
    SQRT3,
    NEIGHBOUR_OFFSETS,
    DIRECTION_ANGLES,
    DIRECTION_X,
    DIRECTION_Y,
)
from .geometry import shifted, neighbour_stack, cell_centres, hex_vertices
from .state import MediumParameters, ModelParameters, FlowState, OutflowPhase, ShiftResult
from .flow import minimization_algorithm, momentum_elevation_correction
from .automaton import SciddicaSS3
from .scenario import make_albano_like, jaccard_fitness
from .visualize import run_headless, animate

__all__ = [
    "GRAVITY", "EPS", "MIN_THICKNESS", "SQRT3",
    "NEIGHBOUR_OFFSETS", "DIRECTION_ANGLES", "DIRECTION_X", "DIRECTION_Y",
    "shifted", "neighbour_stack", "cell_centres", "hex_vertices",
    "MediumParameters", "ModelParameters", "FlowState", "OutflowPhase", "ShiftResult",
    "minimization_algorithm", "momentum_elevation_correction",
    "SciddicaSS3", "make_albano_like", "jaccard_fitness",
    "run_headless", "animate",
]
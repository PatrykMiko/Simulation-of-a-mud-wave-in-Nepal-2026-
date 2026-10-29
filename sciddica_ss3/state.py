#!/usr/bin/env python3
"""Parameters and model state structures."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .constants import GRAVITY, SQRT3


@dataclass(frozen=True)
class MediumParameters:
    """Parameters that may differ between subaerial and subaqueous flow."""

    adh: float = 0.02
    pf: float = 0.05
    pmt: float = 2.0
    per: float = 0.05
    ptd: float = 0.50
    ped: float = 0.30
    shift_friction: float = 0.05

    def __post_init__(self) -> None:
        non_negative = (
            self.adh,
            self.pf,
            self.pmt,
            self.per,
            self.shift_friction,
        )

        if not all(np.isfinite(value) and value >= 0 for value in non_negative):
            raise ValueError(
                "adh, pf, pmt, per and shift_friction must be "
                "finite and non-negative"
            )

        if not all(np.isfinite(value) and 0 <= value <= 1 for value in (self.ptd, self.ped)):
            raise ValueError("ptd and ped must be finite and in [0, 1]")


@dataclass(frozen=True)
class ModelParameters:
    """Global SCIDDICA-SS3 parameters."""

    pc: float = 12.0
    pt: float = 2.0
    rho: float = 2000.0

    air: MediumParameters = field(default_factory=MediumParameters)

    water: MediumParameters = field(
        default_factory=lambda: MediumParameters(
            adh=0.02,
            pf=0.05,
            pmt=2.0,
            per=0.02,
            ptd=0.70,
            ped=0.30,
            shift_friction=0.10,
        )
    )

    pml: float = 0.25
    pwr: float = 0.10
    rho_water: float = 1000.0
    alpha_scale: float = 0.25

    def __post_init__(self) -> None:
        if not all(np.isfinite(value) and value > 0 for value in (self.pc, self.pt, self.rho)):
            raise ValueError("pc, pt and rho must be finite and positive")

        if not np.isfinite(self.pml) or not 0 <= self.pml <= 1:
            raise ValueError("pml must be finite and in [0, 1]")

        if not np.isfinite(self.pwr) or self.pwr < 0:
            raise ValueError("pwr must be finite and non-negative")

        if not np.isfinite(self.rho_water) or self.rho_water <= 0:
            raise ValueError("rho_water must be finite and positive")

        if not np.isfinite(self.alpha_scale) or self.alpha_scale < 0:
            raise ValueError("alpha_scale must be finite and non-negative")

    @property
    def cell_area(self) -> float:
        """Area of a regular hexagon with side length pc [m²]."""
        return 1.5 * SQRT3 * self.pc**2

    @property
    def centre_spacing(self) -> float:
        """Distance between neighbouring hexagonal cell centres [m]."""
        return SQRT3 * self.pc

    @property
    def submerged_gravity(self) -> float:
        """Effective submerged gravity derived from a buoyancy relation."""
        density_contrast = max(self.rho - self.rho_water, 0.0)
        return GRAVITY * density_contrast / self.rho

    @property
    def alpha(self) -> float:
        """Empirical coefficient for the SS3 fictitious-slope correction."""
        return self.alpha_scale * 2.0 * self.pt / self.centre_spacing


@dataclass
class FlowState:
    """Cellular-automaton substates."""

    altitude: np.ndarray
    thickness: np.ndarray
    soil_depth: np.ndarray
    kinetic_head: np.ndarray
    momentum_x: np.ndarray
    momentum_y: np.ndarray
    barycentre_x: np.ndarray
    barycentre_y: np.ndarray


@dataclass
class OutflowPhase:
    """Mass, energy and centroid data for one outflow phase."""

    requested: np.ndarray
    external: np.ndarray
    received: np.ndarray
    kinetic_head: np.ndarray
    shifted_kinetic_head: np.ndarray
    fraction: np.ndarray
    travelled: np.ndarray
    internal: np.ndarray
    internal_barycentre_x: np.ndarray
    internal_barycentre_y: np.ndarray
    external_barycentre_x: np.ndarray
    external_barycentre_y: np.ndarray


@dataclass
class ShiftResult:
    """Geometry and kinematics from shifting six directional cylinders."""

    external_fraction: np.ndarray
    travelled: np.ndarray
    kinetic_head: np.ndarray
    internal_barycentre_x: np.ndarray
    internal_barycentre_y: np.ndarray
    external_barycentre_x: np.ndarray
    external_barycentre_y: np.ndarray

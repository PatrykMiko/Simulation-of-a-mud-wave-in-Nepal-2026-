#!/usr/bin/env python3
"""Shifting, cylinder splitting, and medium utilities."""

from __future__ import annotations

import numpy as np

from .constants import EPS, GRAVITY, DIRECTION_X, DIRECTION_Y
from .state import ModelParameters, ShiftResult


def _medium_arrays(
    submerged: np.ndarray,
    params: ModelParameters,
) -> tuple[np.ndarray, ...]:
    """Expand air/water parameter selection to cellwise arrays."""
    names = (
        "adh",
        "pf",
        "pmt",
        "per",
        "ptd",
        "ped",
        "shift_friction",
    )

    values: list[np.ndarray] = []

    for name in names:
        values.append(
            np.where(
                submerged,
                getattr(params.water, name),
                getattr(params.air, name),
            )
        )

    return tuple(values)


def _circle_segment(
    z: np.ndarray,
    radius: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return circular-cap fraction and centroid offsets."""
    zc = np.clip(z, -radius, radius)

    root = np.sqrt(np.maximum(radius * radius - zc * zc, 0.0))

    cap_area = radius * radius * np.arccos(zc / radius) - zc * root

    total_area = np.pi * radius * radius
    internal_area = np.maximum(total_area - cap_area, 0.0)

    fraction = np.clip(cap_area / total_area, 0.0, 1.0)

    first_moment = (2.0 / 3.0) * np.maximum(radius * radius - zc * zc, 0.0) ** 1.5

    external_centroid = np.divide(
        first_moment,
        cap_area,
        out=np.zeros_like(cap_area),
        where=cap_area > EPS,
    )

    internal_centroid = -np.divide(
        first_moment,
        internal_area,
        out=np.zeros_like(internal_area),
        where=internal_area > EPS,
    )

    external_centroid = np.where(fraction >= 1.0 - 1e-12, 0.0, external_centroid)
    internal_centroid = np.where(fraction <= 1e-12, 0.0, internal_centroid)

    return fraction, external_centroid, internal_centroid


def _shift_fraction(
    kinetic_head: np.ndarray,
    source_surface: np.ndarray,
    neighbour_surface: np.ndarray,
    barycentre_x: np.ndarray,
    barycentre_y: np.ndarray,
    shift_friction: np.ndarray,
    submerged: np.ndarray,
    params: ModelParameters,
) -> ShiftResult:
    """Apply inherited SS2 shifting and cylinder-splitting geometry."""
    kh = np.maximum(np.asarray(kinetic_head, dtype=float), 0.0)

    bx = np.asarray(barycentre_x, dtype=float)
    by = np.asarray(barycentre_y, dtype=float)
    fc = np.asarray(shift_friction, dtype=float)
    sub = np.asarray(submerged, dtype=bool)

    expected_shape = (6, *kh.shape)

    if neighbour_surface.shape != expected_shape:
        raise ValueError("neighbour_surface must have shape (6, H, W)")

    if (
        bx.shape != kh.shape
        or by.shape != kh.shape
        or fc.shape != kh.shape
        or sub.shape != kh.shape
    ):
        raise ValueError("shift inputs must share shape (H, W)")

    fractions = np.zeros(expected_shape, dtype=float)
    travelled = np.zeros(expected_shape, dtype=float)

    heads = np.broadcast_to(kh, expected_shape).copy()
    internal_bx = np.broadcast_to(bx, expected_shape).copy()
    internal_by = np.broadcast_to(by, expected_shape).copy()

    external_bx = np.zeros(expected_shape, dtype=float)
    external_by = np.zeros(expected_shape, dtype=float)

    spacing = params.centre_spacing
    apothem = spacing / 2.0

    initial_speed = np.sqrt(2.0 * GRAVITY * kh)
    submerged_g = params.submerged_gravity

    for direction in range(6):
        ux = DIRECTION_X[direction]
        uy = DIRECTION_Y[direction]

        projected_barycentre = bx * ux + by * uy

        dx = spacing * ux - bx
        dy = spacing * uy - by

        horizontal_distance = np.maximum(np.hypot(dx, dy), EPS)

        geometric_drop = source_surface - neighbour_surface[direction]

        theta = np.arctan2(geometric_drop, horizontal_distance)

        can_move = (geometric_drop + kh) >= fc * horizontal_distance

        acceleration_gravity = np.where(sub, submerged_g, GRAVITY)

        acceleration = acceleration_gravity * (np.sin(theta) - fc * np.cos(theta))

        speed0 = initial_speed
        dt = params.pt

        # Subaerial constant-acceleration shift
        air_speed = np.maximum(speed0 + acceleration * dt, 0.0)
        air_shift = speed0 * dt + 0.5 * acceleration * dt * dt

        air_stops = (acceleration < 0.0) & ((speed0 + acceleration * dt) < 0.0)
        stopping_time = np.divide(
            -speed0,
            acceleration,
            out=np.zeros_like(speed0),
            where=acceleration < -EPS,
        )
        stopping_time = np.clip(stopping_time, 0.0, dt)
        stopping_shift = speed0 * stopping_time + 0.5 * acceleration * stopping_time**2
        air_shift = np.where(air_stops, stopping_shift, air_shift)
        air_speed = np.where(air_stops, 0.0, air_speed)

        # Subaqueous motion under linear drag
        wr = params.pwr

        if wr > EPS:
            decay = np.exp(-wr * dt)
            water_speed = speed0 * decay + (acceleration / wr) * (1.0 - decay)
            water_shift = (
                speed0 * (1.0 - decay) / wr
                + (acceleration / wr) * (dt - (1.0 - decay) / wr)
            )
            stop_denominator = wr * speed0 - acceleration
            stop_ratio = np.divide(
                -acceleration,
                stop_denominator,
                out=np.ones_like(speed0),
                where=stop_denominator > EPS,
            )
            water_stopping_time = -np.log(np.clip(stop_ratio, EPS, 1.0)) / wr
            water_stops = (
                (acceleration < 0.0)
                & (speed0 >= 0.0)
                & (water_stopping_time < dt)
            )
            stop_decay = np.exp(-wr * water_stopping_time)
            water_stopping_shift = (
                speed0 * (1.0 - stop_decay) / wr
                + (acceleration / wr)
                * (water_stopping_time - (1.0 - stop_decay) / wr)
            )
            water_shift = np.where(water_stops, water_stopping_shift, water_shift)
            water_speed = np.where(water_stops, 0.0, water_speed)
        else:
            water_speed = np.maximum(speed0 + acceleration * dt, 0.0)
            water_shift = speed0 * dt + 0.5 * acceleration * dt * dt
            water_stops = (acceleration < 0.0) & ((speed0 + acceleration * dt) < 0.0)
            water_stopping_time = np.divide(
                -speed0,
                acceleration,
                out=np.zeros_like(speed0),
                where=acceleration < -EPS,
            )
            water_stopping_time = np.clip(water_stopping_time, 0.0, dt)
            water_stopping_shift = (
                speed0 * water_stopping_time + 0.5 * acceleration * water_stopping_time**2
            )
            water_shift = np.where(water_stops, water_stopping_shift, water_shift)
            water_speed = np.where(water_stops, 0.0, water_speed)

        new_speed = np.where(sub, np.maximum(water_speed, 0.0), air_speed)
        raw_shift = np.where(sub, water_shift, air_shift)
        raw_shift = np.maximum(raw_shift, 0.0)
        raw_shift = np.where(can_move, raw_shift, 0.0)

        remaining_to_centre = np.maximum(spacing - projected_barycentre, 0.0)
        shift = np.minimum(raw_shift, remaining_to_centre)
        shift = np.where(can_move, shift, 0.0)

        signed_chord = apothem - projected_barycentre - shift
        frac, cap_centroid, inner_centroid = _circle_segment(signed_chord, apothem)
        frac = np.where(can_move, frac, 0.0)
        cap_centroid = np.where(can_move, cap_centroid, 0.0)
        inner_centroid = np.where(can_move, inner_centroid, 0.0)

        perpendicular_x = bx - projected_barycentre * ux
        perpendicular_y = by - projected_barycentre * uy

        inner_along = projected_barycentre + shift + inner_centroid
        outer_along = projected_barycentre + shift + cap_centroid - spacing

        fractions[direction] = frac
        travelled[direction] = shift
        heads[direction] = np.where(
            can_move,
            new_speed * new_speed / (2.0 * GRAVITY),
            kh,
        )
        internal_bx[direction] = np.where(
            can_move,
            perpendicular_x + inner_along * ux,
            bx,
        )
        internal_by[direction] = np.where(
            can_move,
            perpendicular_y + inner_along * uy,
            by,
        )
        external_bx[direction] = perpendicular_x + outer_along * ux
        external_by[direction] = perpendicular_y + outer_along * uy

    return ShiftResult(
        external_fraction=fractions,
        travelled=travelled,
        kinetic_head=np.maximum(heads, 0.0),
        internal_barycentre_x=internal_bx,
        internal_barycentre_y=internal_by,
        external_barycentre_x=external_bx,
        external_barycentre_y=external_by,
    )


def _external_flow_after_water_loss(
    external: np.ndarray,
    submerged: np.ndarray,
    neighbour_submerged: np.ndarray,
    params: ModelParameters,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply air-to-water matter loss, returning received and lost thickness."""
    enters_water = (~submerged)[None, :, :] & neighbour_submerged

    loss_fraction = np.where(enters_water, params.pml, 0.0)

    received = external * (1.0 - loss_fraction)
    lost_depth = external - received

    return received, lost_depth


def _internal_barycentre(
    mass_before: np.ndarray,
    barycentre_x: np.ndarray,
    barycentre_y: np.ndarray,
    requested: np.ndarray,
    external: np.ndarray,
    internal_barycentre_x: np.ndarray,
    internal_barycentre_y: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Mass-weight stationary and shifted internal-cylinder centroids."""
    stationary = np.maximum(mass_before - requested.sum(axis=0), 0.0)

    internal = np.maximum(requested - external, 0.0)

    sum_x = stationary * barycentre_x
    sum_y = stationary * barycentre_y

    for direction in range(6):
        sum_x += internal[direction] * internal_barycentre_x[direction]
        sum_y += internal[direction] * internal_barycentre_y[direction]

    remaining = np.maximum(mass_before - external.sum(axis=0), 0.0)

    bx = np.divide(
        sum_x,
        remaining,
        out=barycentre_x.copy(),
        where=remaining > EPS,
    )
    by = np.divide(
        sum_y,
        remaining,
        out=barycentre_y.copy(),
        where=remaining > EPS,
    )

    return bx, by, remaining

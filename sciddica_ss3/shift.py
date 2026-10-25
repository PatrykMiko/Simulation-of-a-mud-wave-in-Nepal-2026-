#!/usr/bin/env python3
"""Shifting, cylinder splitting, and medium utilities."""

from __future__ import annotations

import numpy as np
from numba import njit, prange

from .constants import EPS, GRAVITY, DIRECTION_X, DIRECTION_Y
from .state import ModelParameters, ShiftResult


def _medium_arrays(
    submerged: np.ndarray,
    params: ModelParameters,
) -> tuple[np.ndarray, ...]:
    """Expand air/water parameter selection to cellwise arrays."""
    names = ("adh", "pf", "pmt", "per", "ptd", "ped", "shift_friction")

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


@njit(cache=True)
def _circle_segment_scalar(z, r, eps):
    zc = min(max(z, -r), r)
    rad = max(r * r - zc * zc, 0.0)
    root = np.sqrt(rad)
    cap = r * r * np.arccos(zc / r) - zc * root
    total = np.pi * r * r
    internal = max(total - cap, 0.0)
    frac = min(max(cap / total, 0.0), 1.0)
    fm = (2.0 / 3.0) * rad ** 1.5
    ec = fm / cap if cap > eps else 0.0
    ic = -fm / internal if internal > eps else 0.0
    if frac >= 1.0 - 1e-12:
        ec = 0.0
    if frac <= 1e-12:
        ic = 0.0
    return frac, ec, ic


@njit(cache=True, parallel=True, error_model="numpy")
def _shift_kernel(
    kh_in, src, nb, bx, by, fc, sub, active,
    spacing, dt, g, g_sub, wr, ux, uy, eps,
    fractions, travelled, heads, ibx, iby, ebx, eby,
):
    N = kh_in.size
    apothem = spacing / 2.0
    for i in prange(N):
        kh = max(kh_in[i], 0.0)

        # inert defaults, identical to the dense version for cells that cannot move
        for d in range(6):
            fractions[d, i] = 0.0
            travelled[d, i] = 0.0
            heads[d, i] = kh
            ibx[d, i] = bx[i]
            iby[d, i] = by[i]
            ebx[d, i] = 0.0
            eby[d, i] = 0.0

        if not active[i]:
            continue

        speed0 = np.sqrt(2.0 * g * kh)
        is_sub = sub[i]

        for d in range(6):
            pb = bx[i] * ux[d] + by[i] * uy[d]
            dx = spacing * ux[d] - bx[i]
            dy = spacing * uy[d] - by[i]
            hdist = max(np.sqrt(dx * dx + dy * dy), eps)

            drop = src[i] - nb[d, i]
            theta = np.arctan2(drop, hdist)
            can_move = (drop + kh) >= fc[i] * hdist

            gacc = g_sub if is_sub else g
            acc = gacc * (np.sin(theta) - fc[i] * np.cos(theta))

            if is_sub and wr > eps:
                decay = np.exp(-wr * dt)
                sp = speed0 * decay + (acc / wr) * (1.0 - decay)
                sh = (speed0 * (1.0 - decay) / wr
                      + (acc / wr) * (dt - (1.0 - decay) / wr))
                den = wr * speed0 - acc
                ratio = -acc / den if den > eps else 1.0
                ratio = min(max(ratio, eps), 1.0)
                ts = -np.log(ratio) / wr
                if acc < 0.0 and speed0 >= 0.0 and ts < dt:
                    sd = np.exp(-wr * ts)
                    sh = (speed0 * (1.0 - sd) / wr
                          + (acc / wr) * (ts - (1.0 - sd) / wr))
                    sp = 0.0
                new_speed = max(sp, 0.0)
                raw = sh
            else:  # air, or water without drag
                sp = max(speed0 + acc * dt, 0.0)
                sh = speed0 * dt + 0.5 * acc * dt * dt
                if acc < 0.0 and (speed0 + acc * dt) < 0.0:
                    ts = -speed0 / acc if acc < -eps else 0.0
                    ts = min(max(ts, 0.0), dt)
                    sh = speed0 * ts + 0.5 * acc * ts * ts
                    sp = 0.0
                new_speed = sp
                raw = sh

            raw = max(raw, 0.0)
            if not can_move:
                raw = 0.0

            remaining = max(spacing - pb, 0.0)
            shift = min(raw, remaining)
            if not can_move:
                shift = 0.0

            frac, cc, ic = _circle_segment_scalar(apothem - pb - shift, apothem, eps)
            if not can_move:
                frac = 0.0
                cc = 0.0
                ic = 0.0

            perp_x = bx[i] - pb * ux[d]
            perp_y = by[i] - pb * uy[d]
            inner_along = pb + shift + ic
            outer_along = pb + shift + cc - spacing

            fractions[d, i] = frac
            travelled[d, i] = shift
            if can_move:
                heads[d, i] = new_speed * new_speed / (2.0 * g)
                ibx[d, i] = perp_x + inner_along * ux[d]
                iby[d, i] = perp_y + inner_along * uy[d]
            ebx[d, i] = perp_x + outer_along * ux[d]
            eby[d, i] = perp_y + outer_along * uy[d]


def _shift_fraction(
    kinetic_head: np.ndarray,
    source_surface: np.ndarray,
    neighbour_surface: np.ndarray,
    barycentre_x: np.ndarray,
    barycentre_y: np.ndarray,
    shift_friction: np.ndarray,
    submerged: np.ndarray,
    params: ModelParameters,
    active: np.ndarray | None = None,
) -> ShiftResult:
    """Shift all cells (or only those in `active`) with the numba kernel.

    Inactive cells carry no outflow, so they receive inert defaults; every
    later use multiplies those values by zero mass.
    """
    shape = np.shape(kinetic_head)

    def flat(a):
        return np.ascontiguousarray(a, dtype=np.float64).ravel()

    kh = flat(kinetic_head)
    n = kh.size

    if active is None:
        act = np.ones(n, dtype=np.bool_)
    else:
        act = np.ascontiguousarray(active, dtype=np.bool_).ravel()

    out = [np.empty((6, n)) for _ in range(7)]

    _shift_kernel(
        kh,
        flat(source_surface),
        np.ascontiguousarray(neighbour_surface, dtype=np.float64).reshape(6, n),
        flat(barycentre_x),
        flat(barycentre_y),
        flat(shift_friction),
        np.ascontiguousarray(submerged, dtype=np.bool_).ravel(),
        act,
        params.centre_spacing,
        params.pt,
        GRAVITY,
        params.submerged_gravity,
        params.pwr,
        np.ascontiguousarray(DIRECTION_X, dtype=np.float64),
        np.ascontiguousarray(DIRECTION_Y, dtype=np.float64),
        EPS,
        *out,
    )

    def grid(a):
        return a.reshape(6, *shape)

    fractions, travelled, heads, ibx, iby, ebx, eby = (grid(a) for a in out)

    return ShiftResult(
        external_fraction=fractions,
        travelled=travelled,
        kinetic_head=np.maximum(heads, 0.0),
        internal_barycentre_x=ibx,
        internal_barycentre_y=iby,
        external_barycentre_x=ebx,
        external_barycentre_y=eby,
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
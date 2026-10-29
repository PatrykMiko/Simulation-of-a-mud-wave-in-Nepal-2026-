#!/usr/bin/env python3
"""Core SCIDDICA-SS3 cellular automaton."""

from __future__ import annotations

import numpy as np

from .constants import EPS, GRAVITY, MIN_THICKNESS, DIRECTION_X, DIRECTION_Y, NEIGHBOUR_OFFSETS
from .geometry import neighbour_stack, shifted
from .flow import minimization_algorithm, momentum_elevation_correction
from .shift import _shift_fraction, _external_flow_after_water_loss, _internal_barycentre, _medium_arrays
from .state import FlowState, ModelParameters, OutflowPhase


class SciddicaSS3:
    """Synchronous SCIDDICA-SS3 cellular automaton on a hexagonal grid."""

    def __init__(
        self,
        altitude: np.ndarray,
        thickness: np.ndarray,
        soil_depth: np.ndarray,
        params: ModelParameters | None = None,
        water_level: float | None = None,
        use_float32: bool = False,
    ) -> None:
        self.params = params or ModelParameters()
        self.use_float32 = bool(use_float32)
        self._dtype = np.float32 if self.use_float32 else float

        if water_level is not None and not np.isfinite(water_level):
            raise ValueError("water_level must be finite or None")

        self.water_level = water_level

        a = np.asarray(altitude, dtype=self._dtype)
        th = np.asarray(thickness, dtype=self._dtype)
        soil = np.asarray(soil_depth, dtype=self._dtype)

        if a.ndim != 2 or th.shape != a.shape or soil.shape != a.shape:
            raise ValueError(
                "altitude, thickness and soil_depth must share a 2-D shape"
            )

        if a.size == 0 or not np.all(np.isfinite(a)):
            raise ValueError("altitude must be a non-empty finite array")

        if not np.all(np.isfinite(th)) or not np.all(np.isfinite(soil)):
            raise ValueError("thickness and soil_depth must be finite")

        if np.any(th < 0) or np.any(soil < 0):
            raise ValueError("thickness and soil_depth cannot be negative")

        zeros = np.zeros_like(a)

        self.state = FlowState(
            altitude=a.copy(),
            thickness=th.copy(),
            soil_depth=soil.copy(),
            kinetic_head=zeros.copy(),
            momentum_x=zeros.copy(),
            momentum_y=zeros.copy(),
            barycentre_x=zeros.copy(),
            barycentre_y=zeros.copy(),
        )

        self.height, self.width = a.shape
        self.time = 0.0
        self.mass_lost = 0.0
        self.steps_completed = 0

        self._inside = neighbour_stack(np.ones_like(a), fill_value=0.0) > 0.5
        # Pre-allocate reusable arrays to reduce allocations
        self._zeros_hw = np.zeros((self.height, self.width), dtype=self._dtype)
        self._zeros_6hw = np.zeros((6, self.height, self.width), dtype=self._dtype)
        self._ones_hw = np.ones((self.height, self.width), dtype=self._dtype)
        self._inf_hw = np.full((self.height, self.width), np.inf, dtype=self._dtype)

    @property
    def altitude(self) -> np.ndarray:
        return self.state.altitude

    @property
    def thickness(self) -> np.ndarray:
        return self.state.thickness

    @property
    def soil_depth(self) -> np.ndarray:
        return self.state.soil_depth

    @property
    def kinetic_head(self) -> np.ndarray:
        return self.state.kinetic_head

    @property
    def momentum_x(self) -> np.ndarray:
        return self.state.momentum_x

    @property
    def momentum_y(self) -> np.ndarray:
        return self.state.momentum_y

    def submerged(self, altitude: np.ndarray | None = None) -> np.ndarray:
        """Return the mask of cells below the configured water level."""
        a = self.altitude if altitude is None else altitude

        if self.water_level is None:
            return np.zeros_like(a, dtype=bool)

        return a < self.water_level

    def total_volume(self) -> float:
        """Current debris volume [m³], excluding uneroded soil."""
        return float(self.thickness.sum() * self.params.cell_area)

    def kinetic_energy(self) -> float:
        """Kinetic energy implied by Qkh [J]."""
        return float(
            self.params.rho
            * GRAVITY
            * self.params.cell_area
            * np.sum(self.thickness * self.kinetic_head)
        )

    def total_energy(self) -> float:
        """Reconstruct total debris energy using the paper's Eq. (1)."""
        cell_energy = (
            self.thickness * (0.5 * self.thickness + self.altitude + self.kinetic_head)
        )

        return float(
            self.params.rho * GRAVITY * self.params.cell_area * cell_energy.sum()
        )

    def _mobilize_soil(self) -> tuple[np.ndarray, ...]:
        """Apply progressive soil mobilization and the SS2 erosion head loss."""
        st = self.state

        old_thickness = st.thickness
        submerged = self.submerged(st.altitude)

        (
            adh,
            pf,
            pmt,
            per,
            ptd,
            ped,
            friction,
        ) = _medium_arrays(submerged, self.params)

        above_threshold = np.maximum(st.kinetic_head - pmt, 0.0)

        erosion = np.minimum(per * above_threshold, st.soil_depth)
        erosion = np.where(old_thickness > EPS, erosion, 0.0)
        erosion = np.maximum(erosion, 0.0)
        active_erosion = erosion > EPS

        head_after_erosion_loss = np.maximum(st.kinetic_head - ped * above_threshold, 0.0)

        dilution = np.divide(
            old_thickness,
            old_thickness + erosion,
            out=np.ones_like(old_thickness),
            where=(old_thickness + erosion) > EPS,
        )

        new_head = np.where(
            active_erosion,
            head_after_erosion_loss * dilution,
            st.kinetic_head,
        )

        mass_ratio = np.divide(
            old_thickness + erosion,
            old_thickness,
            out=np.ones_like(old_thickness),
            where=old_thickness > EPS,
        )

        velocity_scale = np.sqrt(
            np.divide(
                new_head,
                st.kinetic_head,
                out=np.ones_like(new_head),
                where=st.kinetic_head > EPS,
            )
        )

        momentum_scale = mass_ratio * velocity_scale

        st.altitude = st.altitude - erosion
        st.soil_depth = np.maximum(st.soil_depth - erosion, 0.0)
        st.thickness = st.thickness + erosion
        st.kinetic_head = new_head

        st.momentum_x = st.momentum_x * momentum_scale
        st.momentum_y = st.momentum_y * momentum_scale

        submerged = self.submerged(st.altitude)

        (
            adh,
            pf,
            pmt,
            per,
            ptd,
            ped,
            friction,
        ) = _medium_arrays(submerged, self.params)

        return submerged, adh, pf, ptd, friction

    def _make_phase(
        self,
        requested: np.ndarray,
        kinetic_head: np.ndarray,
        source_surface_for_shift: np.ndarray,
        neighbour_surface_for_shift: np.ndarray,
        submerged: np.ndarray,
        neighbour_submerged: np.ndarray,
        friction: np.ndarray,
        barycentre_x: np.ndarray | None = None,
        barycentre_y: np.ndarray | None = None,
    ) -> OutflowPhase:
        """Apply SS2 shifting and the SS3 air–water matter-loss rule."""
        st = self.state

        shifted_result = _shift_fraction(
            kinetic_head=kinetic_head,
            source_surface=source_surface_for_shift,
            neighbour_surface=neighbour_surface_for_shift,
            barycentre_x=st.barycentre_x if barycentre_x is None else barycentre_x,
            barycentre_y=st.barycentre_y if barycentre_y is None else barycentre_y,
            shift_friction=friction,
            submerged=submerged,
            params=self.params,
        )

        external = requested * shifted_result.external_fraction

        received, _lost_depth = _external_flow_after_water_loss(
            external,
            submerged,
            neighbour_submerged,
            self.params,
        )

        (
            _adh,
            _pf,
            _pmt,
            _per,
            ptd,
            _ped,
            _fc,
        ) = _medium_arrays(submerged, self.params)

        external_head = shifted_result.kinetic_head * (1.0 - ptd)[None, :, :]
        external_head = np.where(external > EPS, external_head, 0.0)

        internal = np.maximum(requested - external, 0.0)

        return OutflowPhase(
            requested=requested,
            external=external,
            received=received,
            kinetic_head=external_head,
            shifted_kinetic_head=shifted_result.kinetic_head,
            fraction=shifted_result.external_fraction,
            travelled=shifted_result.travelled,
            internal=internal,
            internal_barycentre_x=shifted_result.internal_barycentre_x,
            internal_barycentre_y=shifted_result.internal_barycentre_y,
            external_barycentre_x=shifted_result.external_barycentre_x,
            external_barycentre_y=shifted_result.external_barycentre_y,
        )

    def _inflows(self, phase1: OutflowPhase, phase2: OutflowPhase) -> tuple[np.ndarray, ...]:
        """Compose incoming mass, kinetic-energy proxy, momentum and centroids."""
        st = self.state
        shape = st.altitude.shape

        incoming_mass = np.zeros(shape, dtype=float)
        incoming_head_mass = np.zeros(shape, dtype=float)
        incoming_px = np.zeros(shape, dtype=float)
        incoming_py = np.zeros(shape, dtype=float)

        incoming_bx_mass = np.zeros(shape, dtype=float)
        incoming_by_mass = np.zeros(shape, dtype=float)

        for direction, (dq, dr) in enumerate(NEIGHBOUR_OFFSETS):
            mass1 = shifted(phase1.received[direction], -dq, -dr, 0.0)
            mass2 = shifted(phase2.received[direction], -dq, -dr, 0.0)

            head_mass1 = shifted(
                phase1.received[direction] * phase1.kinetic_head[direction],
                -dq,
                -dr,
                0.0,
            )
            head_mass2 = shifted(
                phase2.received[direction] * phase2.kinetic_head[direction],
                -dq,
                -dr,
                0.0,
            )

            mass_k = mass1 + mass2
            head_mass_k = head_mass1 + head_mass2

            incoming_mass += mass_k
            incoming_head_mass += head_mass_k

            mean_head = np.divide(
                head_mass_k,
                mass_k,
                out=np.zeros_like(mass_k),
                where=mass_k > EPS,
            )

            speed_k = np.sqrt(2.0 * GRAVITY * mean_head)
            momentum_k = (
                self.params.rho * self.params.cell_area * mass_k * speed_k
            )

            incoming_px += momentum_k * DIRECTION_X[direction]
            incoming_py += momentum_k * DIRECTION_Y[direction]

            source_bx_mass = (
                phase1.received[direction] * phase1.external_barycentre_x[direction]
                + phase2.received[direction] * phase2.external_barycentre_x[direction]
            )
            source_by_mass = (
                phase1.received[direction] * phase1.external_barycentre_y[direction]
                + phase2.received[direction] * phase2.external_barycentre_y[direction]
            )

            incoming_bx_mass += shifted(source_bx_mass, -dq, -dr, 0.0)
            incoming_by_mass += shifted(source_by_mass, -dq, -dr, 0.0)

        return (
            incoming_mass,
            incoming_head_mass,
            incoming_px,
            incoming_py,
            incoming_bx_mass,
            incoming_by_mass,
        )

    def step(self) -> None:
        """Advance one step, computing only on a window around the debris."""
        full_state = self.state
        full_inside = self._inside

        occupied = full_state.thickness > 0.0
        rows = np.flatnonzero(occupied.any(axis=1))

        if rows.size == 0:
            self._step_full()          # nothing moving: keeps time/counters consistent
            return

        cols = np.flatnonzero(occupied.any(axis=0))
        margin = 3                     # >= 2 is required; 3 leaves headroom
        r0 = max(int(rows[0]) - margin, 0)
        r1 = min(int(rows[-1]) + margin + 1, self.height)
        q0 = max(int(cols[0]) - margin, 0)
        q1 = min(int(cols[-1]) + margin + 1, self.width)

        # Window covers most of the grid: slicing would only add copies.
        if (r1 - r0) * (q1 - q0) > 0.6 * self.height * self.width:
            self._step_full()
            return

        window = (slice(r0, r1), slice(q0, q1))

        self.state = FlowState(
            **{name: value[window].copy() for name, value in vars(full_state).items()}
        )
        self._inside = full_inside[:, window[0], window[1]]

        try:
            self._step_full()
            for name, value in vars(self.state).items():
                getattr(full_state, name)[window] = value
        finally:
            self.state = full_state
            self._inside = full_inside

    def _step_full(self) -> None:
        """Advance the grid by one synchronous cellular-automaton step."""
        st = self.state
        p = self.params

        submerged, adh, pf, ptd, friction = self._mobilize_soil()

        a = st.altitude
        th = st.thickness
        kh = st.kinetic_head

        # Active set: cells with thickness above threshold or in neighborhood of active cells
        active_mask = th > (MIN_THICKNESS * 10.0)
        if np.any(active_mask):
            # Expand to include neighbors (1-ring)
            active_nb = neighbour_stack(active_mask.astype(float), fill_value=0.0)
            active_region = active_mask | (np.any(active_nb > 0.5, axis=0))
        else:
            active_region = np.zeros_like(th, dtype=bool)

        sub_nb = neighbour_stack(submerged.astype(float), fill_value=0.0) > 0.5
        a_nb = neighbour_stack(a, fill_value=np.inf)
        th_nb = neighbour_stack(th, fill_value=0.0)

        bed_surface = a + th
        neighbour_bed_surface = a_nb + th_nb

        # Phase 1
        correction = momentum_elevation_correction(
            thickness=th,
            momentum_x=st.momentum_x,
            momentum_y=st.momentum_y,
            kinetic_head=kh,
            barycentre_x=st.barycentre_x,
            barycentre_y=st.barycentre_y,
            params=p,
        )
        phase1_neighbour_heights = neighbour_bed_surface + correction

        allowed1 = self._inside & ((bed_surface[None, :, :] - phase1_neighbour_heights) > pf[None, :, :])
        if np.any(active_region):
            allowed1 = allowed1 & active_region[None, :, :]

        requested1 = minimization_algorithm(
            movable_mass=np.maximum(th - adh, 0.0),
            central_unmovable_height=a + adh + kh,
            neighbour_heights=phase1_neighbour_heights,
            allowed=allowed1,
        )

        phase1 = self._make_phase(
            requested=requested1,
            kinetic_head=kh,
            source_surface_for_shift=bed_surface,
            neighbour_surface_for_shift=neighbour_bed_surface,
            submerged=submerged,
            neighbour_submerged=sub_nb,
            friction=friction,
        )

        th_after1 = np.maximum(th - phase1.external.sum(axis=0), 0.0)

        bx1, by1, _ = _internal_barycentre(
            mass_before=th,
            barycentre_x=st.barycentre_x,
            barycentre_y=st.barycentre_y,
            requested=phase1.requested,
            external=phase1.external,
            internal_barycentre_x=phase1.internal_barycentre_x,
            internal_barycentre_y=phase1.internal_barycentre_y,
        )

        head1 = self._head_after_internal_flow(
            mass_before=th,
            mass_after=th_after1,
            kinetic_head=kh,
            phase=phase1,
            ptd=ptd,
            apply_turbulence=False,
        )

        px1, py1 = self._momentum_after_internal_flow(
            mass_before=th,
            momentum_x=st.momentum_x,
            momentum_y=st.momentum_y,
            phase=phase1,
        )

        # Phase 2: reconstruct the intermediate surface after Phase 1.
        phase1_inflow = np.zeros_like(th)

        for direction, (dq, dr) in enumerate(NEIGHBOUR_OFFSETS):
            phase1_inflow += shifted(
                phase1.received[direction],
                -dq,
                -dr,
                0.0,
            )

        surface_after1 = a + th_after1 + phase1_inflow

        neighbour_surface_after1 = neighbour_stack(
            surface_after1,
            fill_value=np.inf,
        )

        allowed2 = self._inside & ((surface_after1[None, :, :] - neighbour_surface_after1) > pf[None, :, :])
        if np.any(active_region):
            allowed2 = allowed2 & active_region[None, :, :]

        residual_mass_after1 = np.maximum(th - phase1.requested.sum(axis=0), 0.0)

        requested2 = minimization_algorithm(
            movable_mass=np.maximum(residual_mass_after1 - adh, 0.0),
            central_unmovable_height=a + adh + head1,
            neighbour_heights=neighbour_surface_after1,
            allowed=allowed2,
        )

        phase2 = self._make_phase(
            requested=requested2,
            kinetic_head=head1,
            source_surface_for_shift=surface_after1,
            neighbour_surface_for_shift=neighbour_surface_after1,
            submerged=submerged,
            neighbour_submerged=sub_nb,
            friction=friction,
            barycentre_x=bx1,
            barycentre_y=by1,
        )

        retained = np.maximum(th_after1 - phase2.external.sum(axis=0), 0.0)

        retained_head = self._head_after_internal_flow(
            mass_before=th_after1,
            mass_after=retained,
            kinetic_head=head1,
            phase=phase2,
            ptd=ptd,
            apply_turbulence=True,
        )

        px_retained, py_retained = self._momentum_after_internal_flow(
            mass_before=th_after1,
            momentum_x=px1,
            momentum_y=py1,
            phase=phase2,
            keep_stationary=False,
        )

        turbulence_momentum_scale = np.sqrt(np.maximum(1.0 - ptd, 0.0))
        px_retained *= turbulence_momentum_scale
        py_retained *= turbulence_momentum_scale

        (
            incoming_mass,
            incoming_head_mass,
            incoming_px,
            incoming_py,
            incoming_barycentre_x_mass,
            incoming_barycentre_y_mass,
        ) = self._inflows(phase1, phase2)

        new_thickness = retained + incoming_mass

        new_head = np.divide(
            retained * retained_head + incoming_head_mass,
            new_thickness,
            out=np.zeros_like(new_thickness),
            where=new_thickness > EPS,
        )
        new_head = np.maximum(new_head, 0.0)

        new_px = px_retained + incoming_px
        new_py = py_retained + incoming_py

        bx2, by2, retained_barycentre_mass = _internal_barycentre(
            mass_before=th_after1,
            barycentre_x=bx1,
            barycentre_y=by1,
            requested=phase2.requested,
            external=phase2.external,
            internal_barycentre_x=phase2.internal_barycentre_x,
            internal_barycentre_y=phase2.internal_barycentre_y,
        )

        bary_x_sum = retained_barycentre_mass * bx2 + incoming_barycentre_x_mass
        bary_y_sum = retained_barycentre_mass * by2 + incoming_barycentre_y_mass

        new_bx = np.divide(
            bary_x_sum,
            new_thickness,
            out=np.zeros_like(new_thickness),
            where=new_thickness > EPS,
        )
        new_by = np.divide(
            bary_y_sum,
            new_thickness,
            out=np.zeros_like(new_thickness),
            where=new_thickness > EPS,
        )

        lost_depth = (phase1.external - phase1.received).sum() + (phase2.external - phase2.received).sum()
        self.mass_lost += float(lost_depth * p.cell_area)

        dead = new_thickness < MIN_THICKNESS

        st.thickness = np.where(dead, 0.0, np.maximum(new_thickness, 0.0))
        st.kinetic_head = np.where(dead, 0.0, new_head)
        st.momentum_x = np.where(dead, 0.0, new_px)
        st.momentum_y = np.where(dead, 0.0, new_py)
        st.barycentre_x = np.where(dead, 0.0, new_bx)
        st.barycentre_y = np.where(dead, 0.0, new_by)

        self.time += p.pt
        self.steps_completed += 1

        self._assert_finite_state()

    def _head_after_internal_flow(
        self,
        mass_before: np.ndarray,
        mass_after: np.ndarray,
        kinetic_head: np.ndarray,
        phase: OutflowPhase,
        ptd: np.ndarray,
        apply_turbulence: bool,
    ) -> np.ndarray:
        """Energy-weight stationary and internal shifted material."""
        stationary = np.maximum(mass_before - phase.requested.sum(axis=0), 0.0)

        energy_proxy = stationary * np.maximum(kinetic_head, 0.0)

        for direction in range(6):
            energy_proxy += phase.internal[direction] * phase.shifted_kinetic_head[direction]

        result = np.divide(
            energy_proxy,
            mass_after,
            out=np.maximum(kinetic_head, 0.0).copy(),
            where=mass_after > EPS,
        )

        if apply_turbulence:
            result *= 1.0 - ptd

        return np.maximum(result, 0.0)

    def _momentum_after_internal_flow(
        self,
        mass_before: np.ndarray,
        momentum_x: np.ndarray,
        momentum_y: np.ndarray,
        phase: OutflowPhase,
        keep_stationary: bool = True,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Compose retained momentum with directional internal-flow momentum.

        keep_stationary=False models the paper's phase 2: debris that cannot
        overcome obstacles hits and bounces, so its momentum disappears.
        """
        stationary = np.maximum(mass_before - phase.requested.sum(axis=0), 0.0)

        stationary_fraction = np.divide(
            stationary,
            mass_before,
            out=np.zeros_like(stationary),
            where=mass_before > EPS,
        )

        if keep_stationary:
            px = momentum_x * stationary_fraction
            py = momentum_y * stationary_fraction
        else:
            px = np.zeros_like(momentum_x)
            py = np.zeros_like(momentum_y)

        for direction in range(6):
            speed = np.sqrt(
                2.0 * GRAVITY
                * np.maximum(phase.shifted_kinetic_head[direction], 0.0)
            )

            internal_momentum = (
                self.params.rho
                * self.params.cell_area
                * phase.internal[direction]
                * speed
            )

            # Do not create directional momentum when the cylinder
            # did not actually shift in this direction.
            moved = phase.travelled[direction] > EPS
            internal_momentum = np.where(
                moved,
                internal_momentum,
                0.0,
            )

            px += internal_momentum * DIRECTION_X[direction]
            py += internal_momentum * DIRECTION_Y[direction]

        return px, py
    
    def run(self, steps: int, callback=None, progress: bool = True) -> None:
        """Run a number of transitions, optionally calling callback(model, i).

        When ``progress`` is true and ``tqdm`` is available, a progress bar is
        shown.  Pass ``progress=False`` to silence it (e.g. in tests).
        """
        if steps < 0:
            raise ValueError("steps must be non-negative")

        iterator = range(steps)
        if progress:
            try:
                from tqdm import tqdm

                iterator = tqdm(range(steps), desc="steps", unit="step")
            except ImportError:
                pass

        for i in iterator:
            self.step()
            if callback is not None:
                callback(self, i)

    def _assert_finite_state(self) -> None:
        """Fail early if state substates contain invalid numerical values."""
        for name, value in vars(self.state).items():
            if not np.all(np.isfinite(value)):
                raise FloatingPointError(f"Non-finite values in state substate {name}")

        if np.any(self.state.thickness < -EPS) or np.any(self.state.soil_depth < -EPS):
            raise FloatingPointError("Negative thickness or soil depth detected")

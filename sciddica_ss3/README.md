# SCIDDICA-SS3 Technical Guide

## 1. Data conversion workflow

### Prerequisites
- Python 3.14+ with project dependencies (`numpy`, `scipy`, `matplotlib`, `rasterio` optional for GeoTIFF).
- Ravaflow input directory containing `elevation`, `release_height`, `reference_extent`, optional `entrainment_height`.

### Steps
1. Validate inputs:
   ```bash
   .venv/bin/python -m ravaflow_tools.check_inputs data/<name> --volume 2e6 --cellsize 30 --epsg 32645
   ```
2. Generate synthetic data (optional):
   ```bash
   .venv/bin/python -m ravaflow_tools.synthetic data/synthetic --format asc
   ```
3. Convert:
   ```bash
   .venv/bin/python -m ravaflow_tools.to_ss3 data/<name> --out data/<name>_ss3 --water-level 300 --volume 2e6 --pc 30
   ```
   Outputs: `ss3_arrays.npz` (altitude, thickness, soil_depth) and `ss3_meta.json`.

## 2. Running the model

### CLI
```bash
.venv/bin/python -m sciddica_ss3.cli <source_dir_or_npz> [--headless] [--steps N] [--report-every K] [--save file.gif]
```

### Parameters
- `source`: converted data directory or `.npz` file
- `--headless`: run without GUI
- `--steps`: number of steps (default 300)
- `--report-every`: print interval (default 50)
- `--save`: save animation to GIF

### Save final mask
```bash
.venv/bin/python -m sciddica_ss3.cli data/synthetic_ss3 --headless --steps 10 --mask-out data/mask.npy
```

# Short model description model description

## 1. Overview

- The terrain is a grid of **hexagonal cells**, each with **six neighbours**.
- The third dimension is a cell property: debris is stored as a **thickness** in metres, so a 2D grid describes a 3D flow.
- All cells are updated **synchronously**. One step represents `pt` seconds.
- The boundary is **fixed (closed)**: nothing flows out of the domain and nothing wraps around.
- The model is **deterministic**. The same input always gives the same output.

## 2. Grid and neighbourhood

Cells use axial coordinates `(q, r)` and are stored in arrays as `array[r, q]` (rows are `r`, columns are `q`). Hexagons are pointy-top, with side length `pc`.

## 3. One step of the transition function

### 3.1 Mobilization (erosion)

Where the kinetic head exceeds the threshold `pmt`, the flow erodes soil in proportion to the excess, limited by the soil depth left:

```
erosion = min(per * max(kinetic_head - pmt, 0), soil_depth)
```

The eroded soil is added to the debris thickness and subtracted from the bed altitude and the soil depth. The new mass dilutes the kinetic head, and a further fraction `ped` of the excess head is lost. Momentum is rescaled for the mass and speed change.

### 3.2 Momentum correction (fictitious slope)

This is what makes the flow inertial. The Minimization Algorithm knows nothing about velocity, so the heights of the neighbours are altered before it runs:

- The momentum direction is split between the two faces whose directions bracket it.
- Privileged faces (where momentum points) get their height **lowered**.
- The opposite faces get their height **raised** by the same amount.
- The size of the correction is proportional to `alpha * thickness * speed`, where `alpha = alpha_scale * 2 * pt / centre_spacing`. A drop of about twice the thickness is enough to move the whole mass one cell.

### 3.3 Phase 1: inertial outflows

The **Minimization Algorithm of the Differences (MA)** distributes the movable debris `thickness - adh` between the cell and its eligible neighbours so that their surface levels equalise.

- The central cell counts at height `altitude + adh + kinetic_head`.
- Neighbours count at `altitude + thickness + correction`.
- A neighbour is eligible only if the cell's surface is higher than the neighbour by more than the friction threshold `pf`, and the neighbour is inside the domain.
- Implementation: heights are sorted per cell, and the largest set of cells whose average level lies above its highest member is filled up to that level (vectorized with `sort` and `cumsum`).

### 3.4 Phase 2: non-inertial outflows

The debris that was not distributed in phase 1 is distributed again on the real terrain, without the momentum correction. Neighbours are counted with the debris they received in phase 1. This part of the flow hits and bounces off obstacles, so its **momentum is lost**, while some kinetic energy remains.

### 3.5 Shift

Each requested outflow is turned into actual motion:

- The geometric height drop to the neighbour sets the slope. A flow can move if the drop plus its kinetic head pays for friction over the distance.
- **Air:** constant-acceleration motion, `acceleration = g * (sin(theta) - friction * cos(theta))`, stopping if the speed reaches zero within the step.
- **Water:** the same with the buoyancy-reduced gravity and linear drag `pwr`.
- The travelled distance is limited to the distance to the neighbour's centre.
- The distance decides what fraction of the debris, modelled as a cylinder, crosses the cell face (the **external** part) and what stays inside (the **internal** part). Circular-segment geometry also gives the new barycentres of both parts.

### 3.6 Air-to-water interface

Debris moving from a dry cell into a submerged cell loses a fraction `pml` of its matter. The lost volume is accumulated in `mass_lost`.

### 3.7 Composition

Each cell combines what it keeps with what its neighbours send it:

- **Thickness:** retained debris plus all inflows.
- **Kinetic head:** mass-weighted average of retained and incoming heads. The outgoing heads have already been reduced by the turbulence loss `ptd` (the retained mass in phase 2 too).
- **Momentum:** retained momentum (zero for the phase-2 part) plus inflows, each along the direction of its flow.
- **Barycentre:** mass-weighted average of retained and incoming barycentres.
- Cells with thickness below `MIN_THICKNESS` are cleared.

### 3.8 Conservation

Debris, erodible soil and the matter lost entering water together stay constant, apart from the removal of cells thinner than `MIN_THICKNESS`. Erosion only moves material from the soil to the debris, and it only ever lowers the bed.

## 4. Cell state (substates)

Every cell holds the following values. All arrays have shape `(H, W)`.

| Substate | Field | Unit | Meaning | Initial value |
|---|---|:-:|---|---|
| `Qa` | `altitude` | m | Altitude of the bed above sea level (bedrock plus soil cover). Falls as soil is eroded. | Terrain, minus the release thickness in the source area |
| `Qth` | `thickness` | m | Thickness of the flowing debris | Release thickness in the source area, zero elsewhere |
| `Qd` | `soil_depth` | m | Erodible soil cover left | Given soil depth |
| `Qkh` | `kinetic_head` | m | Kinetic head `v^2 / 2g`. A measure of the debris speed. | 0 |
| `Qpx` | `momentum_x` | kg m/s | x component of the momentum (used for its direction) | 0 |
| `Qpy` | `momentum_y` | kg m/s | y component of the momentum | 0 |
| `Qx` | `barycentre_x` | m | x position of the debris barycentre, relative to the cell centre | 0 |
| `Qy` | `barycentre_y` | m | y position of the debris barycentre, relative to the cell centre | 0 |

Derived per cell, not stored:

| Quantity | Formula |
|---|---|
| Surface height | `altitude + thickness` |
| Speed | `sqrt(2 * g * kinetic_head)` |
| Momentum magnitude | `rho * A * thickness * speed` |
| Submerged | `altitude < water_level` (if a water level is set) |
| Total energy | `rho * g * A * thickness * (thickness/2 + altitude + kinetic_head)` |

### Per-step, per-direction working values

During a step, each cell also holds six values for each of its neighbour directions (arrays of shape `(6, H, W)`):

| Name | Meaning |
|---|---|
| `requested` | Debris thickness the MA wants to send to each neighbour |
| `external` | Part of it that actually crosses the cell face after the shift |
| `received` | Part that arrives after matter loss at the air-water interface |
| `internal` | Part that is shifted but stays inside the cell (`requested - external`) |
| `fraction` | Fraction of the shifted cylinder that crosses the face |
| `travelled` | Shift distance in m |
| `kinetic_head` | Kinetic head of the outflow after shifting and turbulence |
| `*_barycentre_x/y` | Barycentres of the internal and external parts |

## 5. Parameters

### 5.1 Global parameters (`ModelParameters`)

| Parameter | Default | Unit | Meaning |
|---|:-:|:-:|---|
| `pc` | 12.0 | m | Side length of a hexagonal cell |
| `pt` | 2.0 | s | Time represented by one step |
| `rho` | 2000 | kg/m^3 | Density of the debris |
| `rho_water` | 1000 | kg/m^3 | Density of water (for submerged gravity) |
| `pml` | 0.25 | - | Fraction of matter lost when debris enters water (0 to 1) |
| `pwr` | 0.10 | 1/s | Linear water drag coefficient (>= 0) |
| `alpha_scale` | 1.0 | - | Scales the strength of the fictitious-slope correction (>= 0). 0 switches the inertial correction off. |
| `air` | see 5.2 | | Parameters for cells above the water level |
| `water` | see 5.2 | | Parameters for submerged cells |

Derived values:

| Quantity | Formula | Default |
|---|---|:-:|
| `cell_area` | `1.5 * sqrt(3) * pc^2` | about 374 m^2 |
| `centre_spacing` | `sqrt(3) * pc` | about 20.8 m |
| `submerged_gravity` | `g * (rho - rho_water) / rho` | about 4.9 m/s^2 |
| `alpha` | `alpha_scale * 2 * pt / centre_spacing` | about 0.19 s/m |

Model parameters are immutable (frozen dataclasses) and are validated at construction.

### 5.2 Medium parameters (`MediumParameters`, one set for air and one for water)

| Parameter | Air | Water | Unit | Meaning |
|---|:-:|:-:|:-:|---|
| `adh` | 0.02 | 0.02 | m | Adherence: debris thickness that cannot move |
| `pf` | 0.05 | 0.05 | m | Minimum height drop for a neighbour to receive debris (friction threshold in the MA) |
| `pmt` | 2.0 | 2.0 | m | Kinetic head above which erosion starts |
| `per` | 0.05 | 0.02 | - | Progressive erosion coefficient |
| `ptd` | 0.50 | 0.70 | - | Fraction of kinetic head lost to turbulence per step (0 to 1) |
| `ped` | 0.30 | 0.30 | - | Fraction of the excess head lost to erosion (0 to 1) |
| `shift_friction` | 0.05 | 0.10 | - | Friction coefficient in the shift, as `tan(friction angle)` |

Validity: `adh`, `pf`, `pmt`, `per` and `shift_friction` must be finite and >= 0. `ptd` and `ped` must lie in [0, 1].

### 5.3 Global constants

| Constant | Value | Meaning |
|---|:-:|---|
| `GRAVITY` | 9.81 m/s^2 | Gravitational acceleration |
| `EPS` | 1e-9 | Numerical threshold for "zero" |
| `MIN_THICKNESS` | 1e-6 m | Cells thinner than this are cleared at the end of a step |

The default parameter values are illustrative. The paper gives no calibrated values, so they have to be tuned against a real event, for example with the overlap indicator below.

## 6. Model inputs

| Input | Shape | Meaning |
|---|---|---|
| `altitude` | `(H, W)` | Terrain altitude, resampled onto the hexagonal cell centres |
| `thickness` | `(H, W)` | Initial debris thickness (release area) |
| `soil_depth` | `(H, W)` | Erodible soil cover |
| `water_level` | scalar or `None` | Cells with `altitude < water_level` are submerged |

All inputs must be finite. `thickness` and `soil_depth` cannot be negative. The terrain must be sampled at the hex cell centres. Copying a rectangular raster into the array directly shears the terrain.

## 7. Outputs and diagnostics

| Output | Meaning |
|---|---|
| `total_volume()` | Debris volume in m^3 |
| `kinetic_energy()` | Kinetic energy in J |
| `total_energy()` | Total energy from Eq. (1) of the paper |
| `mass_lost` | Matter lost entering water in m^3 |
| `time`, `steps_completed` | Simulated time in s and number of steps |
| `jaccard_fitness(R, S)` | Overlap `\|R and S\| / \|R or S\|` between real and simulated footprints. The paper uses 0.7 as the limit for an acceptable simulation. |

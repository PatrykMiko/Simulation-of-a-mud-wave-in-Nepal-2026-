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
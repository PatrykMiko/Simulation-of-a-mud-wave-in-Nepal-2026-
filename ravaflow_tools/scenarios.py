"""Scenariusze "co by było, gdyby" (S1-S3) i serie przebiegów z jednego pliku YAML.

Każdy wariant to modyfikacja konfiguracji bazowej:
  scale_release: 0.5            # S2: mnożnik grubości masy (objętość skaluje się tak samo)
  barrier: {x1, y1, x2, y2, height, width}   # S1: przeszkoda (tama) dodana do NMT, współrzędne mapy [m]
  params: {rhrelease1: 0.8}     # S3 i inne: nadpisanie parametrów r.avaflow

Przykład pliku: configs/scenarios_synthetic.yaml
Użycie:
    python -m ravaflow_tools.scenarios configs/scenarios_synthetic.yaml          # tylko przygotuj
    python -m ravaflow_tools.scenarios configs/scenarios_synthetic.yaml --run    # przygotuj i uruchom

Warianty uruchamiane są po kolei (jeden mapset GRASS = jeden przebieg naraz).
Uwaga: nowe rastry powstają w <workdir>/scenarios/<nazwa>/ i zachowują nazwy rastrów GRASS.
"""
from __future__ import annotations

import argparse
import copy
import csv
from pathlib import Path

import numpy as np
import yaml

from . import metrics as M
from .runner import load_config, run
from .score import score_run


def _raster_path(cfg: dict, grass_name: str) -> Path:
    return Path(cfg["workdir"]).expanduser().resolve() / cfg["inputs"][grass_name]


def _release_names(cfg: dict) -> list[str]:
    """Nazwy rastrów masy (parametry hrelease, hrelease1, ...) obecne w `inputs`."""
    names = []
    for key, val in cfg["params"].items():
        if key.startswith("hrelease") and isinstance(val, str) and val in cfg["inputs"]:
            names.append(val)
    return names


def _swap_input(cfg: dict, grass_name: str, new_file: Path, variant_dir: Path) -> None:
    workdir = Path(cfg["workdir"]).expanduser().resolve()
    cfg["inputs"][grass_name] = str(new_file.relative_to(workdir))


def scale_release(cfg: dict, factor: float, variant_dir: Path) -> dict:
    """S2: mnoży rastry masy przez `factor`. Zwraca nową konfigurację."""
    cfg = copy.deepcopy(cfg)
    names = _release_names(cfg)
    if not names:
        raise ValueError("Brak rastrów masy (hrelease*) wśród inputs - nie ma czego skalować "
                         "(jeśli powstają w pre_commands, skaluj je tam).")
    variant_dir.mkdir(parents=True, exist_ok=True)
    for name in names:
        src = _raster_path(cfg, name)
        arr, meta = M.read_raster(src)
        out = variant_dir / src.name
        valid = arr != meta["nodata"]
        M.write_raster(out, np.where(valid, arr * factor, arr), meta)
        _swap_input(cfg, name, out, variant_dir)
    return cfg


def add_barrier(cfg: dict, barrier: dict, variant_dir: Path) -> dict:
    """S1: podnosi NMT o `height` [m] w pasie szerokości `width` [m] wokół odcinka
    (x1,y1)-(x2,y2) we współrzędnych mapy. Zwraca nową konfigurację."""
    cfg = copy.deepcopy(cfg)
    elev_name = cfg["params"]["elevation"]
    src = _raster_path(cfg, elev_name)
    z, meta = M.read_raster(src)
    cs = meta["cellsize"]
    cols, rows = np.meshgrid(np.arange(meta["ncols"]), np.arange(meta["nrows"]))
    px = meta["xll"] + (cols + 0.5) * cs
    py = meta["yll"] + (meta["nrows"] - rows - 0.5) * cs
    x1, y1, x2, y2 = (float(barrier[k]) for k in ("x1", "y1", "x2", "y2"))
    dx, dy = x2 - x1, y2 - y1
    seg2 = dx * dx + dy * dy
    t = np.clip(((px - x1) * dx + (py - y1) * dy) / seg2, 0.0, 1.0) if seg2 > 0 else np.zeros_like(px)
    dist = np.hypot(px - (x1 + t * dx), py - (y1 + t * dy))
    wall = dist <= float(barrier["width"]) / 2.0
    valid = z != meta["nodata"]
    variant_dir.mkdir(parents=True, exist_ok=True)
    out = variant_dir / src.name
    M.write_raster(out, np.where(wall & valid, z + float(barrier["height"]), z), meta)
    _swap_input(cfg, elev_name, out, variant_dir)
    return cfg


def make_variant(base: dict, spec: dict) -> dict:
    """Buduje konfigurację wariantu wg specyfikacji (name + scale_release/barrier/params)."""
    name = spec["name"]
    workdir = Path(base["workdir"]).expanduser().resolve()
    vdir = workdir / "scenarios" / name
    cfg = copy.deepcopy(base)
    if "scale_release" in spec:
        cfg = scale_release(cfg, float(spec["scale_release"]), vdir)
    if "barrier" in spec:
        cfg = add_barrier(cfg, spec["barrier"], vdir)
    for k, v in (spec.get("params") or {}).items():
        cfg["params"][k] = v
    cfg["params"]["prefix"] = f"{base['params']['prefix']}_{name}"
    return cfg


def run_variants(base: dict, specs: list[dict], do_run: bool = False) -> list[dict]:
    """Przygotowuje (i opcjonalnie uruchamia) warianty. Zwraca listę wierszy podsumowania."""
    rows = []
    for spec in specs:
        cfg = make_variant(base, spec)
        row = {"variant": spec["name"], "prefix": cfg["params"]["prefix"]}
        if do_run:
            res = run(cfg)
            row.update(model_seconds=res.get("model_seconds"), wall_seconds=res["wall_seconds"])
            if "reference_mask" in cfg:
                sc = score_run(cfg)
                row.update({k: sc[k] for k in ("jaccard", "precision", "recall", "area_ratio", "hausdorff_m")})
        rows.append(row)
    if do_run and rows:
        out = Path(base["workdir"]).expanduser().resolve() / "scenarios_summary.csv"
        keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("variant", "prefix"), k))
        with open(out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            w.writerows(rows)
        print(f"Podsumowanie: {out}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("scenarios_yaml")
    ap.add_argument("--run", action="store_true", help="uruchom przebiegi (domyślnie tylko przygotowanie)")
    a = ap.parse_args()
    with open(a.scenarios_yaml, encoding="utf-8") as f:
        sc = yaml.safe_load(f)
    base_cfg = load_config(sc["base"])
    for r in run_variants(base_cfg, sc["variants"], do_run=a.run):
        print(r)

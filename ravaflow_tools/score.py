"""Ocena przebiegu r.avaflow względem zasięgu rzeczywistego.

Użycie:  python -m ravaflow_tools.score configs/ca_exp2.yaml
Wymaga w konfiguracji: reference_mask (plik .tif/.asc z maską zasięgu).
Wyniki: <workdir>/<prefix>_score/{metrics.json, extent_comparison.png, hmax.png}
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

from . import metrics as M
from . import viz
from .runner import load_config


def score_run(cfg: dict) -> dict:
    workdir = Path(cfg["workdir"]).expanduser().resolve()
    prefix = cfg["params"]["prefix"]
    ascii_dir = workdir / f"{prefix}_results" / f"{prefix}_ascii"
    threshold = float(cfg.get("threshold", 0.1))

    cands = sorted(ascii_dir.glob(f"{prefix}_hflow_max*.asc"))
    if not cands:
        raise FileNotFoundError(f"Nie znaleziono {prefix}_hflow_max*.asc w {ascii_dir}. "
                                f"Dostępne pliki bez numerów kroków: "
                                f"{[p.name for p in ascii_dir.iterdir() if not any(ch.isdigit() for ch in p.stem[-4:])][:20]}")
    hmax, meta = M.read_ascii(cands[0])
    elev, _ = M.read_ascii(ascii_dir / f"{prefix}_elev.asc")

    ref_path = workdir / cfg["reference_mask"]
    ref, ref_meta = M.read_raster(ref_path)
    if not M.same_grid(meta, ref_meta):
        diffs = M.grid_differences(meta, ref_meta)
        ref = M.sample_to_grid(ref, ref_meta, meta)  # najbliższy sąsiad na siatkę modelu
        print("[uwaga] siatki modelu i referencji różnią się (model vs referencja): "
              + "; ".join(diffs)
              + ". Referencję próbkuję najbliższym sąsiadem na siatkę modelu; komórki "
                "poza obszarem rastra referencji są pomijane w ocenie.")
    obs, valid = M.reference_masks(ref, {**meta, "nodata": ref_meta["nodata"]})
    model = M.model_extent(hmax, meta, threshold)

    result = M.score(model, obs, valid, meta["cellsize"])
    result["threshold_m"] = threshold

    out = workdir / f"{prefix}_score"
    out.mkdir(exist_ok=True)
    (out / "metrics.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    viz.plot_extent_comparison(elev, model, obs, valid, meta, out / "extent_comparison.png",
                               title=f"{prefix}: zasięg modelu vs obserwowany (próg {threshold} m)",
                               metrics=result)
    viz.plot_hmax(elev, hmax, meta, out / "hmax.png", title=f"{prefix}: maksymalna grubość przepływu",
                  threshold=threshold)
    return result


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Użycie: python -m ravaflow_tools.score configs/plik.yaml")
    r = score_run(load_config(sys.argv[1]))
    print(f"Jaccard e = {r['jaccard']:.3f} | precision {r['precision']:.3f} | recall {r['recall']:.3f} "
          f"| F1 {r['f1']:.3f} | pole model/obs {r['area_ratio']:.2f} | Hausdorff {r['hausdorff_m']:.0f} m")

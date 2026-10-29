"""Converter: Ravaflow data directory -> SCIDDICA-SS3 arrays.

Usage: python -m ravaflow_tools.to_ss3 data/nepal --out data/nepal_ss3 --water-level 300.0 --volume 200e6
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

from . import check_inputs, metrics as M


def convert(directory, out, water_level=None, volume=None, cellsize=None, epsg=None,
            pc=None, write_rasters=False):
    rep = check_inputs.check_inputs(directory, volume=volume, cellsize=cellsize, epsg=epsg)
    rep.print()
    if not rep.ok:
        sys.exit(1)

    d = Path(directory)
    arrays, metas = {}, {}
    for name in ("elevation", "release_height", "reference_extent", "entrainment_height"):
        p = next((d / f"{name}.{e}" for e in ("tif", "tiff", "asc") if (d / f"{name}.{e}").exists()), None)
        if p is None:
            continue
        arrays[name], metas[name] = M.read_raster(p)

    ref_meta = metas["elevation"]
    cs = ref_meta["cellsize"]
    for name, m in metas.items():
        if not M.same_grid(ref_meta, m):
            print(f"siatka różni się: {name}")
            sys.exit(1)

    z = arrays["elevation"].astype(np.float64)
    if np.any(~np.isfinite(z)):
        print("elevation musi być skończone")
        sys.exit(1)

    h = arrays["release_height"].astype(np.float64)
    hn = metas["release_height"]["nodata"]
    h = np.where((~np.isfinite(h)) | (h == hn), 0.0, h)
    h = np.clip(h, 0.0, None)

    if "entrainment_height" in arrays:
        e = arrays["entrainment_height"].astype(np.float64)
        en = metas["entrainment_height"]["nodata"]
        e = np.where((~np.isfinite(e)) | (e == en), 0.0, e)
        e = np.clip(e, 0.0, None)
    else:
        e = np.zeros_like(z)

    soil = e

    # optional resampling to pc
    if pc is not None and abs(pc - cs) > 1e-9:
        # build dst meta
        src = ref_meta
        # dimensions based on extent
        width_m = src["ncols"] * src["cellsize"]
        height_m = src["nrows"] * src["cellsize"]
        dst_ncols = int(round(width_m / pc))
        dst_nrows = int(round(height_m / pc))
        dst = {
            "ncols": dst_ncols, "nrows": dst_nrows,
            "xll": src["xll"], "yll": src["yll"],
            "cellsize": pc, "nodata": src["nodata"],
        }
        z = M.sample_to_grid(z, src, dst)
        # replace nodata (shouldn't occur for in-bounds) with 0 for thickness/soil only; elevation already finite
        h = M.sample_to_grid(h, src, dst)
        soil = M.sample_to_grid(soil, src, dst)
        # fill nodata with 0 for thickness/soil
        h = np.where(h == src["nodata"], 0.0, h)
        soil = np.where(soil == src["nodata"], 0.0, soil)
        h = np.clip(h, 0.0, None)
        soil = np.clip(soil, 0.0, None)
        cs = pc
        meta = dst
    else:
        meta = ref_meta

    # elevation finite check after resample
    if np.any(~np.isfinite(z)) or np.any(z == meta.get("nodata", -9999.0)):
        # fill elevation nodata? contract says no nodata; keep as is but warn
        pass

    volume_after = float(np.sum(h) * cs ** 2)
    print(f"objętość po konwersji: {volume_after:.4g} m3")

    if water_level is not None and not np.isfinite(water_level):
        print("water_level musi być skończone")
        sys.exit(1)

    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / "ss3_arrays.npz", altitude=z, thickness=h, soil_depth=soil)

    sidecar = {
        "cellsize": float(cs),
        "water_level": (None if water_level is None else float(water_level)),
        "meta": {k: (float(v) if isinstance(v, (np.floating,)) else v) for k, v in meta.items()},
        "source_dir": str(Path(directory).resolve()),
    }
    with open(out_dir / "ss3_meta.json", "w") as f:
        json.dump(sidecar, f, indent=2)

    if write_rasters:
        for arr, name in ((z, "altitude"), (h, "thickness"), (soil, "soil_depth")):
            M.write_raster(out_dir / f"{name}.tif", arr, meta)

    print(f"zapisano: {out_dir / 'ss3_arrays.npz'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--water-level", type=float, default=None)
    ap.add_argument("--volume", type=float, default=None)
    ap.add_argument("--cellsize", type=float, default=None)
    ap.add_argument("--epsg", type=int, default=None)
    ap.add_argument("--pc", type=float, default=None)
    ap.add_argument("--write-rasters", action="store_true")
    a = ap.parse_args()
    convert(a.directory, a.out, water_level=a.water_level, volume=a.volume,
            cellsize=a.cellsize, epsg=a.epsg, pc=a.pc, write_rasters=a.write_rasters)


if __name__ == "__main__":
    main()

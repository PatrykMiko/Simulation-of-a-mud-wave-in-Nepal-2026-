"""Converter: Ravaflow data directory -> SCIDDICA-SS3 arrays on a hexagonal grid.

The rasters are resampled onto the centres of hexagonal cells (no shear), and
cells outside the data or holding a constant fill value become walls.

Usage:
    python -m ravaflow_tools.to_ss3 data/nepal --out data/nepal_ss3 \
        --water-level 300.0 --volume 200e6 --fill-value 0
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from scipy.ndimage import distance_transform_edt

from . import check_inputs, metrics as M
from .hexgrid import to_hex

HEX_AREA_FACTOR = 1.5 * np.sqrt(3.0)          # hexagon area = factor * side^2
DEFAULT_PC_RATIO = 0.62                       # hexagon of the same area as one raster pixel
WALL_MARGIN = 1000.0                          # walls are this far above the highest terrain [m]


def _find_raster(directory: Path, name: str):
    return next(
        (directory / f"{name}.{ext}" for ext in ("tif", "tiff", "asc")
         if (directory / f"{name}.{ext}").exists()),
        None,
    )


def _clean(array, nodata, clip: bool = True) -> np.ndarray:
    """Float array with non-finite and nodata cells set to 0 (and negatives clipped)."""
    a = np.asarray(array, dtype=np.float64)
    bad = ~np.isfinite(a)
    if nodata is not None and np.isfinite(nodata):
        bad |= a == nodata
    a = np.where(bad, 0.0, a)
    return np.clip(a, 0.0, None) if clip else a


def _json_default(obj):
    if hasattr(obj, "item"):                   # numpy scalar
        return obj.item()
    return str(obj)


def convert(directory, out, water_level=None, volume=None, cellsize=None, epsg=None,
            pc=None, write_rasters=False, fill_value=None):
    rep = check_inputs.check_inputs(directory, volume=volume, cellsize=cellsize, epsg=epsg)
    rep.print()
    if not rep.ok:
        sys.exit(1)

    d = Path(directory)
    arrays, metas = {}, {}
    for name in ("elevation", "release_height", "reference_extent", "entrainment_height"):
        path = _find_raster(d, name)
        if path is None:
            continue
        arrays[name], metas[name] = M.read_raster(path)

    ref_meta = metas["elevation"]
    cs = float(ref_meta["cellsize"])
    for name, m in metas.items():
        if not M.same_grid(ref_meta, m):
            print(f"siatka różni się: {name}")
            sys.exit(1)

    z = arrays["elevation"].astype(np.float64)
    if np.any(~np.isfinite(z)):
        print("elevation musi być skończone")
        sys.exit(1)

    if water_level is not None and not np.isfinite(water_level):
        print("water_level musi być skończone")
        sys.exit(1)

    h = _clean(arrays["release_height"], metas["release_height"].get("nodata"))

    if "entrainment_height" in arrays:
        soil = _clean(arrays["entrainment_height"], metas["entrainment_height"].get("nodata"))
    else:
        soil = np.zeros_like(z)

    reference = None
    if "reference_extent" in arrays:
        reference = _clean(arrays["reference_extent"],
                           metas["reference_extent"].get("nodata")) > 0

    source_volume = float(h.sum() * cs ** 2)

    # ---- cells that carry no data: constant fill value or the raster's nodata ----
    blocked_src = np.zeros(z.shape, dtype=bool)
    if fill_value is not None:
        blocked_src |= z == fill_value
    nodata = ref_meta.get("nodata")
    if nodata is not None and np.isfinite(nodata):
        blocked_src |= z == nodata
    if blocked_src.all():
        print("cała siatka jest wypełnieniem (fill) - brak danych terenu")
        sys.exit(1)

    if blocked_src.any():
        print(f"komórki bez danych: {blocked_src.mean():.1%} siatki -> ściany")
        # Replace them by the nearest real altitude so interpolation makes no ramps
        # at the edge of the data.
        nearest = np.asarray(
            distance_transform_edt(blocked_src, return_distances=False, return_indices=True),
            dtype=np.intp,
        )
        z = z[tuple(nearest)]
        del nearest

    # ---- resample everything onto the hexagonal cell centres ----------------------
    pc = float(pc) if pc is not None else DEFAULT_PC_RATIO * cs
    if not np.isfinite(pc) or pc <= 0:
        print("pc musi być dodatnie")
        sys.exit(1)

    z_hex, valid, xoff = to_hex(z, ref_meta, pc)
    h_hex, _, _ = to_hex(h, ref_meta, pc)
    s_hex, _, _ = to_hex(soil, ref_meta, pc)

    if blocked_src.any():
        blocked_hex = to_hex(blocked_src.astype(np.float64), ref_meta, pc, order=0)[0] > 0.5
    else:
        blocked_hex = np.zeros_like(valid)

    usable = valid & ~blocked_hex                   # cells that are real terrain
    wall = float(z.max()) + WALL_MARGIN

    altitude = np.where(usable, z_hex, wall)
    thickness = np.where(usable, np.clip(h_hex, 0.0, None), 0.0)
    soil_depth = np.where(usable, np.clip(s_hex, 0.0, None), 0.0)

    arrays_out = {
        "altitude": altitude,
        "thickness": thickness,
        "soil_depth": soil_depth,
        "domain": usable,
    }
    if reference is not None:
        ref_hex = to_hex(reference.astype(np.float64), ref_meta, pc, order=0)[0] > 0.5
        arrays_out["reference"] = ref_hex & usable

    # ---- report ---------------------------------------------------------------------
    hex_area = HEX_AREA_FACTOR * pc ** 2
    volume_after = float(thickness.sum() * hex_area)
    print(f"siatka heksagonalna: {altitude.shape[0]} x {altitude.shape[1]}, "
          f"pc = {pc:.2f} m (bok), pole komórki = {hex_area:.1f} m2")
    print(f"objętość źródłowa:     {source_volume:.4g} m3")
    print(f"objętość po konwersji: {volume_after:.4g} m3")
    if source_volume > 0:
        diff = (volume_after - source_volume) / source_volume
        print(f"różnica objętości: {diff:+.2%}")
        if abs(diff) > 0.05:
            print("UWAGA: objętość zmieniła się o więcej niż 5% (próbkowanie lub ściany w obszarze źródła)")
    if volume is not None and volume > 0:
        print(f"objętość zadana (--volume): {volume:.4g} m3, "
              f"różnica po konwersji: {(volume_after - volume) / volume:+.2%}")

    # ---- save -----------------------------------------------------------------------
    out_dir = Path(out)
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / "ss3_arrays.npz", **arrays_out)

    height_m = ref_meta["nrows"] * cs
    meta = {k: v for k, v in ref_meta.items()}
    meta.update(
        hex=True,
        xoff=float(xoff),
        origin_e=float(ref_meta["xll"] - xoff),
        origin_n=float(ref_meta["yll"] + height_m),
    )

    sidecar = {
        "cellsize": float(pc),                      # hexagon side length, used as pc by the model
        "source_cellsize": float(cs),
        "hex": True,
        "water_level": None if water_level is None else float(water_level),
        "meta": meta,
        "source_dir": str(d.resolve()),
    }
    with open(out_dir / "ss3_meta.json", "w") as f:
        json.dump(sidecar, f, indent=2, default=_json_default)

    if write_rasters:
        print("--write-rasters: pominięto, tablice heksagonalne nie są regularnym rastrem")

    print(f"zapisano: {out_dir / 'ss3_arrays.npz'}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory")
    ap.add_argument("--out", required=True)
    ap.add_argument("--water-level", type=float, default=None)
    ap.add_argument("--volume", type=float, default=None)
    ap.add_argument("--cellsize", type=float, default=None)
    ap.add_argument("--epsg", type=int, default=None)
    ap.add_argument("--pc", type=float, default=None,
                    help="hexagon side length in m (default: 0.62 x raster cell size)")
    ap.add_argument("--fill-value", type=float, default=None,
                    help="altitude used as 'no data' outside the DEM; those cells become walls")
    ap.add_argument("--write-rasters", action="store_true")
    a = ap.parse_args()
    convert(a.directory, a.out, water_level=a.water_level, volume=a.volume,
            cellsize=a.cellsize, epsg=a.epsg, pc=a.pc, write_rasters=a.write_rasters,
            fill_value=a.fill_value)


if __name__ == "__main__":
    main()
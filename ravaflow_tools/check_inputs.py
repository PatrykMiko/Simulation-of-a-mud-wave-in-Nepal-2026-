"""Walidator danych wejściowych (kontrakt danych projektu).

Użycie:
    python -m ravaflow_tools.check_inputs data/nepal --volume 200e6 --cellsize 30 --epsg 32645

Sprawdza: obecność plików, identyczną siatkę, układ współrzędnych (GeoTIFF),
brak wartości brakujących w elevation i release_height, release_height = 0 poza
obszarem odrywania, zgodność objętości, konwencję maski zasięgu (1/0/ujemne).
Kod wyjścia 1, jeśli są błędy (ostrzeżenia nie psują wyniku).
"""
from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from . import metrics as M

REQUIRED = ("elevation", "release_height", "reference_extent")
OPTIONAL = ("entrainment_height",)


@dataclass
class Report:
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    info: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.errors

    def print(self) -> None:
        for tag, items in (("BŁĄD", self.errors), ("UWAGA", self.warnings), ("INFO", self.info)):
            for m in items:
                print(f"[{tag}] {m}")
        print("\nWYNIK:", "OK - dane zgodne z kontraktem" if self.ok else f"{len(self.errors)} błęd(ów) do poprawy")


def _find(directory: Path, name: str) -> Path | None:
    for ext in ("tif", "tiff", "asc"):
        p = directory / f"{name}.{ext}"
        if p.exists():
            return p
    return None


def _crs_info(path: Path):
    """(is_projected, linear_unit_metre, epsg) dla GeoTIFF; None dla .asc."""
    if path.suffix.lower() == ".asc":
        return None
    import rasterio

    with rasterio.open(path) as ds:
        crs = ds.crs
        if crs is None:
            return (False, False, None)
        unit = (crs.linear_units or "").lower()
        return (crs.is_projected, unit in ("metre", "meter", "m"), crs.to_epsg())


def check_inputs(directory: str | Path, volume: float | None = None, cellsize: float | None = None,
                 epsg: int | None = None, volume_tol: float = 0.01) -> Report:
    rep = Report()
    d = Path(directory)
    arrays, metas, paths = {}, {}, {}
    for name in REQUIRED + OPTIONAL:
        p = _find(d, name)
        if p is None:
            if name in REQUIRED:
                rep.errors.append(f"brak pliku {name}.tif (lub .asc) w {d}")
            continue
        try:
            arrays[name], metas[name] = M.read_raster(p)
            paths[name] = p
        except Exception as e:  # noqa: BLE001
            rep.errors.append(f"{p.name}: nie można wczytać ({e})")
    if "elevation" not in arrays:
        return rep

    ref_meta = metas["elevation"]
    cs = ref_meta["cellsize"]
    rep.info.append(f"siatka: {ref_meta['ncols']} x {ref_meta['nrows']} komórek, rozmiar {cs} m")

    # 1. ta sama siatka
    for name, m in metas.items():
        if not M.same_grid(ref_meta, m):
            rep.errors.append(f"{paths[name].name}: siatka różni się od elevation "
                              f"({m['ncols']}x{m['nrows']}, {m['cellsize']} m, xll={m['xll']}, yll={m['yll']})")
    # 2. rozmiar komórki
    if cellsize is not None and abs(cs - cellsize) > 1e-6:
        rep.errors.append(f"rozmiar komórki {cs} m, oczekiwano {cellsize} m")
    # 3. układ współrzędnych
    for name, p in paths.items():
        info = _crs_info(p)
        if info is None:
            rep.warnings.append(f"{p.name}: format .asc nie niesie układu współrzędnych - nie sprawdzono")
            continue
        projected, metre, found = info
        if not projected or not metre:
            rep.errors.append(f"{p.name}: układ musi być rzutowany w metrach (jest: "
                              f"{'brak CRS' if found is None and not projected else 'geograficzny lub inna jednostka'})")
        elif epsg is not None and found != epsg:
            rep.errors.append(f"{p.name}: EPSG:{found}, oczekiwano EPSG:{epsg}")

    nd = ref_meta["nodata"]
    # 4. elevation
    z = arrays["elevation"]
    bad = (~np.isfinite(z)) | (z == nd)
    if bad.any():
        rep.errors.append(f"elevation: {int(bad.sum())} komórek bez wartości (NaN/nodata) - uzupełnijcie dziury")
    zv = z[~bad]
    if zv.size:
        rep.info.append(f"elevation: {zv.min():.0f}..{zv.max():.0f} m n.p.m.")
        if zv.min() < -500 or zv.max() > 9000:
            rep.warnings.append("elevation: wartości poza zakresem -500..9000 m - sprawdźcie jednostki")

    # 5. release_height
    if "release_height" in arrays:
        h = arrays["release_height"]
        hn = metas["release_height"]["nodata"]
        hbad = (~np.isfinite(h)) | (h == hn)
        if hbad.any():
            rep.errors.append(f"release_height: {int(hbad.sum())} komórek bez wartości - poza obszarem odrywania "
                              f"ma być dokładnie 0, a nie nodata")
        hv = np.where(hbad, 0.0, h)
        if (hv < 0).any():
            rep.errors.append("release_height: wartości ujemne")
        pos = hv > 0
        if not pos.any():
            rep.errors.append("release_height: brak komórek z masą (same zera)")
        else:
            vol = float(hv[pos].sum() * cs ** 2)
            rep.info.append(f"release_height: {int(pos.sum())} komórek, maks. grubość {hv.max():.1f} m, "
                            f"objętość {vol:.4g} m3")
            if volume is not None:
                if abs(vol - volume) / volume > volume_tol:
                    rep.errors.append(f"objętość z rastra {vol:.4g} m3 różni się od zadanej {volume:.4g} m3 "
                                      f"o {100 * abs(vol - volume) / volume:.1f}% (tolerancja {100 * volume_tol:.0f}%)")
            else:
                rep.warnings.append("nie podano --volume: porównajcie objętość z rastra z założoną")
            if hv.max() > 1000:
                rep.warnings.append("release_height: grubość > 1000 m - sprawdźcie jednostki")
            tiny = (hv > 0) & (hv < 0.01)
            if tiny.any():
                rep.warnings.append(f"release_height: {int(tiny.sum())} komórek z grubością < 1 cm "
                                    f"(szum? poza obszarem odrywania ma być 0)")
            if zv.size and "elevation" in arrays and bad[pos].any():
                rep.errors.append("elevation ma dziury pod obszarem odrywania")

    # 6. reference_extent
    if "reference_extent" in arrays:
        r = arrays["reference_extent"]
        rn = metas["reference_extent"]["nodata"]
        valid = (r != rn) & (r >= 0) & np.isfinite(r)
        vals = np.unique(r[valid])
        if not set(np.round(vals, 6)).issubset({0.0, 1.0}):
            rep.errors.append(f"reference_extent: dozwolone 1 (zasięg), 0 (brak), ujemne/nodata (nieznane); "
                              f"znaleziono wartości {vals[:6].tolist()}")
        n_pos = int((r[valid] > 0).sum())
        if n_pos == 0:
            rep.errors.append("reference_extent: brak komórek z zasięgiem (same 0)")
        else:
            rep.info.append(f"reference_extent: zasięg {n_pos} komórek = {n_pos * cs ** 2 / 1e6:.3f} km2; "
                            f"{int((~valid).sum())} komórek bez informacji")
        if "release_height" in arrays and 'pos' in locals() and pos.any() and n_pos:
            inside = float(((r > 0) & pos).sum()) / float(pos.sum())
            if inside < 0.5:
                rep.warnings.append(f"tylko {100 * inside:.0f}% obszaru odrywania leży w zasięgu rzeczywistym - "
                                    f"sprawdźcie, czy rastry są wyrównane")

    # 7. entrainment_height (opcjonalny)
    if "entrainment_height" in arrays:
        e = arrays["entrainment_height"]
        en = metas["entrainment_height"]["nodata"]
        if ((~np.isfinite(e)) | (e == en)).any():
            rep.errors.append("entrainment_height: wartości brakujące - poza obszarem ma być 0")
        elif (e < 0).any():
            rep.errors.append("entrainment_height: wartości ujemne")
    return rep


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("directory")
    ap.add_argument("--volume", type=float, default=None, help="założona objętość [m3], np. 200e6")
    ap.add_argument("--cellsize", type=float, default=None, help="oczekiwany rozmiar komórki [m]")
    ap.add_argument("--epsg", type=int, default=None, help="oczekiwany kod EPSG, np. 32645")
    a = ap.parse_args()
    rep = check_inputs(a.directory, a.volume, a.cellsize, a.epsg)
    rep.print()
    sys.exit(0 if rep.ok else 1)

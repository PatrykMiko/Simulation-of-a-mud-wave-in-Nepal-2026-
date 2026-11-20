"""
Kontrola wyników prepare_nepal_data.py (zgodność z kontraktem danych i spójność DEM-ów).

Uruchomienie (z głównego folderu repo), te same --version i --volume co przy budowaniu:
    python -m ravaflow_tools.verify_derived --version v1 --volume 140e6
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
NEPAL = ROOT / "data" / "nepal"
DEMS = ("copernicus", "alos", "srtm")
FILES = ("elevation", "release_height", "reference_extent", "entrainment_height")
CELL, EPSG, MIN_ELEV, H_WARN = 30.0, 32645, 100.0, 150.0


def check_arrays(a: dict, volume: float, cell: float = CELL, mode: str = "level", base_m: float = 0.0):
    """a: nazwa -> tablica; mode: "level" (pozioma tafla) lub "hybrid" (baza base_m + pozioma tafla).
    Zwraca listę (status, opis), status: OK / BŁĄD / UWAGA."""
    out = []

    def add(ok, msg, warn=False):
        out.append(("OK" if ok else ("UWAGA" if warn else "BŁĄD"), msg))

    elev, rel = a["elevation"], a["release_height"]
    ref, ent = a["reference_extent"], a["entrainment_height"]
    for k, arr in a.items():
        add(bool(np.isfinite(arr).all()), f"{k}: brak NaN/inf")

    dom = ref >= 0
    wet = rel > 0
    add(bool(np.isin(ref, [-1.0, 0.0, 1.0]).all()), "reference_extent: tylko wartości -1, 0, 1")
    add(bool(rel.min() >= 0), "release_height: brak wartości ujemnych")
    vol = float(rel.sum()) * cell * cell
    add(abs(vol - volume) <= 1e-3 * volume,
        f"release_height: objętość {vol / 1e6:.3f} mln m3 (cel {volume / 1e6:.3f})")
    add(bool(wet.any()), f"release_height: {int(wet.sum())} komórek z masą")
    if wet.any():
        add(bool(dom[wet].all()), "release_height: cała masa w obszarze roboczym")
        if mode == "level":
            spread = float((elev + rel)[wet].std())
            add(spread < 0.01, f"powierzchnia masy pozioma (odch. std {spread:.4f} m)")
        else:
            add(float(rel[wet].min()) >= base_m - 0.01,
                f"warstwa bazowa: najmniejsza grubość w polygonie {float(rel[wet].min()):.1f} m (baza {base_m:g} m)")
            pool = rel > base_m + 0.01
            add(bool(pool.any()), f"tafla ponad bazą: {int(pool.sum())} komórek")
            if pool.any():
                spread = float((elev + rel)[pool].std())
                add(spread < 0.01, f"powierzchnia tafli pozioma (odch. std {spread:.4f} m)")
        hmax = float(rel.max())
        add(hmax <= H_WARN, f"maks. grubość masy {hmax:.1f} m (ostrzeżenie od {H_WARN:.0f} m)", warn=True)
    add(bool(dom.any()), f"obszar roboczy: {int(dom.sum())} komórek")
    if dom.any():
        emin = float(elev[dom].min())
        add(emin >= MIN_ELEV, f"elevation: min w obszarze {emin:.0f} m (próg {MIN_ELEV:.0f} m)")
    add(bool((elev[~dom] == 0).all()), "elevation: poza obszarem dokładnie 0")
    add(bool(np.isin(ent, [0.0, 0.5, 1.5, 3.0]).all()), "entrainment_height: tylko wartości 0, 0.5, 1.5, 3.0")
    add(bool((ent[~dom] == 0).all()), "entrainment_height: poza obszarem 0")
    add(bool((rel[~dom] == 0).all()), "release_height: poza obszarem dokładnie 0")
    return out


def load(dem_dir: Path):
    import rasterio

    arrs, metas = {}, {}
    for k in FILES:
        p = dem_dir / f"{k}.tif"
        if not p.exists():
            return None, None, f"brak pliku {p}"
        with rasterio.open(p) as f:
            arrs[k] = f.read(1).astype("float64")
            epsg = f.crs.to_epsg() if f.crs else None
            metas[k] = (tuple(f.transform)[:6], epsg, f.res, f.nodata)
    return arrs, metas, None


def main() -> int:
    import rasterio

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--version", required=True)
    ap.add_argument("--volume", required=True, type=float)
    args = ap.parse_args()

    root = NEPAL / "derived" / args.version
    n_err = n_warn = 0
    refs = {}

    def show(status, msg):
        nonlocal n_err, n_warn
        n_err += status == "BŁĄD"
        n_warn += status == "UWAGA"
        print(f"  [{status}] {msg}")

    for name in DEMS:
        print(f"\n=== {name} ===")
        arrs, metas, err = load(root / name)
        if err:
            show("BŁĄD", err)
            continue

        # siatka: identyczna we wszystkich plikach, 30 m, EPSG:32645
        tr = {m[0] for m in metas.values()}
        shapes = {a.shape for a in arrs.values()}
        show("OK" if len(tr) == 1 and len(shapes) == 1 else "BŁĄD",
             f"identyczna siatka we wszystkich 4 plikach {sorted(shapes)[0] if len(shapes) == 1 else shapes}")
        res_ok = all(abs(m[2][0] - CELL) < 1e-6 and abs(m[2][1] - CELL) < 1e-6 for m in metas.values())
        epsg_ok = all(m[1] == EPSG for m in metas.values())
        show("OK" if res_ok else "BŁĄD", f"rozmiar komórki {CELL:g} m")
        show("OK" if epsg_ok else "BŁĄD", f"układ EPSG:{EPSG}")
        if len(shapes) != 1:
            continue

        mf = root / name / "manifest.json"
        manifest = json.loads(mf.read_text(encoding="utf-8")) if mf.exists() else {}
        mode = manifest.get("release_mode", "level")
        print(f"  [INFO] tryb masy źródłowej: {mode}")
        for st, msg in check_arrays(arrs, args.volume, mode=mode, base_m=manifest.get("base_m") or 0.0):
            show(st, msg)

        ref, rel, elev, ent = (arrs[k] for k in ("reference_extent", "release_height",
                                                 "elevation", "entrainment_height"))
        dom, wet = ref >= 0, rel > 0

        # zgodność z dem_original.tif
        p = NEPAL / name / "dem_original.tif"
        if p.exists():
            with rasterio.open(p) as f:
                d = f.read(1).astype("float64")
            same = d.shape == elev.shape and np.array_equal(d[dom], elev[dom])
            show("OK" if same else "BŁĄD", "elevation w obszarze = dem_original.tif (bez odejmowania masy)")

        # informacje do oceny "na oko"
        if wet.any():
            rows, cols = np.nonzero(wet)
            t = metas["elevation"][0]
            x = t[2] + (cols.mean() + 0.5) * t[0]
            y = t[5] + (rows.mean() + 0.5) * t[4]
            print(f"  [INFO] środek źródła (UTM 45N): x = {x:.0f}, y = {y:.0f}; "
                  f"pole {wet.sum() * CELL * CELL / 1e6:.3f} km2")
            print(f"  [INFO] komórek źródła w zasięgu fali (reference_extent = 1): "
                  f"{100 * (ref[wet] == 1).mean():.0f}%")
        print(f"  [INFO] zasięg fali: {int((ref == 1).sum())} komórek "
              f"({100 * (ref == 1).sum() / max(dom.sum(), 1):.1f}% obszaru roboczego)")
        print(f"  [INFO] maks. objętość do porwania: {ent.sum() * CELL * CELL / 1e6:.1f} mln m3")

        if manifest:
            z = manifest.get("z_level_m")
            print(f"  [INFO] manifest: {'Z = %.1f m n.p.m., ' % z if z is not None else ''}"
                  f"h_max = {manifest['h_max_m']:.1f} m")
        refs[name] = (ref, ent, metas["elevation"][0])

    # spójność między DEM-ami: te same rastry zasięgu i porywania, ta sama siatka
    print("\n=== spójność między DEM-ami ===")
    if len(refs) == len(DEMS):
        base = refs[DEMS[0]]
        for name in DEMS[1:]:
            r, e, t = refs[name]
            show("OK" if np.array_equal(r, base[0]) and np.array_equal(e, base[1]) and t == base[2]
                 else "BŁĄD", f"{name} ma ten sam reference_extent, entrainment_height i siatkę co {DEMS[0]}")
    else:
        show("BŁĄD", "nie wszystkie DEM-y zostały wczytane")

    print(f"\nPodsumowanie: błędów {n_err}, uwag {n_warn}")
    return 1 if n_err else 0


if __name__ == "__main__":
    sys.exit(main())

"""
Przygotowanie danych wejściowych do symulacji fali błotnej (Nepal, 26.08.2026).

Skrypt robi automatycznie to, co wcześniej robiono ręcznie w QGIS:
  1. wczytuje surowe DEM-y (Copernicus COP30, ALOS AW3D30, SRTM GL1),
  2. buduje obszar roboczy: bufor wokół osi rzeki (warstwa `rzeka`, 2 km z każdej strony;
     jeśli w pliku jest gotowa warstwa `bufor_rzeki`, bierze ją),
  3. przelicza każdy DEM na wspólną siatkę (EPSG:32645, 30 m), poza obszarem daje -9999,
  4. wypełnia dziury w danych (interpolacja, a gdy nie wystarczy - wartości z Copernicusa),
  5. buduje pliki zgodne z kontraktem danych: elevation, release_height, reference_extent,
     entrainment_height - osobno dla każdego DEM-u, w data/nepal/derived/<wersja>/<dem>/.

Uruchomienie (z głównego folderu repo):
    python -m ravaflow_tools.prepare_nepal_data --version v1 --volume 140e6

Wymaga: numpy, rasterio, geopandas (shapely).
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

# ============================================================== KONFIGURACJA
ROOT = Path(__file__).resolve().parents[1]
NEPAL = ROOT / "data" / "nepal"

# Surowe pliki (te z prefiksem output_, sprzed jakiejkolwiek obróbki w QGIS)
DOWNLOADS = Path(r"C:\Users\huber\Downloads\cale_dane_nepal")
RAW = {
    "copernicus": DOWNLOADS / "rasters_COP30" / "output_hh.tif",
    "alos": DOWNLOADS / "rasters_AW3D30" / "output_AW3D30.tif",
    "srtm": DOWNLOADS / "rasters_SRTMGL1" / "output_SRTMGL1.tif",
}
DEMS = ("copernicus", "alos", "srtm")  # Copernicus pierwszy: baza dla zboczy i łatania dziur

# Warstwy wektorowe wyeksportowane z QGIS do jednego GeoPackage
VECTORS = NEPAL / "vector" / "nepal_vectors.gpkg"
LAYER_RIVER = "rzeka"          # linia: oś rzeki
LAYER_BUFFER = "bufor_rzeki"   # opcjonalnie: gotowy bufor (jeśli jest, ma pierwszeństwo)
LAYER_EXTENT = "zasieg"        # poligon: rzeczywisty zasięg fali ze zdjęć satelitarnych
LAYER_EROSION = None           # None = cały obszar roboczy (jak w pierwotnych danych); albo nazwa warstwy
SOURCE_LAYER = "zrodlo"        # w pliku source_<wersja>.gpkg

EPSG = 32645
CELL = 30.0           # m
BUFFER_M = 2000.0     # m z każdej strony osi rzeki
NODATA = -9999.0
MIN_ELEV = 100.0      # w tym rejonie nie ma terenu poniżej 100 m: niższe wartości to dziury/zera
FILL_DIST = 30        # maks. odległość interpolacji dziur [komórki]
MAX_HOLE_FRACTION = 0.02   # więcej dziur niż 2% obszaru = coś jest nie tak z surowym plikiem
H_MAX_WARN = 150.0    # m; ostrzeżenie, gdy masa źródłowa wychodzi grubsza
PATCH_OUTLIERS = True # łatanie artefaktów w ALOS/SRTM wartościami z Copernicusa
PATCH_THR = 200.0     # m: różnica względem Copernicusa uznawana za artefakt
PATCH_AGREE = 100.0   # m: drugi DEM musi zgadzać się z Copernicusem w tym zakresie
PATCH_GROW = 3        # komórek: powiększenie łaty (artefakt jest większy niż skrajne komórki)


# ===================================================== FUNKCJE LICZBOWE (numpy)
def slope_horn(z: np.ndarray, cell: float) -> np.ndarray:
    """Nachylenie [stopnie], metoda Horna (jak w QGIS). NaN tam, gdzie brakuje sąsiada."""
    p = np.pad(z.astype("float64"), 1, mode="constant", constant_values=np.nan)
    z1, z2, z3 = p[:-2, :-2], p[:-2, 1:-1], p[:-2, 2:]
    z4, z6 = p[1:-1, :-2], p[1:-1, 2:]
    z7, z8, z9 = p[2:, :-2], p[2:, 1:-1], p[2:, 2:]
    dzdx = ((z3 + 2 * z6 + z9) - (z1 + 2 * z4 + z7)) / (8 * cell)
    dzdy = ((z7 + 2 * z8 + z9) - (z1 + 2 * z2 + z3)) / (8 * cell)
    return np.degrees(np.arctan(np.hypot(dzdx, dzdy)))


def entrainment_classes(slope_deg: np.ndarray) -> np.ndarray:
    """Miąższość erodowalnego materiału [m]: >35 st. 0.5; 20-35 st. 1.5; <20 st. 3.0."""
    s = np.nan_to_num(slope_deg, nan=0.0)
    return np.where(s > 35, 0.5, np.where(s >= 20, 1.5, 3.0)).astype("float32")


def level_pool_release(dem: np.ndarray, mask: np.ndarray, volume: float, cell_area: float):
    """Masa uwolnienia jako pozioma tafla: h = max(0, Z - DEM) w poligonie.

    Poziom Z dobierany bisekcją tak, by objętość = volume [m3]. Wynik: (h, Z, objętość).
    DEM nie może mieć NaN w obrębie maski.
    """
    if not mask.any():
        raise ValueError("polygon źródła nie obejmuje żadnej komórki")
    dem = dem.astype("float64")
    z_lo = float(dem[mask].min())
    z_hi = float(dem[mask].max()) + 5000.0
    h, v, z = None, 0.0, z_hi
    for _ in range(200):
        z = 0.5 * (z_lo + z_hi)
        h = np.where(mask, np.clip(z - dem, 0.0, None), 0.0)
        v = float(h.sum()) * cell_area
        if abs(v - volume) <= 1e-6 * volume:
            break
        if v < volume:
            z_lo = z
        else:
            z_hi = z
    return h.astype("float32"), z, v


def shrink(mask: np.ndarray) -> np.ndarray:
    """Erozja maski o 1 komórkę (komórka zostaje tylko, gdy ma wszystkich 8 sąsiadów w masce)."""
    p = np.pad(mask, 1, constant_values=False)
    out = mask.copy()
    h, w = mask.shape
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            out &= p[1 + dy:1 + dy + h, 1 + dx:1 + dx + w]
    return out


def fill_holes(a: np.ndarray, domain: np.ndarray):
    """Wypełnia dziury (wartości < MIN_ELEV, w tym -9999) w obszarze `domain` interpolacją
    (do FILL_DIST komórek). Zwraca (tablica, maska dziur, maska dziur niewypełnionych)."""
    from rasterio.fill import fillnodata

    hole = domain & (a < MIN_ELEV)          # NODATA = -9999 też spełnia ten warunek
    valid = domain & ~hole
    work = np.where(valid, a, 0.0).astype("float32")
    filled = fillnodata(work, mask=valid.astype("uint8"),
                        max_search_distance=FILL_DIST, smoothing_iterations=0)
    left = hole & (filled < MIN_ELEV)       # interpolacja nie dosięgła
    return filled, hole, left


def grow(mask: np.ndarray, n: int) -> np.ndarray:
    """Powiększenie maski o n komórek (sąsiedztwo 8-komórkowe)."""
    out = mask.copy()
    h, w = mask.shape
    for _ in range(n):
        p = np.pad(out, 1, constant_values=False)
        nxt = out.copy()
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                nxt |= p[1 + dy:1 + dy + h, 1 + dx:1 + dx + w]
        out = nxt
    return out


def patch_outliers(dem: np.ndarray, other: np.ndarray, cop: np.ndarray, domain: np.ndarray):
    """Artefakty `dem`: różni się od Copernicusa o > PATCH_THR, a `other` zgadza się z Copernicusem
    (< PATCH_AGREE). Łata powiększona o PATCH_GROW komórek jest zastępowana Copernicusem.
    Zwraca (DEM, liczba flagowanych komórek, liczba załatanych komórek)."""
    flag = domain & (np.abs(dem - cop) > PATCH_THR) & (np.abs(other - cop) < PATCH_AGREE)
    patch = grow(flag, PATCH_GROW) & domain
    out = dem.copy()
    out[patch] = cop[patch]
    return out, int(flag.sum()), int(patch.sum())


def hybrid_release(dem: np.ndarray, mask: np.ndarray, volume: float, cell_area: float, base_m: float):
    """Masa uwolnienia = warstwa bazowa `base_m` w całym polygonie + pozioma tafla na tej bazie.

    Objętość tafli to volume - base_m * powierzchnia polygonu. Tafla wypełnia obniżenia terenu
    (DEM + baza) do wspólnego poziomu Z, więc tam, gdzie głębiej, jest więcej wody, a wyżej
    zostaje sama baza. Zwraca (h, Z, objętość, objętość bazy, objętość tafli).
    """
    if not mask.any():
        raise ValueError("polygon źródła nie obejmuje żadnej komórki")
    v_base = base_m * float(mask.sum()) * cell_area
    v_pool = volume - v_base
    if v_pool <= 0:
        raise ValueError(f"objętość {volume / 1e6:.1f} mln m3 jest nie większa niż sama warstwa bazowa "
                         f"({v_base / 1e6:.1f} mln m3 = {base_m:g} m na całym polygonie); "
                         "zwiększ --volume albo zmniejsz --base-m")
    base = np.where(mask, base_m, 0.0)
    pool, z, _ = level_pool_release(dem.astype("float64") + base, mask, v_pool, cell_area)
    h = base + pool.astype("float64")
    return h.astype("float32"), z, float(h.sum()) * cell_area, v_base, v_pool


# ============================================================ POMOCNICZE (I/O)
def read_layer(path: Path, layer: str, required: bool = True, polygon: bool = False):
    import geopandas as gpd

    try:
        gdf = gpd.read_file(path, layer=layer)
    except Exception as e:  # noqa: BLE001
        if required:
            raise SystemExit(f"Nie mogę wczytać warstwy '{layer}' z {path}:\n  {e}")
        return None
    if gdf.crs is None:
        raise SystemExit(f"Warstwa '{layer}' nie ma układu współrzędnych.")
    gdf = gdf.to_crs(epsg=EPSG)
    gdf = gdf[gdf.geometry.notna() & ~gdf.geometry.is_empty].copy()
    if polygon:
        gdf["geometry"] = gdf.geometry.buffer(0)   # naprawia drobne błędy geometrii
    if len(gdf) == 0:
        if required:
            raise SystemExit(f"Warstwa '{layer}' jest pusta.")
        return None
    return gdf


def burn(geoms, transform, shape) -> np.ndarray:
    from rasterio.features import rasterize

    return rasterize([(g, 1) for g in geoms], out_shape=shape, transform=transform,
                     fill=0, dtype="uint8").astype(bool)


def write_tif(path: Path, arr: np.ndarray, transform, nodata=None) -> None:
    import rasterio
    from rasterio.crs import CRS

    path.parent.mkdir(parents=True, exist_ok=True)
    with rasterio.open(path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1],
                       count=1, dtype="float32", crs=CRS.from_epsg(EPSG), transform=transform,
                       nodata=nodata, compress="deflate") as dst:
        dst.write(arr.astype("float32"), 1)


def read_tif(path: Path):
    import rasterio

    with rasterio.open(path) as f:
        return f.read(1).astype("float64"), f.transform


def warp_to_grid(raw_path: Path, transform, shape):
    """Surowy DEM -> wspólna siatka (EPSG, CELL), interpolacja dwuliniowa, brak danych = NODATA.

    Zwraca (DEM, maska pokrycia): pokrycie = komórki siatki, których środek leży w zasięgu surowego pliku."""
    import rasterio
    from rasterio.crs import CRS
    from rasterio.transform import array_bounds
    from rasterio.warp import Resampling, reproject, transform_bounds

    crs = CRS.from_epsg(EPSG)
    if not raw_path.exists():
        raise SystemExit(f"Brak pliku: {raw_path}")
    with rasterio.open(raw_path) as src:
        print(f"  surowy: {raw_path.name} | CRS {src.crs} | {src.width}x{src.height} | "
              f"nodata {src.nodata}")
        a = src.read(1).astype("float32")
        if src.nodata is not None:
            a[a == src.nodata] = NODATA
        a[a < MIN_ELEV] = NODATA                  # zera, -9999, -32768: nie mogą rozmywać sąsiadów
        l, b, r, t = transform_bounds(crs, src.crs, *array_bounds(shape[0], shape[1], transform))
        sl, sb, sr, st = src.bounds
        if l < sl or b < sb or r > sr or t > st:
            print("  UWAGA: surowy plik może nie pokrywać całego obszaru roboczego "
                  "(sprawdzimy to po liczbie dziur)")
        dst = np.full(shape, NODATA, dtype="float32")
        reproject(source=a, destination=dst, src_transform=src.transform, src_crs=src.crs,
                  src_nodata=NODATA, dst_transform=transform, dst_crs=crs, dst_nodata=NODATA,
                  resampling=Resampling.bilinear)
        # pokrycie: przeliczamy "jedynki" metodą najbliższego sąsiada
        cov = np.zeros(shape, dtype="uint8")
        reproject(source=np.ones((src.height, src.width), dtype="uint8"), destination=cov,
                  src_transform=src.transform, src_crs=src.crs, dst_transform=transform,
                  dst_crs=crs, dst_nodata=0, resampling=Resampling.nearest)
    return dst, cov.astype(bool)


# ================================================================ GŁÓWNY PRZEBIEG
def main() -> None:
    import geopandas as gpd
    from rasterio.transform import from_origin

    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("--version", required=True, help="nazwa wersji źródła, np. v1")
    ap.add_argument("--volume", required=True, type=float, help="objętość masy źródłowej [m3], np. 140e6")
    ap.add_argument("--release", choices=("hybrid", "level"), default="hybrid",
                    help="hybrid = warstwa bazowa + pozioma tafla (domyślnie); level = sama pozioma tafla")
    ap.add_argument("--base-m", type=float, default=50.0,
                    help="dla --release hybrid: grubość warstwy bazowej w całym polygonie [m]")
    ap.add_argument("--source-gpkg", type=Path, default=None,
                    help="domyślnie data/nepal/source/source_<wersja>.gpkg")
    ap.add_argument("--reuse-dem", action="store_true",
                    help="nie buduj DEM-ów od nowa, użyj istniejących dem_original.tif i domain_mask.tif")
    args = ap.parse_args()
    if args.volume <= 0:
        raise SystemExit("--volume musi być dodatnie")

    out_root = NEPAL / "derived" / args.version
    if out_root.exists():
        raise SystemExit(f"{out_root} już istnieje - nie nadpisuję. Użyj nowej --version.")
    source_gpkg = args.source_gpkg or (NEPAL / "source" / f"source_{args.version}.gpkg")

    # ---------- wektory
    print("[1] Warstwy wektorowe")
    extent_gdf = read_layer(VECTORS, LAYER_EXTENT, polygon=True)
    erosion_gdf = read_layer(VECTORS, LAYER_EROSION, polygon=True) if LAYER_EROSION else None
    source_gdf = read_layer(source_gpkg, SOURCE_LAYER, polygon=True)

    buf = read_layer(VECTORS, LAYER_BUFFER, required=False, polygon=True)
    if buf is not None:
        domain_geoms = list(buf.geometry)
        domain_note = f"gotowa warstwa '{LAYER_BUFFER}'"
    else:
        river = read_layer(VECTORS, LAYER_RIVER)
        domain_geoms = list(river.geometry.buffer(BUFFER_M))
        domain_note = f"bufor {BUFFER_M:.0f} m wokół warstwy '{LAYER_RIVER}'"
    print(f"  obszar roboczy: {domain_note}")

    # ---------- siatka i obszar roboczy
    if args.reuse_dem:
        _, transform = read_tif(NEPAL / "copernicus" / "dem_original.tif")
        dom_arr, _ = read_tif(NEPAL / "domain_mask.tif")
        domain = dom_arr == 1
        shape = domain.shape
    else:
        minx, miny, maxx, maxy = gpd.GeoSeries(domain_geoms, crs=EPSG).total_bounds
        x0 = np.floor(minx / CELL) * CELL
        y1 = np.ceil(maxy / CELL) * CELL
        width = int(np.ceil((maxx - x0) / CELL))
        height = int(np.ceil((y1 - miny) / CELL))
        transform = from_origin(x0, y1, CELL, CELL)
        shape = (height, width)
        domain = burn(domain_geoms, transform, shape)
    print(f"  siatka {shape[1]}x{shape[0]}, {CELL:.0f} m, komórek w obszarze (bufor): {int(domain.sum())}")

    # ---------- DEM-y
    dem_res, dem_stats = {}, {}
    if args.reuse_dem:
        print("[2] DEM-y: używam istniejących dem_original.tif")
        for name in DEMS:
            d, _ = read_tif(NEPAL / name / "dem_original.tif")
            dem_res[name] = np.where(domain, d, NODATA).astype("float32")
            dem_stats[name] = {"holes": "-", "from_interp": "-", "from_copernicus": "-", "patched": "-"}
    else:
        print("[2] DEM-y")
        warped, cover = {}, np.ones(shape, dtype=bool)
        for name in DEMS:
            print(f" {name}")
            warped[name], cov = warp_to_grid(RAW[name], transform, shape)
            print(f"  komórek bufora poza zasięgiem surowego pliku: {int((domain & ~cov).sum())}")
            cover &= cov
        # obszar roboczy = bufor ∩ pokrycie wszystkich trzech DEM-ów (z marginesem 1 komórki)
        n_before = int(domain.sum())
        domain = domain & shrink(cover)
        print(f" obszar roboczy po przycięciu do pokrycia DEM-ów: {int(domain.sum())} komórek "
              f"(odcięto {n_before - int(domain.sum())})")
        write_tif(NEPAL / "domain_mask.tif", domain.astype("float32"), transform)

        for name in DEMS:
            print(f" {name}")
            a = warped[name]
            filled, hole, left = fill_holes(a, domain)
            print(f"  dziur w obszarze (wartości < {MIN_ELEV:.0f} m): {int(hole.sum())} "
                  f"({100 * hole.sum() / domain.sum():.2f}%), po interpolacji zostaje {int(left.sum())}")
            dbg = np.zeros(shape, dtype="float32")
            dbg[hole] = 1.0                   # 1 = dziura wypełniona interpolacją
            dbg[left] = 2.0                   # 2 = dziura niewypełniona
            write_tif(NEPAL / "debug" / f"{name}_dziury.tif", dbg, transform)
            if hole.sum() > MAX_HOLE_FRACTION * domain.sum():
                raise SystemExit("  Za dużo dziur w obszarze pokrytym przez surowy plik - sprawdź "
                                 f"data/nepal/debug/{name}_dziury.tif w QGIS.")
            n_cop = 0
            if left.any():
                rows, cols = np.nonzero(left)
                print(f"  niewypełnione: x {transform.c + cols.min() * CELL:.0f}-"
                      f"{transform.c + (cols.max() + 1) * CELL:.0f}, "
                      f"y {transform.f - (rows.max() + 1) * CELL:.0f}-{transform.f - rows.min() * CELL:.0f}")
                if name == "copernicus":
                    raise SystemExit("  Copernicus ma dziury, których nie da się wypełnić interpolacją "
                                     f"(maska: data/nepal/debug/{name}_dziury.tif).")
                filled[left] = dem_res["copernicus"][left]
                n_cop = int(left.sum())
            st = {"holes": int(hole.sum()), "from_interp": int((hole & ~left).sum()),
                  "from_copernicus": n_cop}
            out = np.where(domain, filled, NODATA).astype("float32")
            st["patched"] = 0
            dem_res[name], dem_stats[name] = out, st
            print(f"  wypełniono: interpolacją {st['from_interp']}, z Copernicusa "
                  f"{st['from_copernicus']} | min {out[domain].min():.0f} m, "
                  f"max {out[domain].max():.0f} m")

        if PATCH_OUTLIERS:
            print(" łatanie artefaktów (DEM odbiega od Copernicusa, a drugi DEM się z nim zgadza)")
            cop_f = dem_res["copernicus"]
            before = {n: dem_res[n].copy() for n in ("alos", "srtm")}
            for name, other in (("alos", "srtm"), ("srtm", "alos")):
                new, n_flag, n_patch = patch_outliers(before[name], before[other], cop_f, domain)
                dem_res[name] = new
                dem_stats[name]["patched"] = n_patch
                print(f"  {name}: artefaktów {n_flag} komórek, załatano {n_patch} komórek (Copernicus)")

        for name in DEMS:
            p = NEPAL / name / "dem_original.tif"
            if p.exists() and not (NEPAL / name / "dem_original_stary.tif").exists():
                shutil.copy(p, NEPAL / name / "dem_original_stary.tif")
            write_tif(p, dem_res[name], transform, nodata=NODATA)

    # ---------- rastry wspólne dla wszystkich DEM-ów
    print("[3] Zasięg, porywanie, źródło")
    extent = burn(list(extent_gdf.geometry), transform, shape)
    erosion = burn(list(erosion_gdf.geometry), transform, shape) if erosion_gdf is not None else domain.copy()
    smask = burn(list(source_gdf.geometry), transform, shape)

    if (extent & ~domain).any():
        print(f"  UWAGA: {int((extent & ~domain).sum())} komórek zasięgu leży poza obszarem roboczym")
    if not domain[smask].all():
        raise SystemExit("Polygon źródła wychodzi poza obszar roboczy - popraw go w QGIS.")
    if smask.sum() < 50:
        print(f"  UWAGA: polygon źródła ma tylko {int(smask.sum())} komórek - siatka 30 m słabo go odwzoruje")

    # reference_extent: 1 = zasięg, 0 = brak zasięgu, -1 = brak informacji (poza obszarem)
    reference_extent = np.where(domain, extent.astype("float32"), -1.0).astype("float32")

    cop_nan = np.where(domain, dem_res["copernicus"], np.nan)
    slope = slope_horn(cop_nan, CELL)
    entrainment = np.where(erosion & domain, entrainment_classes(slope), 0.0).astype("float32")
    print(f"  komórek: zasięg {int(extent.sum())}, strefa erozji {int(erosion.sum())}, "
          f"źródło {int(smask.sum())}")
    print(f"  maks. objętość do porwania: "
          f"{float(entrainment.sum()) * CELL * CELL / 1e6:.1f} mln m3")

    # ---------- pliki kontraktu dla każdego DEM-u
    print("[4] Pliki kontraktu")
    polygon_hash = hashlib.sha256(source_gpkg.read_bytes()).hexdigest()
    cell_area = CELL * CELL
    summary = {}
    for name in DEMS:
        dem = dem_res[name].astype("float64")
        elev = np.where(domain, dem, 0.0)                      # nie odejmujemy masy (patrz NOTES)
        v_base = v_pool = None
        if args.release == "level":
            h, z_level, vol = level_pool_release(elev, smask, args.volume, cell_area)
        else:
            try:
                h, z_level, vol, v_base, v_pool = hybrid_release(elev, smask, args.volume, cell_area,
                                                                 args.base_m)
            except ValueError as e:
                raise SystemExit(f"  {e}")
        h_max = float(h.max())
        n_dry = int((smask & (h <= 0)).sum())
        d = out_root / name
        write_tif(d / "elevation.tif", elev, transform)
        write_tif(d / "release_height.tif", h, transform)
        write_tif(d / "reference_extent.tif", reference_extent, transform)
        write_tif(d / "entrainment_height.tif", entrainment, transform)
        manifest = {
            "version": args.version, "dem": name, "raw_file": str(RAW[name]),
            "created": dt.datetime.now().isoformat(timespec="seconds"),
            "target_volume_m3": args.volume, "actual_volume_m3": vol,
            "release_mode": args.release,
            "base_m": args.base_m if args.release == "hybrid" else None,
            "base_volume_m3": v_base, "pool_volume_m3": v_pool, "dry_source_cells": n_dry,
            "z_level_m": z_level, "h_max_m": h_max, "source_cells": int(smask.sum()),
            "source_polygon_sha256": polygon_hash, "source_file": str(source_gpkg),
            "cell_m": CELL, "epsg": EPSG, "domain": domain_note,
        }
        (d / "manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False),
                                         encoding="utf-8")
        summary[name] = manifest
        flag = "  <-- UWAGA: bardzo gruba masa, rozszerz polygon w górę doliny" if h_max > H_MAX_WARN else ""
        split = (f", baza {v_base / 1e6:.1f} + tafla {v_pool / 1e6:.1f} mln m3" if v_base is not None else "")
        print(f"  {name}: Z = {z_level:.1f} m n.p.m., h_max = {h_max:.1f} m, objętość = "
              f"{vol / 1e6:.3f} mln m3{split}; komórek polygonu bez masy: {n_dry}{flag}")

    # ---------- kontrola zgodności DEM-ów
    print("[5] Porównanie DEM-ów (względem Copernicusa, w obszarze roboczym)")
    cop = dem_res["copernicus"].astype("float64")
    for name in ("alos", "srtm"):
        diff = dem_res[name].astype("float64") - cop
        big = domain & (np.abs(diff) > 200)
        dd = diff[domain]
        print(f"  {name}: mediana {np.median(dd):.1f} m, 95% |różnicy| "
              f"{np.percentile(np.abs(dd), 95):.1f} m, >200 m: {int(big.sum())}, "
              f"w zasięgu fali: {int((big & extent).sum())}")

    # ---------- NOTES.md
    rows = "\n".join(
        f"| {n} | {dem_stats[n]['holes']} | {dem_stats[n]['from_interp']} | "
        f"{dem_stats[n]['from_copernicus']} | {dem_stats[n]['patched']} | "
        f"{summary[n]['h_max_m']:.1f} | {summary[n]['actual_volume_m3'] / 1e6:.3f} |"
        for n in DEMS)
    if args.release == "level":
        release_txt = "pozioma tafla w polygonie źródła, h = max(0, Z - DEM); poziom Z dobrany bisekcją."
    else:
        release_txt = (f"warstwa bazowa {args.base_m:g} m w całym polygonie źródła + pozioma tafla (poziom Z) "
                       "wypełniająca obniżenia terenu (DEM + baza) resztą objętości.")
    patch_txt = (f"Artefakty ALOS/SRTM (różnica od Copernicusa > {PATCH_THR:.0f} m, gdy drugi DEM zgadza się "
                 f"z Copernicusem w {PATCH_AGREE:.0f} m) zastąpiono Copernicusem, łata powiększona o {PATCH_GROW} komórki."
                 if PATCH_OUTLIERS else "Artefaktów nie łatano.")
    notes = f"""# Dane wejściowe, wersja źródła `{args.version}`

Wygenerowano: {dt.date.today().isoformat()} skryptem `ravaflow_tools/prepare_nepal_data.py`.

## Źródła
- DEM-y (surowe pliki): {", ".join(f"{n}: `{RAW[n].name}`" for n in DEMS)}
- Polygon źródła: `{source_gpkg.name}`, warstwa `{SOURCE_LAYER}` (SHA-256 {polygon_hash[:12]}...)
- Zasięg rzeczywisty: warstwa `{LAYER_EXTENT}`; strefa erozji: {'warstwa `' + LAYER_EROSION + '`' if LAYER_EROSION else 'cały obszar roboczy'} (plik `{VECTORS.name}`)
- Obszar roboczy: {domain_note}

## Siatka
GeoTIFF, EPSG:{EPSG}, komórka {CELL:.0f} m, {shape[1]}x{shape[0]}, identyczna dla wszystkich plików i DEM-ów.
Poza obszarem roboczym: `elevation` = 0, `release_height` = 0, `entrainment_height` = 0,
`reference_extent` = -1 (brak informacji).

## Jak powstały pliki
- `elevation.tif`: DEM sprzed zdarzenia, przeliczony na siatkę (interpolacja dwuliniowa), dziury
  (wartości < {MIN_ELEV:.0f} m) wypełnione interpolacją (do {FILL_DIST} komórek), a gdy nie wystarczyła,
  wartościami z Copernicusa. {patch_txt}
- `release_height.tif`: {release_txt} Objętość = {args.volume:.4g} m3. Poza polygonem dokładnie 0.
- `reference_extent.tif`: 1 w zasięgu, 0 w obszarze roboczym poza zasięgiem, -1 poza obszarem roboczym.
- `entrainment_height.tif`: z nachylenia Copernicusa (metoda Horna): >35 st. 0.5 m, 20-35 st. 1.5 m,
  <20 st. 3.0 m; tylko w strefie erozji, poza nią 0. Wartości miąższości są założeniem empirycznym.

## Odstępstwo od opisu kontraktu
Kontrakt mówi, że `elevation + release_height` = teren sprzed zdarzenia (źródło w bliźnie oderwania).
Tu źródło jest w dnie doliny (masa nadeszła z góry), więc `elevation` = teren sprzed zdarzenia bez
odejmowania masy, a `release_height` leży na nim. Odejmowanie wykopałoby w dolinie rów.

## Statystyki
| DEM | dziur | interpolacja | z Copernicusa | załatane artefakty | h_max [m] | objętość [mln m3] |
|---|---|---|---|---|---|---|
{rows}
"""
    (out_root / "NOTES.md").write_text(notes, encoding="utf-8")

    print("\nGotowe. Kontrola kontraktu dla każdego DEM-u:")
    for name in DEMS:
        print(f"  python -m ravaflow_tools.check_inputs data/nepal/derived/{args.version}/{name} "
              f"--volume {args.volume:g} --cellsize {CELL:g} --epsg {EPSG}")


if __name__ == "__main__":
    main()
"""Generator syntetycznego zestawu danych w formacie kontraktu danych.

Zestaw nie ma sensu fizycznego dla Nepalu - służy do testu całego łańcucha
(walidator -> r.avaflow -> metryki -> wizualizacja -> scenariusze) przed
dostarczeniem prawdziwych danych. Współrzędne udają UTM 45N (okolice Langtang),
żeby sprawdzić, że duże wartości współrzędnych nie psują importu do GRASS.

Użycie:
    python -m ravaflow_tools.synthetic data/synthetic                 # GeoTIFF
    python -m ravaflow_tools.synthetic data/synthetic --format asc    # bez rasterio
Pliki: elevation, release_height, reference_extent, entrainment_height.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .metrics import write_raster

NODATA = -9999.0


def generate(out_dir: str | Path, fmt: str = "tif", volume: float = 2.0e6,
             cellsize: float = 30.0, epsg: int = 32645, seed: int = 0,
             xll: float = 340000.0, yll: float = 3110000.0,
             nrows: int = 300, ncols: int = 200) -> dict:
    """Tworzy zestaw w out_dir. Zwraca meta siatki i słownik ścieżek."""
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    meta = {"ncols": ncols, "nrows": nrows, "xll": xll, "yll": yll, "cellsize": cellsize, "nodata": NODATA}

    rows, cols = np.mgrid[0:nrows, 0:ncols]
    s = rows * cellsize                       # odległość od góry (północ) wzdłuż doliny [m]
    x = (cols + 0.5) * cellsize               # współrzędna poprzeczna [m]
    xc = ncols * cellsize / 2 + 400.0 * np.sin(s / 1500.0)   # oś doliny (lekko meandruje)
    d = np.abs(x - xc)                        # odległość od osi doliny

    # profil podłużny: stromo - umiarkowanie - płasko (wybieg)
    slope = np.where(s < 3000, 0.45, np.where(s < 6000, 0.25, 0.06))
    z_long = 5500.0 - np.cumsum(slope * cellsize, axis=0) + slope[0, 0] * cellsize
    rng = np.random.default_rng(seed)
    noise = rng.normal(0, 1.0, size=(nrows, ncols))
    from scipy.ndimage import gaussian_filter
    noise = gaussian_filter(noise, 3) * 6.0   # łagodne zaburzenia terenu
    terrain = z_long + np.minimum(0.30 * d, 700.0) + noise   # teren przed zdarzeniem

    # masa: eliptyczny blok u góry doliny, grubość paraboliczna, normowana do zadanej objętości
    s0, s_len, half_w = 450.0, 500.0, 260.0
    inside = ((s - s0) / (s_len / 2)) ** 2 + (d / half_w) ** 2
    shape = np.clip(1.0 - inside, 0.0, None)
    shape = np.where(shape < 0.01, 0.0, shape)    # bez śladowych grubości na brzegu bloku
    h = shape * volume / (shape.sum() * cellsize ** 2)
    elevation = terrain - h                   # powierzchnia PODSTAWY masy

    entr = np.where(d < 150.0, 2.0, 0.0)      # opcjonalna warstwa do unoszenia

    ref = ((d < 220.0) & (s > 150.0) & (s < 5200.0)).astype(float)
    ref[:10, :10] = -1.0                      # fragment "brak informacji" (test konwencji)

    ext = "tif" if fmt == "tif" else "asc"
    paths = {}
    for name, arr in (("elevation", elevation), ("release_height", h),
                      ("reference_extent", ref), ("entrainment_height", entr)):
        p = out / f"{name}.{ext}"
        write_raster(p, arr, meta, epsg=epsg if ext == "tif" else None)
        paths[name] = p
    return {"meta": meta, "paths": paths, "volume": float(h.sum() * cellsize ** 2)}


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("out_dir")
    ap.add_argument("--format", choices=["tif", "asc"], default="tif")
    ap.add_argument("--volume", type=float, default=2.0e6, help="objętość masy [m3]")
    ap.add_argument("--cellsize", type=float, default=30.0)
    a = ap.parse_args()
    r = generate(a.out_dir, a.format, a.volume, a.cellsize)
    print(f"Zapisano {len(r['paths'])} plików w {a.out_dir}; objętość {r['volume']:.3e} m3")

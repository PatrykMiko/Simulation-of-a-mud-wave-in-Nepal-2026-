"""Metryki zgodności zasięgu - niezależne od modelu (to samo dla r.avaflow i SS3).

Siatka rastra opisana słownikiem meta:
    {"ncols", "nrows", "xll", "yll", "cellsize", "nodata"}
Wiersz 0 = północ (jak w ASCII i GeoTIFF north-up). Komórki zakładamy kwadratowe.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
from scipy import ndimage


# ---------------------------------------------------------------- wczytywanie
def read_ascii(path: str | Path) -> tuple[np.ndarray, dict]:
    """Czyta ESRI ASCII grid (.asc). Zwraca (tablica float, meta)."""
    header: dict[str, float] = {}
    with open(path, encoding="utf-8") as f:
        for _ in range(6):
            pos = f.tell()
            parts = f.readline().split()
            if len(parts) != 2 or not parts[0][0].isalpha():
                f.seek(pos)  # brak nodata w nagłówku
                break
            header[parts[0].lower()] = float(parts[1])
        data = np.loadtxt(f, dtype=float, ndmin=2)
    ncols, nrows = int(header["ncols"]), int(header["nrows"])
    cell = header["cellsize"]
    xll = header.get("xllcorner", header.get("xllcenter", 0.0) - cell / 2)
    yll = header.get("yllcorner", header.get("yllcenter", 0.0) - cell / 2)
    nodata = header.get("nodata_value", -9999.0)
    if data.shape != (nrows, ncols):
        raise ValueError(f"{path}: rozmiar danych {data.shape} != nagłówek {(nrows, ncols)}")
    meta = {"ncols": ncols, "nrows": nrows, "xll": xll, "yll": yll, "cellsize": cell, "nodata": nodata}
    return data, meta


def read_raster(path: str | Path) -> tuple[np.ndarray, dict]:
    """Czyta .asc (własny parser) lub .tif (rasterio)."""
    path = Path(path)
    if path.suffix.lower() == ".asc":
        return read_ascii(path)
    import rasterio  # tylko dla GeoTIFF

    with rasterio.open(path) as ds:
        data = ds.read(1).astype(float)
        t = ds.transform
        if t.b != 0 or t.d != 0 or t.e > 0 or abs(abs(t.a) - abs(t.e)) > 1e-9:
            raise ValueError(f"{path}: wymagany raster north-up z kwadratowymi komórkami")
        meta = {
            "ncols": ds.width, "nrows": ds.height, "xll": t.c,
            "yll": t.f + t.e * ds.height, "cellsize": t.a,
            "nodata": ds.nodata if ds.nodata is not None else -9999.0,
        }
    return data, meta


def write_ascii(path: str | Path, arr: np.ndarray, meta: dict) -> None:
    """Zapisuje ESRI ASCII grid (.asc)."""
    nd = meta.get("nodata", -9999.0)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"ncols {arr.shape[1]}\nnrows {arr.shape[0]}\n"
                f"xllcorner {meta['xll']}\nyllcorner {meta['yll']}\n"
                f"cellsize {meta['cellsize']}\nNODATA_value {nd}\n")
        np.savetxt(f, arr, fmt="%.4f")


def write_raster(path: str | Path, arr: np.ndarray, meta: dict, epsg: int | None = None) -> None:
    """Zapisuje .asc (bez rasterio) albo GeoTIFF (rasterio) w zależności od rozszerzenia."""
    path = Path(path)
    if path.suffix.lower() == ".asc":
        write_ascii(path, arr, meta)
        return
    import rasterio
    from rasterio.transform import from_origin

    top = meta["yll"] + meta["nrows"] * meta["cellsize"]
    kwargs = {}
    if epsg is not None:
        kwargs["crs"] = f"EPSG:{epsg}"
    with rasterio.open(
        path, "w", driver="GTiff", height=arr.shape[0], width=arr.shape[1], count=1,
        dtype="float32", nodata=meta.get("nodata", -9999.0),
        transform=from_origin(meta["xll"], top, meta["cellsize"], meta["cellsize"]), **kwargs,
    ) as ds:
        ds.write(arr.astype("float32"), 1)


def sample_to_grid(arr: np.ndarray, src: dict, dst: dict) -> np.ndarray:
    """Próbkowanie najbliższym sąsiadem: arr (siatka src) -> siatka dst.

    Komórki dst leżące poza src dostają wartość nodata src.
    """
    cols = np.arange(dst["ncols"])
    rows = np.arange(dst["nrows"])
    x = dst["xll"] + (cols + 0.5) * dst["cellsize"]
    y = dst["yll"] + (dst["nrows"] - rows - 0.5) * dst["cellsize"]
    ci = np.floor((x - src["xll"]) / src["cellsize"]).astype(int)
    ri = np.floor((src["yll"] + src["nrows"] * src["cellsize"] - y) / src["cellsize"]).astype(int)
    out = np.full((dst["nrows"], dst["ncols"]), src["nodata"], dtype=float)
    rmask = (ri >= 0) & (ri < src["nrows"])
    cmask = (ci >= 0) & (ci < src["ncols"])
    out[np.ix_(rmask, cmask)] = arr[np.ix_(ri[rmask], ci[cmask])]
    return out


def grid_differences(a: dict, b: dict, tol: float = 1e-6) -> list[str]:
    """Czytelna lista różnic między siatkami a i b (np. model i referencja).

    Pusta lista oznacza identyczne siatki. Opisuje dokładnie te cechy, które porównuje
    same_grid: rozmiar, rozmiar komórki i położenie lewego dolnego rogu.
    """
    diffs = []
    if a["ncols"] != b["ncols"] or a["nrows"] != b["nrows"]:
        diffs.append(f"rozmiar (kolumny x wiersze): {a['ncols']}x{a['nrows']} vs {b['ncols']}x{b['nrows']}")
    same_cell = abs(a["cellsize"] - b["cellsize"]) < tol
    if not same_cell:
        diffs.append(f"rozmiar komórki: {a['cellsize']} m vs {b['cellsize']} m")
    if abs(a["xll"] - b["xll"]) >= tol or abs(a["yll"] - b["yll"]) >= tol:
        text = f"róg siatki: ({a['xll']}, {a['yll']}) vs ({b['xll']}, {b['yll']})"
        if same_cell:  # przesunięcie o ułamek komórki (z pominięciem całych komórek)
            cs = a["cellsize"]
            dx = (a["xll"] - b["xll"] + cs / 2) % cs - cs / 2
            dy = (a["yll"] - b["yll"] + cs / 2) % cs - cs / 2
            text += f"; przesunięcie o ułamek komórki: dx={dx:+.1f} m, dy={dy:+.1f} m"
        diffs.append(text)
    return diffs


def same_grid(a: dict, b: dict, tol: float = 1e-6) -> bool:
    return not grid_differences(a, b, tol)

# ------------------------------------------------------------- maski zasięgu
def model_extent(hmax: np.ndarray, meta: dict, threshold: float = 0.1) -> np.ndarray:
    """Zasięg modelu = maksymalna grubość przepływu >= próg (domyślnie 0.1 m,
    jak pierwsza wartość parametru `thresholds` w r.avaflow)."""
    valid = hmax != meta["nodata"]
    return valid & (hmax >= threshold)


def reference_masks(ref: np.ndarray, meta: dict) -> tuple[np.ndarray, np.ndarray]:
    """Zasięg rzeczywisty wg konwencji r.avaflow: >0 = zasięg, 0 = brak,
    wartości ujemne / nodata = brak informacji. Zwraca (obserwowany, ważny)."""
    valid = (ref != meta["nodata"]) & (ref >= 0)
    return valid & (ref > 0), valid


# ------------------------------------------------------------------- metryki
def confusion(model: np.ndarray, obs: np.ndarray, valid: np.ndarray | None = None) -> dict:
    if valid is None:
        valid = np.ones_like(obs, dtype=bool)
    m, o = model & valid, obs & valid
    tp = int(np.sum(m & o))
    fp = int(np.sum(m & ~o))
    fn = int(np.sum(~m & o & valid))
    tn = int(np.sum(~m & ~o & valid))
    return {"tp": tp, "fp": fp, "fn": fn, "tn": tn}


def hausdorff_m(a: np.ndarray, b: np.ndarray, cellsize: float) -> float:
    """Symetryczna odległość Hausdorffa między zbiorami komórek (w metrach)."""
    if not a.any() or not b.any():
        return float("nan")
    d_to_a = ndimage.distance_transform_edt(~a)
    d_to_b = ndimage.distance_transform_edt(~b)
    return float(max(d_to_a[b].max(), d_to_b[a].max()) * cellsize)


def score(model: np.ndarray, obs: np.ndarray, valid: np.ndarray, cellsize: float) -> dict:
    """Metryki z planu: Jaccard e (iloraz, jak w SS3 i CSI w r.avaflow),
    precision, recall, F1, stosunek pól (= FoC w r.avaflow), różnica pól, Hausdorff."""
    c = confusion(model, obs, valid)
    tp, fp, fn = c["tp"], c["fp"], c["fn"]
    union = tp + fp + fn
    precision = tp / (tp + fp) if tp + fp else float("nan")
    recall = tp / (tp + fn) if tp + fn else float("nan")
    f1 = 2 * tp / (2 * tp + fp + fn) if (2 * tp + fp + fn) else float("nan")
    cell_area = cellsize ** 2
    area_model = (tp + fp) * cell_area
    area_obs = (tp + fn) * cell_area
    return {
        **c,
        "jaccard": tp / union if union else float("nan"),
        "precision": precision,
        "recall": recall,
        "f1": f1,
        "area_model_m2": area_model,
        "area_obs_m2": area_obs,
        "area_ratio": area_model / area_obs if area_obs else float("nan"),
        "area_diff_m2": area_model - area_obs,
        "hausdorff_m": hausdorff_m(model & valid, obs & valid, cellsize),
    }

"""Wizualizacje wyników: porównanie zasięgu, maksymalna grubość, animacja."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, ListedColormap
from matplotlib.patches import Patch


def hillshade(elev: np.ndarray, cellsize: float, nodata: float = -9999.0,
              azimuth: float = 315.0, altitude: float = 45.0) -> np.ndarray:
    """Prosty hillshade 0-1 (do tła map)."""
    z = np.where(elev == nodata, np.nan, elev)
    z = np.where(np.isnan(z), np.nanmean(z), z)
    dzdy, dzdx = np.gradient(z, cellsize)
    slope = np.arctan(np.hypot(dzdx, dzdy))
    aspect = np.arctan2(dzdy, -dzdx)
    az, alt = np.radians(360.0 - azimuth + 90.0), np.radians(altitude)
    hs = np.sin(alt) * np.cos(slope) + np.cos(alt) * np.sin(slope) * np.cos(az - aspect)
    return np.clip(hs, 0, 1)


def _extent(meta: dict) -> list[float]:
    cs = meta["cellsize"]
    return [meta["xll"], meta["xll"] + meta["ncols"] * cs,
            meta["yll"], meta["yll"] + meta["nrows"] * cs]


def plot_extent_comparison(elev, model, obs, valid, meta, out_png, title="", metrics=None):
    """Mapa: trafienia (TP), nadmiar modelu (FP), niedobór (FN) na hillshade."""
    fig, ax = plt.subplots(figsize=(9, 7), dpi=130)
    ext = _extent(meta)
    ax.imshow(hillshade(elev, meta["cellsize"], meta["nodata"]), cmap="gray",
              extent=ext, vmin=0, vmax=1.2)
    cat = np.zeros(model.shape, dtype=float)
    cat[model & obs & valid] = 1   # TP
    cat[model & ~obs & valid] = 2  # FP
    cat[~model & obs & valid] = 3  # FN
    layer = np.ma.masked_where(cat == 0, cat)
    cmap = ListedColormap(["#2e9d4f", "#e8590c", "#1c7ed6"])
    ax.imshow(layer, cmap=cmap, vmin=1, vmax=3, extent=ext, alpha=0.65, interpolation="nearest")
    ax.legend(handles=[Patch(color="#2e9d4f", label="zgodne (TP)"),
                       Patch(color="#e8590c", label="nadmiar modelu (FP)"),
                       Patch(color="#1c7ed6", label="niedobór modelu (FN)")],
              loc="upper right", fontsize=8)
    if metrics:
        txt = (f"Jaccard e = {metrics['jaccard']:.3f}\nprecision = {metrics['precision']:.3f}\n"
               f"recall = {metrics['recall']:.3f}\nF1 = {metrics['f1']:.3f}\n"
               f"pole model/obs = {metrics['area_ratio']:.2f}")
        ax.text(0.01, 0.01, txt, transform=ax.transAxes, fontsize=8, va="bottom",
                bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_hmax(elev, hmax, meta, out_png, title="", threshold=0.1, vmax=None):
    """Mapa maksymalnej grubości przepływu (skala logarytmiczna) na hillshade."""
    fig, ax = plt.subplots(figsize=(9, 7), dpi=130)
    ext = _extent(meta)
    ax.imshow(hillshade(elev, meta["cellsize"], meta["nodata"]), cmap="gray",
              extent=ext, vmin=0, vmax=1.2)
    h = np.ma.masked_where((hmax == meta["nodata"]) | (hmax < threshold), hmax)
    im = ax.imshow(h, cmap="turbo", norm=LogNorm(vmin=threshold, vmax=vmax or max(h.max(), threshold * 10)),
                   extent=ext, alpha=0.8, interpolation="nearest")
    fig.colorbar(im, ax=ax, label="maks. grubość przepływu [m]", shrink=0.8)
    ax.set_title(title)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def make_gif(frame_paths: list[str | Path], out_gif: str | Path, fps: int = 5):
    """Składa klatki PNG w animację GIF."""
    import imageio.v2 as imageio

    frames = [imageio.imread(p) for p in frame_paths]
    imageio.mimsave(out_gif, frames, duration=1.0 / fps, loop=0)


def animate_hflow(ascii_dir, prefix, elev, meta, out_gif, threshold=0.1, every=1, fps=5):
    """Animacja grubości przepływu z klatek <prefix>_hflowNNNN.asc."""
    from .metrics import read_ascii
    import tempfile

    files = sorted(Path(ascii_dir).glob(f"{prefix}_hflow[0-9][0-9][0-9][0-9].asc"))[::every]
    if not files:
        raise FileNotFoundError(f"Brak klatek {prefix}_hflowNNNN.asc w {ascii_dir}")
    arrays = [read_ascii(f)[0] for f in files]
    vmax = max(a.max() for a in arrays)
    with tempfile.TemporaryDirectory() as tmp:
        pngs = []
        for i, a in enumerate(arrays):
            p = Path(tmp) / f"f{i:04d}.png"
            plot_hmax(elev, a, meta, p, title=f"{prefix}: krok {i * every}", threshold=threshold, vmax=vmax)
            pngs.append(p)
        make_gif(pngs, out_gif, fps=fps)

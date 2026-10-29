#!/usr/bin/env python3
"""Headless diagnostics and visualization."""

from __future__ import annotations

import numpy as np

from .automaton import SciddicaSS3
from .constants import SQRT3


def run_headless(
    model: SciddicaSS3,
    steps: int = 300,
    report_every: int = 50,
    progress: bool = True,
) -> None:
    """Run the model and report debris, soil and volume conservation."""
    if steps < 0:
        raise ValueError("steps must be non-negative")

    if report_every < 1:
        raise ValueError("report_every must be positive")

    area = model.params.cell_area

    initial_debris = model.total_volume()

    initial_solid = float((model.thickness.sum() + model.soil_depth.sum()) * area)

    print(
        f"initial debris={initial_debris:.2f} m³; "
        f"debris+erodible soil={initial_solid:.2f} m³; "
        f"grid={model.height}x{model.width}; "
        f"dt={model.params.pt:g} s"
    )

    iterator = range(steps)
    if progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(range(steps), desc="steps", unit="step", leave=False)
        except ImportError:
            pass

    for i in iterator:
        model.step()

        if i % report_every == 0 or i == steps - 1:
            current_solid = float((model.thickness.sum() + model.soil_depth.sum()) * area)

            accounted_solid = current_solid + model.mass_lost

            print(
                f"step={i + 1:4d} "
                f"t={model.time:7.1f} s "
                f"debris={model.total_volume():10.2f} m³ "
                f"max_th={model.thickness.max():6.3f} m "
                f"max_kh={model.kinetic_head.max():7.3f} m "
                f"active={(model.thickness > 0.05).sum():5d} cells "
                f"lost={model.mass_lost:9.2f} m³ "
                f"solid-balance-error={accounted_solid - initial_solid:+.3e} m³"
            )

    final_solid = float((model.thickness.sum() + model.soil_depth.sum()) * area)

    print(
        f"finished: "
        f"initial debris={initial_debris:.2f} m³, "
        f"final debris={model.total_volume():.2f} m³, "
        f"remaining soil={model.soil_depth.sum() * area:.2f} m³, "
        f"water-entry loss={model.mass_lost:.2f} m³, "
        f"solid balance error={final_solid + model.mass_lost - initial_solid:+.3e} m³"
    )


def animate(
    model: SciddicaSS3,
    x,
    y,
    steps_per_frame: int = 2,
    frames: int = 250,
    interval: int = 40,
    save: str | None = None,
    steps: int | None = None,
    progress: bool = True,
    max_pixels: int = 800,
    fps: int = 20,
    dpi: int = 100,
) -> object:
    """Animate debris thickness over a static hillshade, drawn as an image."""
    import shutil

    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.colors import ListedColormap

    if steps_per_frame < 1 or frames < 1 or interval < 1:
        raise ValueError("steps_per_frame, frames and interval must be positive")

    if steps is not None and steps >= 0:
        frames = max((steps + steps_per_frame - 1) // steps_per_frame, 1)

    p = model.params

    # North is up: plotted northing = -y, as in plot_summary.
    extent = (
        float(x.min()) - p.pc, float(x.max()) + p.pc,
        -float(y.max()) - p.pc, -float(y.min()) + p.pc,
    )
    lookup = _hex_lookup(model.altitude.shape, p.pc, (0.0, 0.0), extent, max_pixels)
    valid = lookup[2]

    # ---- static layers, computed once ------------------------------------
    shade = _to_image(hillshade(model.altitude + model.thickness, p.pc), lookup)
    base = np.ma.masked_invalid(0.30 + 0.70 * shade)

    lake = None
    if model.water_level is not None:
        in_lake = valid & (_to_image(model.altitude, lookup) < model.water_level)
        if in_lake.any():
            lake = np.ma.masked_where(~in_lake, np.ones(valid.shape))

    fig, ax = plt.subplots(figsize=(8, 8))

    land_cmap = plt.get_cmap("Greys").copy()
    land_cmap.set_bad((1, 1, 1, 0))
    ax.imshow(base, extent=extent, origin="upper", cmap=land_cmap,
              vmin=0, vmax=1, interpolation="nearest", aspect="equal")

    if lake is not None:
        ax.imshow(lake, extent=extent, origin="upper",
                  cmap=ListedColormap([(0.55, 0.75, 0.95, 0.6)]),
                  interpolation="nearest", aspect="equal")

    # ---- dynamic layer ----------------------------------------------------
    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad((0, 0, 0, 0))

    def debris_image():
        img = _to_image(model.thickness, lookup, 0.0)
        return np.ma.masked_less(img, 0.02)

    debris = ax.imshow(
        debris_image(), extent=extent, origin="upper", cmap=cmap,
        vmin=0.0, vmax=max(4.0, float(model.thickness.max())),
        interpolation="nearest", aspect="equal",
    )
    fig.colorbar(debris, ax=ax, shrink=0.7, label="debris thickness [m]")

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])

    # Inside the axes, so blitting refreshes it.
    label = ax.text(
        0.02, 0.97, "t = 0 s", transform=ax.transAxes, va="top", fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="none"),
    )

    pbar = None
    if progress and save:
        try:
            from tqdm import tqdm

            pbar = tqdm(total=frames, desc="animating", unit="frame")
        except ImportError:
            pass

    def init():
        return debris, label

    def update(_frame: int):
        for _ in range(steps_per_frame):
            model.step()

        debris.set_data(debris_image())
        label.set_text(f"t = {model.time:.0f} s   debris volume = {model.total_volume():.0f} m³")

        if pbar is not None:
            pbar.update(1)

        return debris, label

    animation = FuncAnimation(
        fig, update, frames=frames, init_func=init, interval=interval,
        blit=save is None, cache_frame_data=False,
    )

    if save:
        streaming = str(save).lower().endswith((".mp4", ".mkv", ".webm"))
        writer = "ffmpeg" if streaming and shutil.which("ffmpeg") else "pillow"
        animation.save(save, writer=writer, fps=fps, dpi=dpi)
        plt.close(fig)
        if pbar is not None:
            pbar.close()
    else:
        plt.show()

    return animation


def hillshade(altitude, side_length, azimuth_deg=315.0, sun_deg=40.0, vert_exag=2.0):
    """Hillshade (0..1) of a terrain stored on the axial hex grid, north up.

    Gradients are taken in axial coordinates and converted to Cartesian ones
    using x = sqrt3*pc*(q + r/2), y = 1.5*pc*r.
    """
    dz_dr, dz_dq = np.gradient(np.asarray(altitude, float))
    dz_dx = dz_dq / (SQRT3 * side_length)
    dz_dy_down = (dz_dr - 0.5 * dz_dq) / (1.5 * side_length)
    dz_dn = -dz_dy_down                      # rows increase southwards
    nx, ny, nz = -vert_exag * dz_dx, -vert_exag * dz_dn, np.ones_like(dz_dx)
    norm = np.sqrt(nx**2 + ny**2 + nz**2)
    az, el = np.deg2rad(azimuth_deg), np.deg2rad(sun_deg)
    lx, ly, lz = np.cos(el) * np.sin(az), np.cos(el) * np.cos(az), np.sin(el)
    return np.clip((nx * lx + ny * ly + nz * lz) / norm, 0.0, 1.0)


def _nice_levels(vmax):
    """Binned colour levels, like the stepped colour bars of the reference."""
    candidates = np.array([0.05, 0.5, 1, 2, 3, 5, 7.5, 10, 15, 20, 30, 50, 75, 100])
    levels = list(candidates[candidates < vmax])
    above = candidates[candidates >= vmax]
    levels.append(float(above[0]) if above.size else float(vmax) * 1.1)
    if len(levels) < 3:
        levels = [0.05, max(vmax / 2, 0.1), max(vmax, 0.2)]
    return np.array(levels)


def _outline(ax, tri, mask, **style):
    """Draw the boundary of a boolean mask (skipped if empty or full)."""
    if mask.any() and not mask.all():
        ax.tricontour(tri, mask.astype(float).ravel(), levels=[0.5], **style)


def _hex_lookup(shape, pc, origin, extent, max_pixels=1100):
    """For every pixel of a regular north-up image, find the hex cell under it.

    Returns (ri, qi, valid, east, north); ri/qi index arrays shaped like the image.
    """
    H, W = shape
    e0, e1, n0, n1 = extent
    span_e, span_n = e1 - e0, n1 - n0
    if span_e >= span_n:
        ncols, nrows = max_pixels, max(int(max_pixels * span_n / span_e), 2)
    else:
        nrows, ncols = max_pixels, max(int(max_pixels * span_e / span_n), 2)

    east = np.linspace(e0, e1, ncols)
    north = np.linspace(n1, n0, nrows)               # first image row = northernmost
    X = east[None, :] - origin[0]
    Y = origin[1] - north[:, None]                   # model y grows with row index

    r_f = Y / (1.5 * pc)
    q_f = X / (SQRT3 * pc) - 0.5 * r_f

    # cube rounding of the fractional axial coordinates
    cx, cz = q_f, r_f
    cy = -cx - cz
    rx, ry, rz = np.rint(cx), np.rint(cy), np.rint(cz)
    dx, dy, dz = np.abs(rx - cx), np.abs(ry - cy), np.abs(rz - cz)
    fix_x = (dx > dy) & (dx > dz)
    fix_y = ~fix_x & (dy > dz)
    qi = np.where(fix_x, -ry - rz, rx).astype(np.int64)
    ri = np.where(~fix_x & ~fix_y, -rx - ry, rz).astype(np.int64)

    valid = (ri >= 0) & (ri < H) & (qi >= 0) & (qi < W)
    return (np.clip(ri, 0, H - 1), np.clip(qi, 0, W - 1), valid,
            np.broadcast_to(east, (nrows, ncols)), np.broadcast_to(north[:, None], (nrows, ncols)))


def _to_image(field, lookup, fill=np.nan):
    ri, qi, valid = lookup[0], lookup[1], lookup[2]
    return np.where(valid, np.asarray(field)[ri, qi], fill)


def plot_summary(
    model,
    x,
    y,
    max_thickness,
    initial_thickness,
    terrain,
    observed_impact=None,
    observed_deposit=None,
    origin=(1_489_000.0, 5_205_500.0),
    save: str | None = None,
    dpi: int = 200,
    max_pixels: int = 1100,
):
    """Two-panel summary (a) maximum flow height, (b) final deposit.

    Fields are resampled to a regular image of at most `max_pixels` on its long
    side, so cost and memory do not grow with the number of cells.
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter

    p = model.params
    extent = (
        origin[0] + float(x.min()) - p.pc, origin[0] + float(x.max()) + p.pc,
        origin[1] - float(y.max()) - p.pc, origin[1] - float(y.min()) + p.pc,
    )
    lookup = _hex_lookup(model.altitude.shape, p.pc, origin, extent, max_pixels)
    east_img, north_img = lookup[3], lookup[4]
    valid = lookup[2]

    shade = _to_image(hillshade(terrain, p.pc), lookup)
    grey = np.ma.masked_invalid(0.30 + 0.70 * shade)
    if model.water_level is None:
        lake = np.ma.masked_all(valid.shape)
    else:
        in_lake = valid & (_to_image(model.altitude, lookup) < model.water_level)
        lake = np.ma.masked_where(~in_lake, np.ones(valid.shape))

    final_thickness = model.thickness
    simulated = observed_impact is None
    impact = (max_thickness > 0.05) if observed_impact is None else observed_impact
    deposit = (final_thickness > 0.05) if observed_deposit is None else observed_deposit
    obs_label = "Simulated" if simulated else "Observation"

    release_img = _to_image(np.asarray(initial_thickness) > 0.05, lookup, 0.0).astype(float)
    impact_img = _to_image(np.asarray(impact), lookup, 0.0).astype(float)
    deposit_img = _to_image(np.asarray(deposit), lookup, 0.0).astype(float)

    fig = plt.figure(figsize=(5.4, 11.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.045], wspace=0.07, hspace=0.30)
    fmt = FuncFormatter(lambda v, _: f"{v:,.0f}".replace(",", " "))
    land_cmap = plt.get_cmap("Greys").copy()
    land_cmap.set_bad((1, 1, 1, 0))
    lake_cmap = ListedColormap([(0.55, 0.75, 0.95, 0.45)])

    def outline(ax, img, **style):
        if np.nanmin(img) < 0.5 < np.nanmax(img):
            ax.contour(east_img, north_img, img, levels=[0.5], **style)

    def panel(row, field, title, cmap_name, cbar_label):
        ax = fig.add_subplot(gs[row, 0])
        cax = fig.add_subplot(gs[row, 1])

        ax.imshow(grey, extent=extent, origin="upper", cmap=land_cmap,
                  vmin=0, vmax=1, interpolation="nearest", aspect="equal")
        if np.ma.count(lake):
            ax.imshow(lake, extent=extent, origin="upper", cmap=lake_cmap,
                      interpolation="nearest", aspect="equal")

        field_img = _to_image(field, lookup)
        vmax = float(np.nanmax(field_img)) if np.isfinite(field_img).any() else 0.0
        levels = _nice_levels(vmax)
        n = len(levels) - 1
        cmap = ListedColormap(plt.get_cmap(cmap_name)(np.linspace(0.25, 0.95, n)))
        cmap.set_bad((0, 0, 0, 0))
        norm = BoundaryNorm(levels, n)
        shown = np.ma.masked_less(np.minimum(np.nan_to_num(field_img), levels[-1] - 1e-9), levels[0])
        ax.imshow(shown, extent=extent, origin="upper", cmap=cmap, norm=norm,
                  alpha=0.85, interpolation="nearest", aspect="equal")

        contour_levels = [lv for lv in levels[1:] if lv < vmax]
        if contour_levels:
            ax.contour(east_img, north_img, np.nan_to_num(field_img),
                       levels=contour_levels, colors="k", linewidths=0.4, alpha=0.65)

        outline(ax, release_img, colors="red", linestyles=":", linewidths=1.3)
        outline(ax, impact_img, colors="red", linestyles="-", linewidths=1.2)
        outline(ax, deposit_img, colors="gold", linestyles="--", linewidths=1.4)

        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
        ax.xaxis.set_major_formatter(fmt)
        ax.yaxis.set_major_formatter(fmt)
        ax.tick_params(labelsize=7)
        plt.setp(ax.get_yticklabels(), rotation=90, va="center")
        ax.set_title(title, fontsize=13, fontweight="bold", pad=6)

        cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=levels)
        cb.ax.set_yticklabels([f"{lv:g}" for lv in levels])
        cb.ax.tick_params(labelsize=7)
        cb.set_label(cbar_label, fontsize=8)

        handles = [
            Line2D([0], [0], color="red", ls=":", lw=1.3, label="Release area"),
            Line2D([0], [0], color="red", ls="-", lw=1.2, label=f"{obs_label}: impact area"),
            Line2D([0], [0], color="gold", ls="--", lw=1.4, label=f"{obs_label}: deposit"),
        ]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07),
                  ncol=3, frameon=False, fontsize=6.5, handlelength=1.8, columnspacing=1.0)

    panel(0, max_thickness, "(a) Maximum flow height", "YlOrBr", "Flow height [m]")
    panel(1, final_thickness, "(b) Height of final deposit", "Oranges", "Deposit thickness [m]")

    if save:
        fig.savefig(save, dpi=dpi)        # no bbox_inches="tight": it renders twice
    return fig
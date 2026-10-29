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
) -> object:
    """Animate terrain, dynamic water mask and debris thickness."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.collections import PolyCollection

    from .geometry import hex_vertices

    if steps_per_frame < 1 or frames < 1 or interval < 1:
        raise ValueError("steps_per_frame, frames and interval must be positive")

    if steps is not None and steps >= 0:
        frames = (steps + steps_per_frame - 1) // steps_per_frame

    p = model.params
    vertices = hex_vertices(x, y, p.pc)
    polygons = list(vertices)

    fig, ax = plt.subplots(figsize=(8, 8))

    land = PolyCollection(
        polygons,
        array=model.altitude.ravel(),
        cmap="terrain",
        edgecolors="none",
    )
    ax.add_collection(land)

    lake = PolyCollection(polygons, edgecolors="none")
    ax.add_collection(lake)

    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad((0, 0, 0, 0))

    debris = PolyCollection(polygons, cmap=cmap, edgecolors="none")

    debris.set_clim(0.0, max(4.0, float(model.thickness.max())))
    debris.set_array(np.ma.masked_less(model.thickness.ravel(), 0.02))

    ax.add_collection(debris)

    fig.colorbar(debris, ax=ax, shrink=0.7, label="debris thickness [m]")

    ax.set_xlim(x.min() - p.pc, x.max() + p.pc)
    ax.set_ylim(y.min() - p.pc, y.max() + p.pc)

    ax.set_aspect("equal")
    ax.invert_yaxis()

    title = ax.set_title("t = 0 s")

    pbar = None
    if progress and save:
        try:
            from tqdm import tqdm

            pbar = tqdm(total=frames, desc="animating", unit="frame")
        except ImportError:
            pass

    def update(_frame: int):
        for _ in range(steps_per_frame):
            model.step()

        land.set_array(model.altitude.ravel())

        water_alpha = model.submerged().ravel().astype(float) * 0.75
        water_colours = np.zeros((water_alpha.size, 4), dtype=float)
        water_colours[:, 0] = 0.55
        water_colours[:, 1] = 0.75
        water_colours[:, 2] = 0.95
        water_colours[:, 3] = water_alpha

        lake.set_facecolor([tuple(color) for color in water_colours])

        debris.set_array(np.ma.masked_less(model.thickness.ravel(), 0.02))

        title.set_text(f"t = {model.time:.0f} s   debris volume = {model.total_volume():.0f} m³")

        if pbar is not None:
            pbar.update(1)

        return land, lake, debris, title

    animation = FuncAnimation(fig, update, frames=frames, interval=interval, blit=False)

    if save:
        animation.save(save, writer="pillow", fps=20)
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


def plot_summary(
    model: SciddicaSS3,
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
):
    """Two-panel figure: (a) maximum flow height, (b) height of final deposit.

    max_thickness      per-cell maximum debris thickness reached during the run
    initial_thickness  thickness at t=0 (release area)
    terrain            pre-event surface (altitude + initial thickness), for hillshade
    observed_*         optional boolean masks; if None, simulated extents
                       (max_thickness > 0.05 and final deposit > 0.05) are drawn instead
    origin             (easting, northing) offset so axes look like UTM coordinates
    """
    import matplotlib.pyplot as plt
    from matplotlib.cm import ScalarMappable
    from matplotlib.collections import PolyCollection
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter
    from matplotlib.tri import Triangulation

    from .geometry import hex_vertices

    p = model.params
    east = origin[0] + x
    north = origin[1] - y                      # north up: top rows = highest northing
    verts = hex_vertices(east, north, p.pc)
    # Delaunay triangulation is translation-invariant, but QHull chokes on
    # the large UTM coordinates; triangulate centered coordinates then restore.
    _cx = east.ravel() - east.ravel().mean()
    _cy = north.ravel() - north.ravel().mean()
    tri = Triangulation(_cx, _cy)
    tri.x = east.ravel().copy()
    tri.y = north.ravel().copy()

    shade = hillshade(terrain, p.pc)
    grey = plt.cm.gray(0.30 + 0.70 * shade).reshape(-1, 4)
    lake_mask = model.submerged().ravel()

    final_thickness = model.thickness.copy()
    simulated = observed_impact is None
    impact = (max_thickness > 0.05) if observed_impact is None else observed_impact
    deposit = (final_thickness > 0.05) if observed_deposit is None else observed_deposit
    obs_label = "Simulated" if simulated else "Observation"

    fig = plt.figure(figsize=(5.4, 11.0))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 0.045], wspace=0.07, hspace=0.30)

    fmt = FuncFormatter(lambda v, _: f"{v:,.0f}".replace(",", " "))

    def panel(row, field, title, cmap_name, cbar_label):
        ax = fig.add_subplot(gs[row, 0])
        cax = fig.add_subplot(gs[row, 1])

        # hillshaded terrain (+ faint lake)
        ax.add_collection(PolyCollection(list(verts), facecolors=grey,
                                         edgecolors="face", linewidths=0.3))
        if lake_mask.any():
            ax.add_collection(PolyCollection(
            list(verts[lake_mask]), facecolors=(0.55, 0.75, 0.95, 0.45),
                edgecolors="none"))

        # binned colour layer
        levels = _nice_levels(float(field.max()))
        n = len(levels) - 1
        cmap = ListedColormap(plt.get_cmap(cmap_name)(np.linspace(0.25, 0.95, n)))
        cmap.set_bad((0, 0, 0, 0))
        norm = BoundaryNorm(levels, n)
        values = np.ma.masked_less(np.minimum(field, levels[-1] - 1e-9), levels[0])
        ax.add_collection(PolyCollection(
            list(verts), array=values.ravel(), cmap=cmap, norm=norm,
            edgecolors="face", linewidths=0.3, alpha=0.85))

        # thickness contours
        contour_levels = [lv for lv in levels[1:] if lv < field.max()]
        if contour_levels:
            ax.tricontour(tri, field.ravel(), levels=contour_levels,
                          colors="k", linewidths=0.4, alpha=0.65)

        # release area (dotted red), impact area (solid red), deposit (dashed yellow)
        _outline(ax, tri, initial_thickness > 0.05,
                 colors="red", linestyles=":", linewidths=1.3)
        _outline(ax, tri, impact, colors="red", linestyles="-", linewidths=1.2)
        _outline(ax, tri, deposit, colors="gold", linestyles="--", linewidths=1.4)

        ax.set_xlim(east.min() - p.pc, east.max() + p.pc)
        ax.set_ylim(north.min() - p.pc, north.max() + p.pc)
        ax.set_aspect("equal")
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
                  ncol=3, frameon=False, fontsize=6.5, handlelength=1.8,
                  columnspacing=1.0)
        return ax

    panel(0, max_thickness, "(a) Maximum flow height", "YlOrBr",
          "Flow height [m]")
    panel(1, final_thickness, "(b) Height of final deposit", "Oranges",
          "Deposit thickness [m]")

    if save:
        fig.savefig(save, dpi=dpi, bbox_inches="tight")
    return fig
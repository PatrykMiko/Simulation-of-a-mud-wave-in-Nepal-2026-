#!/usr/bin/env python3
"""Headless diagnostics and visualization."""

from __future__ import annotations

import numpy as np

from .automaton import SciddicaSS3


def run_headless(
    model: SciddicaSS3,
    steps: int = 300,
    report_every: int = 50,
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

    for i in range(steps):
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
) -> object:
    """Animate terrain, dynamic water mask and debris thickness."""
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.collections import PolyCollection

    from .geometry import hex_vertices

    if steps_per_frame < 1 or frames < 1 or interval < 1:
        raise ValueError("steps_per_frame, frames and interval must be positive")

    p = model.params
    vertices = hex_vertices(x, y, p.pc)

    fig, ax = plt.subplots(figsize=(8, 8))

    land = PolyCollection(
        vertices,
        array=model.altitude.ravel(),
        cmap="Greys",
        edgecolors="none",
    )
    ax.add_collection(land)

    lake = PolyCollection(vertices, edgecolors="none")
    ax.add_collection(lake)

    cmap = plt.get_cmap("inferno").copy()
    cmap.set_bad((0, 0, 0, 0))

    debris = PolyCollection(vertices, cmap=cmap, edgecolors="none")

    debris.set_clim(0.0, max(4.0, float(model.thickness.max())))
    debris.set_array(np.ma.masked_less(model.thickness.ravel(), 0.02))

    ax.add_collection(debris)

    fig.colorbar(debris, ax=ax, shrink=0.7, label="debris thickness [m]")

    ax.set_xlim(x.min() - p.pc, x.max() + p.pc)
    ax.set_ylim(y.min() - p.pc, y.max() + p.pc)

    ax.set_aspect("equal")
    ax.invert_yaxis()

    title = ax.set_title("t = 0 s")

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

        lake.set_facecolors(water_colours)

        debris.set_array(np.ma.masked_less(model.thickness.ravel(), 0.02))

        title.set_text(f"t = {model.time:.0f} s   debris volume = {model.total_volume():.0f} m³")

        return land, lake, debris, title

    animation = FuncAnimation(fig, update, frames=frames, interval=interval, blit=False)

    if save:
        animation.save(save, writer="pillow", fps=20)
        plt.close(fig)
    else:
        plt.show()

    return animation

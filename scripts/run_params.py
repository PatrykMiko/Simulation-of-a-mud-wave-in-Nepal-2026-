#!/usr/bin/env python3
"""Run SCIDDICA-SS3 with a fixed parameter set and save an animation.

Usage:
    python scripts/run_params.py [SCENE_DIR] --save out.gif
    python scripts/run_params.py [SCENE_DIR] --live        # interactive window

The parameters below are the best Optuna trial (#380, IoU 0.200) for the
Nepal ALOS scene.  They are passed explicitly to ModelParameters so the CLI
defaults (adh=0.02, pml=0.25, pwr=0.10, alpha_scale=0.25, ...) are overridden.

The simulation runs to stability (same early-stopping rule as optimize.run_model):
mass change and extent change both drop below 1e-4 for 5 consecutive checks,
or max_steps is reached.  The animation then replays the same number of steps.
"""

from __future__ import annotations

import argparse
from pathlib import Path

from sciddica_ss3 import SciddicaSS3, ModelParameters, cell_centres, animate
from sciddica_ss3.constants import MIN_THICKNESS
from sciddica_ss3.scenario import jaccard_fitness
from sciddica_ss3.state import MediumParameters

# Best-trial parameters (trial #380, IoU 0.200387).
BEST_PARAMS = dict(
    adh=0.0243959613216525,
    alpha_scale=1.46399702349686,
    ped=0.344222113524263,
    per=0.195483644405376,
    pf=0.198663995707226,
    pml=0.0745233160899741,
    pmt=1.15032241934698,
    ptd=0.259957200943208,
    pwr=0.506162705978809,
    shift_friction=0.0219516603518123,
)


def build_params(pc: float) -> ModelParameters:
    return ModelParameters(
        pc=pc,
        air=MediumParameters(
            adh=BEST_PARAMS["adh"],
            pf=BEST_PARAMS["pf"],
            pmt=BEST_PARAMS["pmt"],
            per=BEST_PARAMS["per"],
            ptd=BEST_PARAMS["ptd"],
            ped=BEST_PARAMS["ped"],
            shift_friction=BEST_PARAMS["shift_friction"],
        ),
        pml=BEST_PARAMS["pml"],
        pwr=BEST_PARAMS["pwr"],
        alpha_scale=BEST_PARAMS["alpha_scale"],
    )


def run_to_stability(model, max_steps=10000, check_every=100,
                     stable_tolerance=1e-4, stable_window=5, log_every=100):
    """Mirror optimize.run_model's early-stopping rule on a live model."""
    previous_mass = float(model.total_volume())
    previous_mask = model.thickness > MIN_THICKNESS
    stable_count = 0
    step = 0
    while step < max_steps:
        model.step()
        step += 1
        if check_every and step % check_every == 0:
            current_mass = float(model.total_volume())
            current_mask = model.thickness > MIN_THICKNESS
            mass_change = abs(current_mass - previous_mass) / max(previous_mass, 1e-12)
            extent_change = 1.0 - jaccard_fitness(previous_mask, current_mask)
            if log_every and step % log_every == 0:
                print(
                    "step=%d mass=%.6f mass_change=%.3e extent_change=%.3e "
                    "active_cells=%d stable=%d/%d" % (
                        step, current_mass, mass_change, extent_change,
                        int(current_mask.sum()), stable_count, stable_window),
                    flush=True,
                )
            if mass_change < stable_tolerance and extent_change < stable_tolerance:
                stable_count += 1
                if stable_count >= stable_window:
                    break
            else:
                stable_count = 0
            previous_mass = current_mass
            previous_mask = current_mask
    return step


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("source", nargs="?", default="data/nepal_ss3/alos",
                    help="path to SS3 .npz or directory with ss3_arrays.npz")
    ap.add_argument("--save", metavar="FILE.gif", help="save animation")
    ap.add_argument("--live", action="store_true",
                    help="show the animation live in a window (default)")
    ap.add_argument("--max-steps", type=int, default=10000,
                    help="hard cap on simulation steps (default: 10000)")
    ap.add_argument("--use-float32", action="store_true")
    ap.add_argument("--use-numba", action="store_true")
    ap.add_argument("--use-cpp", action="store_true")
    ap.add_argument("--release-scale", type=float, default=None)
    ap.add_argument("--release-volume", type=float, default=None)
    ap.add_argument("--frames", type=int, default=250,
                    help="max animation frames (default: 250)")
    ap.add_argument("--steps-per-frame", type=int, default=2,
                    help="steps per animation frame (default: 2)")
    a = ap.parse_args()

    if a.release_scale is not None and a.release_volume is not None:
        raise SystemExit("use --release-scale or --release-volume, not both")

    import json
    import numpy as np

    p = Path(a.source)
    if p.is_dir():
        p = p / "ss3_arrays.npz"
    data = np.load(p)
    meta_path = p.parent / "ss3_meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    side = float(meta.get("cellsize", 12.0))

    params = build_params(side)
    model = SciddicaSS3(
        data["altitude"], data["thickness"], data["soil_depth"], params,
        water_level=meta.get("water_level"),
        use_float32=a.use_float32,
        use_numba=a.use_numba,
        use_cpp=a.use_cpp,
    )
    x, y = cell_centres(model.height, model.width, side)

    if a.release_scale is not None:
        model.scale_release(factor=a.release_scale)
    elif a.release_volume is not None:
        model.scale_release(volume=a.release_volume)

    print("params:", BEST_PARAMS)
    print("pc=%.3f  grid=%dx%d  volume=%.2f m3" % (
        side, model.height, model.width, model.total_volume()))

    p = model.params
    extent = (
        float(x.min()) - p.pc, float(x.max()) + p.pc,
        -float(y.max()) - p.pc, -float(y.min()) + p.pc,
    )

    if a.save:
        import matplotlib
        matplotlib.use("Agg")

    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.colors import ListedColormap
    from sciddica_ss3.visualize import hillshade, _hex_lookup, _to_image

    lookup = _hex_lookup(model.altitude.shape, p.pc, (0.0, 0.0), extent, 1100)
    valid = lookup[2]
    shade = _to_image(hillshade(model.altitude + model.thickness, p.pc), lookup)
    base = np.ma.masked_invalid(0.30 + 0.70 * shade)

    fig, ax = plt.subplots(figsize=(8, 8))
    land_cmap = plt.get_cmap("Greys").copy()
    land_cmap.set_bad((1, 1, 1, 0))
    ax.imshow(base, extent=extent, origin="upper", cmap=land_cmap,
              vmin=0, vmax=1, interpolation="nearest", aspect="equal")

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

    label = ax.text(
        0.02, 0.97, "t = 0 s", transform=ax.transAxes, va="top", fontsize=10,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.75, edgecolor="none"),
    )

    previous_mass = float(model.total_volume())
    previous_mask = model.thickness > MIN_THICKNESS
    stable_count = 0
    step = 0
    max_steps = a.max_steps
    check_every = 100
    stable_tolerance = 1e-4
    stable_window = 5

    def init():
        return debris, label

    def update(_frame):
        nonlocal step, previous_mass, previous_mask, stable_count
        for _ in range(a.steps_per_frame):
            model.step()
            step += 1
            if step >= max_steps:
                break
            if check_every and step % check_every == 0:
                current_mass = float(model.total_volume())
                current_mask = model.thickness > MIN_THICKNESS
                mass_change = abs(current_mass - previous_mass) / max(previous_mass, 1e-12)
                extent_change = 1.0 - jaccard_fitness(previous_mask, current_mask)
                if mass_change < stable_tolerance and extent_change < stable_tolerance:
                    stable_count += 1
                    if stable_count >= stable_window:
                        return debris, label
                else:
                    stable_count = 0
                previous_mass = current_mass
                previous_mask = current_mask
        debris.set_data(debris_image())
        label.set_text(
            "t = %.0f s   step=%d   volume = %.0f m3" % (
                model.time, step, model.total_volume()))
        return debris, label

    anim = FuncAnimation(
        fig, update, frames=a.frames, init_func=init, interval=40,
        blit=a.save is None, cache_frame_data=False,
    )

    if a.save:
        import shutil
        streaming = str(a.save).lower().endswith((".mp4", ".mkv", ".webm"))
        writer = "ffmpeg" if streaming and shutil.which("ffmpeg") else "pillow"
        anim.save(a.save, writer=writer, fps=20, dpi=100)
        plt.close(fig)
        print("saved animation to %s (%d steps)" % (a.save, step))
    else:
        print("live simulation: %d steps, t=%.0f s, volume=%.0f m3" % (
            step, model.time, model.total_volume()))
        plt.show()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

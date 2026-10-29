"""Run the SCIDDICA-SS3 simulation on a synthetic or converted scenario."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import numpy as np

from . import make_albano_like, run_headless, animate, plot_summary
from .automaton import SciddicaSS3
from .geometry import cell_centres
from .state import ModelParameters


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument(
        "source",
        nargs="?",
        help="path to converted SS3 data (.npz or directory with ss3_arrays.npz)",
    )
    parser.add_argument(
        "--headless",
        action="store_true",
        help="run without opening a window",
    )
    parser.add_argument(
        "--save",
        metavar="FILE.gif",
        help="save animation as a GIF",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=80,
        help="grid width and height (default: 80)",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=300,
        help="number of steps for --headless, --summary and --mask-out (default: 300)",
    )
    parser.add_argument(
        "--report-every",
        type=int,
        default=50,
        help="headless report interval",
    )
    parser.add_argument(
        "--mask-out",
        metavar="FILE",
        help="save final debris mask (thickness>0.1) to .npy or .tif",
    )
    parser.add_argument(
        "--summary",
        metavar="FILE.png",
        help="run --steps steps and save the two-panel summary figure",
    )
    parser.add_argument(
        "--alpha-scale",
        type=float,
        default=0.3,
        help="momentum elevation correction scale (0 disables it)",
    )
    parser.add_argument(
        "--release-scale",
        type=float,
        metavar="F",
        help="multiply the initial debris thickness by F",
    )
    parser.add_argument(
        "--release-volume",
        type=float,
        metavar="M3",
        help="rescale the initial debris to this volume in m³ (hexagon cell area)",
    )

    return parser.parse_args(argv)


def save_mask(mask: np.ndarray, path: str) -> None:
    """Save a boolean mask as .npy (default) or .tif/.tiff (via Pillow)."""
    target = Path(path)

    if target.suffix.lower() in {".tif", ".tiff"}:
        from PIL import Image

        Image.fromarray(mask.astype(np.uint8) * 255).save(target)
    else:
        if target.suffix.lower() != ".npy":
            target = target.with_suffix(".npy")
        np.save(target, mask)

    print(f"saved mask to {target}")


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)

    if args.release_scale is not None and args.release_volume is not None:
        raise SystemExit("use --release-scale or --release-volume, not both")

    if args.headless or args.save or args.summary or args.mask_out:
        import matplotlib

        matplotlib.use("Agg")

    if args.source:
        p = Path(args.source)
        if p.is_dir():
            p = p / "ss3_arrays.npz"
        data = np.load(p)
        meta_path = p.parent / "ss3_meta.json"
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        side = float(meta.get("cellsize", 12.0))
        params = ModelParameters(
            pc=side,
            alpha_scale=args.alpha_scale,
        )
        model = SciddicaSS3(
            data["altitude"],
            data["thickness"],
            data["soil_depth"],
            params,
            water_level=meta.get("water_level"),
        )
        x, y = cell_centres(model.height, model.width, side)
    else:
        params = ModelParameters(alpha_scale=args.alpha_scale)
        model, x, y = make_albano_like(args.n, params=params)

    keep_surface = not args.source

    if args.release_scale is not None:
        factor = model.scale_release(factor=args.release_scale, keep_surface=keep_surface)
        print(f"release scaled by {factor:.4g}; volume = {model.total_volume():.4g} m³")
    elif args.release_volume is not None:
        factor = model.scale_release(volume=args.release_volume, keep_surface=keep_surface)
        print(f"release scaled by {factor:.4g}; volume = {model.total_volume():.4g} m³")

    if args.summary:
        initial = model.thickness.copy()
        terrain = model.altitude + initial          # pre-event surface
        max_thickness = initial.copy()

        def track(m: SciddicaSS3, _i: int) -> None:
            np.maximum(max_thickness, m.thickness, out=max_thickness)

        model.run(args.steps, callback=track)
        plot_summary(model, x, y, max_thickness, initial, terrain, save=args.summary)
        print(f"saved summary figure to {args.summary}")

    elif args.headless:
        run_headless(model, steps=args.steps, report_every=args.report_every)

    elif args.save:
        animate(model, x, y, save=args.save, steps=args.steps)

    elif args.mask_out:
        # Only a mask was requested: just run the model, no window.
        model.run(args.steps)

    else:
        animate(model, x, y)

    if args.mask_out:
        save_mask(model.thickness > 0.1, args.mask_out)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
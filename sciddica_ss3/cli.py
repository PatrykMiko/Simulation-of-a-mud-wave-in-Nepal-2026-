"""Run the SCIDDICA-SS3 simulation on a synthetic or converted scenario."""

from __future__ import annotations

import argparse
from typing import Iterable

import json
from pathlib import Path

import numpy as np

from sciddica_ss3.geometry import cell_centres
from sciddica_ss3.state import ModelParameters

from . import make_albano_like, run_headless, animate
from .automaton import SciddicaSS3


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
        help="headless steps (default: 300)",
    )
    parser.add_argument(
        "--report-every",
        type=int,
        default=50,
        help="headless report interval",
    )

    return parser.parse_args(argv)


def main(argv: Iterable[str] | None = None) -> int:
    args = parse_args(argv)

    if args.headless or args.save:
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
        params = ModelParameters(pc=side)
        model = SciddicaSS3(
            data["altitude"],
            data["thickness"],
            data["soil_depth"],
            params,
            water_level=meta.get("water_level"),
        )
        x, y = cell_centres(model.height, model.width, side)
    else:
        model, x, y = make_albano_like(args.n)
        
    if args.headless:
        run_headless(model, steps=args.steps, report_every=args.report_every)
    else:
        animate(model, x, y, save=args.save)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
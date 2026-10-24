from __future__ import annotations

import argparse
from typing import Iterable

from . import make_albano_like, run_headless, animate


def parse_args(argv: Iterable[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)

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

    model, x, y = make_albano_like(args.n)

    if args.headless:
        run_headless(model, steps=args.steps, report_every=args.report_every)
    else:
        animate(model, x, y, save=args.save)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
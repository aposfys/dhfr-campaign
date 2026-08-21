"""Command line entry point: ``python -m dhfrcamp.cli`` or ``dhfrcamp``."""

from __future__ import annotations

import argparse
from pathlib import Path

from dhfrcamp import __version__


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dhfrcamp",
        description="Retrospective virtual screening campaign against human DHFR",
    )
    parser.add_argument("--version", action="version", version=f"dhfrcamp {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("prepare", help="fetch structures and known actives, strip additives")

    decoys = sub.add_parser("decoys", help="generate property-matched decoys")
    decoys.add_argument("--per-active", type=int, default=50)

    screen = sub.add_parser("screen", help="co-folding and affinity scoring (resumable)")
    screen.add_argument("--batch-size", type=int, default=8)
    screen.add_argument(
        "--resume",
        action="store_true",
        default=True,
        help="skip ligands already scored; on by default because this step is expensive",
    )

    generate = sub.add_parser("generate", help="REINVENT 4 expansion around confirmed hits")
    generate.add_argument("--n-designs", type=int, default=1000)

    evaluate = sub.add_parser("evaluate", help="enrichment against the matched-decoy null")
    evaluate.add_argument("--fraction", type=float, default=0.01)
    evaluate.add_argument("--alpha", type=float, default=20.0)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    raise SystemExit(f"'{args.command}' is not implemented yet; see README milestones")


if __name__ == "__main__":
    raise SystemExit(main())

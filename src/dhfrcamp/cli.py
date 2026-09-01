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

    campaign = sub.add_parser(
        "campaign", help="the decoy-bias experiment: one screen, two decoy sets"
    )
    campaign.add_argument("--catalog", type=Path, default=Path("data/catalog.sqlite"))
    campaign.add_argument("--pool-size", type=int, default=60_000)
    campaign.add_argument("--per-active", type=int, default=50)
    campaign.add_argument("--max-actives", type=int, default=150)
    campaign.add_argument("--seed", type=int, default=0)

    evaluate = sub.add_parser("evaluate", help="enrichment against the matched-decoy null")
    evaluate.add_argument("--fraction", type=float, default=0.01)
    evaluate.add_argument("--alpha", type=float, default=20.0)

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "prepare":
        from dhfrcamp.prepare import COMPLEXES, binding_site

        for pdb_id in COMPLEXES:
            site = sorted(binding_site(pdb_id, data_dir=args.data_dir / "pdb"))
            print(f"{pdb_id} ({COMPLEXES[pdb_id]}): {len(site)} residues -> {site}")
        return 0

    if args.command == "campaign":
        from dhfrcamp.campaign import run

        run(
            args.catalog,
            args.results_dir,
            pool_size=args.pool_size,
            decoys_per_active=args.per_active,
            max_actives=args.max_actives,
            seed=args.seed,
        )
        print(f"wrote {args.results_dir / 'findings.json'}")
        return 0

    raise SystemExit(
        f"'{args.command}' is not implemented: the co-folding screen needs a GPU. "
        "Run 'campaign' for the decoy-bias experiment that does not."
    )


if __name__ == "__main__":
    raise SystemExit(main())

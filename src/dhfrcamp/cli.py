"""Command line entry point: ``python -m dhfrcamp.cli`` or ``dhfrcamp``.

Two rules apply to every subcommand here.

A missing prerequisite is reported, not raised. Running ``campaign`` without a catalogue
used to surface ``sqlite3.OperationalError: unable to open database file``, which tells the
user nothing about what to do next; preconditions are now checked up front and reported with
the command that fixes them.

A subcommand is refused only for the reason it is actually blocked. ``screen`` and
``generate`` need a GPU. ``decoys`` and ``evaluate`` do not, and used to be refused with a
GPU message anyway even though both are implemented.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from dhfrcamp import __version__, statistics

#: Subcommands that genuinely require hardware this project has never had.
GPU_GATED = {
    "screen": "Boltz-2 co-folding with an affinity head",
    "generate": "REINVENT 4 generation",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dhfrcamp",
        description="Retrospective virtual screening campaign against human DHFR",
    )
    parser.add_argument("--version", action="version", version=f"dhfrcamp {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("prepare", help="fetch structures, strip additives, define the site")

    decoys = sub.add_parser(
        "decoys", help="generate property-matched decoys and report the match"
    )
    decoys.add_argument("--catalog", type=Path, default=Path("data/catalog.sqlite"))
    decoys.add_argument("--per-active", type=int, default=50)
    decoys.add_argument("--pool-size", type=int, default=60_000)
    decoys.add_argument("--max-actives", type=int, default=150)
    decoys.add_argument("--seed", type=int, default=0)

    screen = sub.add_parser("screen", help="co-folding and affinity scoring (needs a GPU)")
    screen.add_argument("--batch-size", type=int, default=8)

    generate = sub.add_parser("generate", help="REINVENT 4 expansion (needs a GPU)")
    generate.add_argument("--n-designs", type=int, default=1000)

    campaign = sub.add_parser(
        "campaign", help="the decoy-bias experiment: one screen, two decoy sets"
    )
    campaign.add_argument("--catalog", type=Path, default=Path("data/catalog.sqlite"))
    campaign.add_argument("--pool-size", type=int, default=60_000)
    campaign.add_argument("--per-active", type=int, default=50)
    campaign.add_argument("--max-actives", type=int, default=150)
    campaign.add_argument("--seed", type=int, default=0)

    sub.add_parser("evaluate", help="print the enrichment table from an existing run")

    return parser


def _require_catalog(path: Path) -> None:
    """Refuse early, with the command that produces the missing file."""
    if not path.exists():
        raise SystemExit(
            f"catalogue not found at {path}.\n"
            "Build it from a ChEMBL chemreps dump:\n"
            "  curl -O https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/"
            "chembl_36/chembl_36_chemreps.txt.gz\n"
            "  python3 tools/build_catalog.py chembl_36_chemreps.txt.gz "
            f"{path}\n"
            "(tools/build_catalog.py ships with this repository.)"
        )


def _run_campaign(args: argparse.Namespace) -> dict:
    from dhfrcamp.campaign import run

    _require_catalog(args.catalog)
    try:
        return run(
            args.catalog,
            args.results_dir,
            pool_size=args.pool_size,
            decoys_per_active=args.per_active,
            max_actives=args.max_actives,
            seed=args.seed,
        )
    except OSError as exc:
        # Network or filesystem: real, expected, and not worth a traceback.
        raise SystemExit(f"campaign could not run: {exc}") from exc


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command in GPU_GATED:
        raise SystemExit(
            f"'{args.command}' is not implemented: {GPU_GATED[args.command]} needs a GPU "
            "this project has not had access to.\n"
            "Run 'dhfrcamp campaign' for the decoy-bias experiment, which does not."
        )

    if args.command == "prepare":
        from dhfrcamp.prepare import COMPLEXES, binding_site

        try:
            for pdb_id in COMPLEXES:
                site = sorted(binding_site(pdb_id, data_dir=args.data_dir / "pdb"))
                print(f"{pdb_id} ({COMPLEXES[pdb_id]}): {len(site)} residues -> {site}")
        except OSError as exc:
            raise SystemExit(f"could not fetch structures from RCSB: {exc}") from exc
        return 0

    if args.command == "decoys":
        # Implemented and CPU-only. This used to be refused with a GPU message.
        findings = _run_campaign(args)
        report = findings["match_report"]
        print(
            f"\n{report['n_decoys']:,} property-matched decoys for {report['n_actives']} "
            f"actives ({len(report['unmatched_actives'])} actives short of a full set)"
        )
        for prop, worst in sorted(report["max_abs_difference"].items()):
            print(f"  {prop:<18} max difference {worst:.3f}")
        return 0

    if args.command == "campaign":
        _run_campaign(args)
        print(f"wrote {args.results_dir / 'findings.json'}")
        return 0

    if args.command == "evaluate":
        # Also implemented and CPU-only: it reads what campaign already produced.
        path = args.results_dir / "findings.json"
        if not path.exists():
            raise SystemExit(f"no findings at {path}. Run 'dhfrcamp campaign' first.")
        findings = json.loads(path.read_text())
        print(
            f"target {findings['configuration']['target']}, "
            f"screen: {findings['configuration']['screen']}"
        )
        print(f"{'arm':<18} {'EF 1%':>8} {'ceiling':>8} {'BEDROC':>8} {'prop AUC':>9}")
        for arm in findings["arms"]:
            print(
                f"{arm['name']:<18} {arm['ef_1pct']:>8.2f} {arm['ef_1pct_ceiling']:>8.2f} "
                f"{arm['bedroc_20']:>8.3f} {arm['property_auc']:>9.3f}"
            )
        report = statistics.bias_report(findings["arms"])
        print("\nproperty-only AUC with uncertainty (Hanley-McNeil 95%):")
        for name, entry in report["per_arm"].items():
            chance = entry["vs_chance"]
            print(
                f"  {name:<18} {entry['auc']:.4f} "
                f"[{entry['ci_lower']:.4f}, {entry['ci_upper']:.4f}]  "
                f"vs chance z={chance['z']:.1f} p={chance['p_value']:.2e}"
            )
        effect = report.get("matching_effect")
        if effect:
            print(
                f"\n  matching removes {effect['difference']:.4f} +- "
                f"{effect['standard_error']:.4f} of property-only AUC "
                f"(z={effect['z']:.2f}, p={effect['p_value']:.2e}), and the "
                f"matched arm still sits "
                f"{report['per_arm']['property_matched']['vs_chance']['excess_over_chance']:.4f} "
                "above chance."
            )
        print(f"\nnot run: {findings['configuration']['screen_not_run']}")
        return 0

    raise SystemExit(f"unknown command {args.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())

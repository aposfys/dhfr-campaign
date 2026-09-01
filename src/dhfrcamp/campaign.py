"""The campaign: the same screen scored against two different denominators.

One method, one set of actives, two decoy sets -- property-matched and unmatched. Nothing
varies except the denominator, so any difference in enrichment is attributable to the
decoys and to nothing else. That is the entire experiment.
"""

from __future__ import annotations

import json
import random
import sqlite3
import time
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path

from dhfrcamp.decoys import (
    MATCHED_PROPERTIES,
    compute_properties,
    generate_pool,
    match_decoys,
    sample_unmatched,
)
from dhfrcamp.evaluate import (
    bedroc,
    enrichment_factor,
    max_enrichment_factor,
    property_separability,
)
from dhfrcamp.screen import fingerprints, similarity_screen

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"

#: Human dihydrofolate reductase.
DHFR_TARGET = "CHEMBL202"

#: pChEMBL at or above which a compound counts as an active. 6.0 is 1 uM, the usual line.
ACTIVE_THRESHOLD = 6.0


@dataclass
class Arm:
    """One decoy set and what the screen scored against it."""

    name: str
    n_actives: int
    n_decoys: int
    ef_1pct: float
    ef_1pct_ceiling: float
    ef_5pct: float
    bedroc_20: float
    property_auc: float


def fetch_actives(out_path: Path, *, target: str = DHFR_TARGET) -> list[dict]:
    """Every ChEMBL compound with a recorded potency against human DHFR."""
    if out_path.exists():
        return json.loads(out_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    query = {
        "target_chembl_id": target,
        "pchembl_value__isnull": "false",
        "limit": "1000",
    }
    url = f"{CHEMBL_API}/activity.json?{urllib.parse.urlencode(query)}"
    while url:
        with urllib.request.urlopen(url, timeout=180) as response:
            payload = json.loads(response.read().decode("utf-8"))
        rows.extend(payload["activities"])
        nxt = payload["page_meta"].get("next")
        url = f"https://www.ebi.ac.uk{nxt}" if nxt else ""

    out_path.write_text(json.dumps(rows, indent=1))
    return rows


def load_catalogue(catalog_db: Path, *, limit: int = 0) -> list[tuple[str, str]]:
    """`(id, smiles)` pairs from the ChEMBL structure catalogue."""
    connection = sqlite3.connect(f"file:{catalog_db}?mode=ro", uri=True)
    try:
        sql = "SELECT id, smiles FROM molecule"
        if limit:
            sql += f" LIMIT {int(limit)}"
        return [(f"CHEMBL{row[0]}", row[1]) for row in connection.execute(sql)]
    finally:
        connection.close()


def run(
    catalog_db: Path,
    results_dir: Path,
    *,
    pool_size: int = 150_000,
    decoys_per_active: int = 50,
    max_actives: int = 200,
    seed: int = 0,
) -> dict:
    """Run both arms and write ``findings.json``."""
    started = time.time()
    results_dir.mkdir(parents=True, exist_ok=True)

    # ---- actives ------------------------------------------------------------------
    raw = fetch_actives(results_dir.parent / "data" / "dhfr_activities.json")
    actives_smiles: dict[str, str] = {}
    any_dhfr_activity: set[str] = set()
    for row in raw:
        molecule_id = row.get("molecule_chembl_id")
        smiles = row.get("canonical_smiles")
        pchembl = row.get("pchembl_value")
        if not molecule_id or not smiles:
            continue
        # Anything with *any* recorded DHFR potency is barred from the decoy pool, not
        # only the compounds strong enough to count as actives here.
        any_dhfr_activity.add(molecule_id)
        if pchembl is not None and float(pchembl) >= ACTIVE_THRESHOLD:
            actives_smiles.setdefault(molecule_id, smiles)

    rng = random.Random(seed)
    active_ids = sorted(actives_smiles)
    rng.shuffle(active_ids)
    active_ids = active_ids[:max_actives]
    actives_smiles = {key: actives_smiles[key] for key in active_ids}
    print(
        f"actives: {len(active_ids)} at pChEMBL >= {ACTIVE_THRESHOLD} "
        f"({len(any_dhfr_activity)} compounds barred from the pool)",
        flush=True,
    )

    active_props = {}
    for key, smiles in list(actives_smiles.items()):
        try:
            active_props[key] = compute_properties(smiles)
        except ValueError:
            actives_smiles.pop(key)
    active_ids = sorted(active_props)

    # ---- candidate pool -----------------------------------------------------------
    catalogue = load_catalogue(catalog_db)
    print(f"catalogue: {len(catalogue):,} structures", flush=True)
    pool_props, pool_smiles = generate_pool(
        catalogue, exclude=any_dhfr_activity, size=pool_size, seed=seed
    )
    print(f"pool: {len(pool_props):,} candidates with properties", flush=True)

    # ---- the two decoy sets -------------------------------------------------------
    assigned, match_report = match_decoys(
        active_props, pool_props, n_per_active=decoys_per_active
    )
    matched_ids = sorted({key for decoys in assigned.values() for key in decoys})
    unmatched_ids = sample_unmatched(
        pool_props, n=len(matched_ids), exclude=matched_ids, seed=seed
    )
    print(
        f"decoys: {len(matched_ids):,} property-matched, {len(unmatched_ids):,} unmatched",
        flush=True,
    )

    # ---- one screen, two denominators ---------------------------------------------
    structures = dict(actives_smiles)
    for key in matched_ids + unmatched_ids:
        structures[key] = pool_smiles[key]
    prints = fingerprints(structures)

    arms = []
    for name, decoy_ids in (("property_matched", matched_ids), ("unmatched", unmatched_ids)):
        ranking = similarity_screen(active_ids, decoy_ids, prints)
        arm = Arm(
            name=name,
            n_actives=sum(ranking.labels),
            n_decoys=len(ranking) - sum(ranking.labels),
            ef_1pct=enrichment_factor(ranking.labels, 0.01),
            ef_1pct_ceiling=max_enrichment_factor(ranking.labels, 0.01),
            ef_5pct=enrichment_factor(ranking.labels, 0.05),
            bedroc_20=bedroc(ranking.labels, 20.0),
            # The direct measurement: can properties alone tell these apart?
            property_auc=property_separability(
                [active_props[key] for key in active_ids],
                [pool_props[key] for key in decoy_ids],
                properties=MATCHED_PROPERTIES,
                seed=seed,
            ),
        )
        arms.append(arm)
        print(
            f"  {name:<17} EF1% {arm.ef_1pct:6.2f} of {arm.ef_1pct_ceiling:.2f} max  "
            f"BEDROC {arm.bedroc_20:.3f}  property-only AUC {arm.property_auc:.3f}",
            flush=True,
        )

    findings = {
        "configuration": {
            "target": DHFR_TARGET,
            "active_threshold_pchembl": ACTIVE_THRESHOLD,
            "decoys_per_active": decoys_per_active,
            "pool_size": len(pool_props),
            "seed": seed,
            "screen": "ligand-based max-Tanimoto to actives, leave-one-out (ECFP4/2048)",
            "screen_not_run": (
                "Boltz-2 co-folding with an affinity head: needs a GPU this repository "
                "has not had access to. No structure-based numbers are reported."
            ),
        },
        "match_report": {
            **asdict(match_report),
            "complete": match_report.complete,
        },
        "arms": [asdict(arm) for arm in arms],
        "elapsed_seconds": round(time.time() - started, 1),
    }
    (results_dir / "findings.json").write_text(json.dumps(findings, indent=1))
    return findings

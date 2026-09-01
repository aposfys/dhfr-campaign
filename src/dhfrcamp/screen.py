"""The screen, and an honest account of which screen this is.

**The co-folding screen described in the design notes has not been run.** Boltz-2 with an
affinity head needs a GPU this repository has never had access to, and reporting numbers
from a screen that did not happen is the one thing this project is organised against.

What *is* run here is a ligand-based screen: each candidate is scored by its maximum
Tanimoto similarity to the known actives, leave-one-out. That is a real virtual screening
method, it is the baseline any structure-based method has to beat, and -- critically for
the question this repo actually asks -- it is *exactly* the method that decoy bias flatters
most. A 2D similarity screen has no idea what a protein is. If it posts a large enrichment
against one decoy set and a small one against another, the difference is a property of the
decoys, not of the method.

So the experiment here is not "how good is this screen". It is "how much of a screen's
apparent enrichment is manufactured by the choice of decoys", and a similarity screen is
the cleanest possible instrument for measuring that.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass


@dataclass(frozen=True)
class Ranking:
    """A scored, ranked candidate list."""

    ids: list[str]
    scores: list[float]
    labels: list[int]

    def __len__(self) -> int:
        return len(self.ids)


def fingerprints(smiles_by_id: Mapping[str, str], *, n_bits: int = 2048, radius: int = 2):
    """ECFP fingerprints keyed by id."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdFingerprintGenerator

    RDLogger.DisableLog("rdApp.*")
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=n_bits)
    out = {}
    for key, smiles in smiles_by_id.items():
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            out[key] = generator.GetFingerprint(mol)
    return out


def similarity_screen(
    active_ids: Sequence[str],
    decoy_ids: Sequence[str],
    prints: Mapping[str, object],
) -> Ranking:
    """Score every candidate by max Tanimoto to the actives, leave-one-out.

    Leave-one-out is not optional. Scoring an active against a reference set that contains
    itself gives it a similarity of 1.0 and guarantees perfect enrichment -- the purest
    form of the circularity this repository is about, and an easy one to commit by
    accident.
    """
    from rdkit import DataStructs

    references = [key for key in active_ids if key in prints]
    scored: list[tuple[str, float, int]] = []

    for key in references:
        others = [prints[other] for other in references if other != key]
        if not others:
            continue
        scored.append((key, max(DataStructs.BulkTanimotoSimilarity(prints[key], others)), 1))

    reference_prints = [prints[key] for key in references]
    for key in decoy_ids:
        if key not in prints:
            continue
        scored.append(
            (key, max(DataStructs.BulkTanimotoSimilarity(prints[key], reference_prints)), 0)
        )

    # Ties broken by label ascending so decoys win ties: an optimistic tie-break is a
    # silent way to inflate early enrichment.
    scored.sort(key=lambda row: (-row[1], row[2]))
    return Ranking(
        ids=[row[0] for row in scored],
        scores=[row[1] for row in scored],
        labels=[row[2] for row in scored],
    )

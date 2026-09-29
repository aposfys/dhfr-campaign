# Results

Run 2026-09-01. Every number here came from `python -m dhfrcamp.cli campaign`; nothing is
estimated.

## What was run, and what was not

**Not run: the structure-based screen.** Boltz-2 co-folding with an affinity head needs a
GPU this repository has never had access to. No structure-based enrichment is reported, and
none should be inferred from what follows.

**Run instead: a ligand-based screen, used as an instrument rather than as a result.** Each
candidate is scored by its maximum Tanimoto to the known actives, leave-one-out. This is not
here because it is a good screen. It is here because it is the cleanest available probe of
the actual question — how much of a reported enrichment is manufactured by the choice of
decoys — since one method scored against two decoy sets isolates the denominator perfectly.

## Setup

| | |
| --- | --- |
| Target | Human DHFR (`CHEMBL202`) |
| Actives | 150 compounds, a seed-0 random sample of the 714 that qualify at pChEMBL ≥ 6.0 (`--max-actives 150`) |
| Barred from the pool | 1,234 compounds — *any* recorded DHFR potency, not just the potent ones |
| Candidate pool | 60,000 ChEMBL structures with computed properties |
| Decoys | 6,440 per arm, 50 per active where the pool allowed |
| Matched on | MW, logP, HBD, HBA, rotatable bonds, formal charge |

A weak binder in a decoy set is still not a non-binder, which is why exclusion is on any
recorded activity rather than on the active threshold.

## The two arms

| Arm | EF 1% | Ceiling | BEDROC (α=20) | **Property-only AUC** |
| --- | ---: | ---: | ---: | ---: |
| Property-matched | 40.61 | 43.93 | 0.924 | **0.841** |
| Unmatched | 43.27 | 43.93 | 0.945 | **0.913** |

## Two findings, one of them inconvenient

**1. Enrichment could not distinguish the arms, because the screen was saturated.** Both
arms sit at 92–98% of the theoretical EF ceiling of 43.93. A metric pinned against its
maximum has no resolution left, so the 40.61-versus-43.27 gap is close to meaningless. This
is why `max_enrichment_factor` exists and why the ceiling is printed beside every EF: an
enrichment factor quoted without it invites exactly the comparison that cannot be made.

DHFR antifolates share the 2,4-diaminopyrimidine head group, so a leave-one-out similarity
screen finds them almost perfectly whatever the decoys look like. The instrument was too
easy for the measurement.

**2. Property matching reduced decoy bias but did not remove it.** The measurement that did
work is the direct one: train a logistic regression on the six matched properties *alone*
and see whether it can tell actives from decoys.

- Unmatched decoys: **AUC 0.913**. A model that never sees a structure gets almost all of
  the way on molecular weight and logP. This is the DUD-E failure mode, reproduced.
- Property-matched decoys: **AUC 0.841**. Better, and still nowhere near the 0.5 that
  "matched" implies.

**Per-active matching under this protocol does not produce distribution-level
indistinguishability.** Each decoy sits inside a ±25 Da, ±1 logP box around its active, and
the aggregate distributions still differ enough to be separated at 0.841. Twenty-seven of
the 150 actives could not be given a full complement of 50 decoys from a 60,000-compound
pool.

Two properties of the protocol bound how far that reading goes. `match_decoys` is
**first-fit**, accepting the first candidate inside tolerance rather than the closest, so
some of the residual separability belongs to the matcher and to the width of the box rather
than to per-active matching in principle. And the direction of the residual offsets is not
measured. The table below records per-property maxima only, so a claim that the offsets run
one way would not be supported by anything in this repository. Best-fit matching, tighter
tolerances and a wider pool are the controls that would settle it, and none of them has been
run.

The honest conclusion is that generating your own property-matched decoys is an improvement
over downloading DUD-E and, at these tolerances and with this matcher, is **not** sufficient
on its own. A decoy set should be reported with its property-only AUC, the same way a
classifier is reported with a baseline, and this one at 0.841 would not yet support a claim
that an enrichment measured against it reflects binding.

## Residual mismatch

Largest absolute per-property difference between an active and one of its assigned decoys:

| Property | Max difference | Tolerance |
| --- | ---: | ---: |
| MW | 25.0 | 25.0 |
| logP | 1.00 | 1.00 |
| HBD | 1.0 | 1.0 |
| HBA | 2.0 | 2.0 |
| Rotatable bonds | 2.0 | 2.0 |
| Formal charge | 0.0 | 0.0 |

Every property is pressed against its tolerance, which is itself the signal that the
tolerances are doing the work rather than the pool supplying genuinely close matches. These
are maxima over assigned pairs, not distributions. `findings.json` does not carry the
assigned decoy ids, so the match is summarised here rather than re-auditable from the
repository alone.

## The binding site reproduces the earlier work

`prepare.binding_site` recomputes the site from the deposited coordinates with a Biopython
KD-tree, additives stripped first. It was required to reproduce
[`protein-ligand-interaction-pymol`](https://github.com/aposfys/protein-ligand-interaction-pymol),
and it does:

```
1HFR (MOT)  7 8 9 22 30 31 34 35 59 60 61 64 67 70 115 121 136
1KMV (LII)  7 8 9 22 30 31 34    56 59 60 61 64        115 121 136
```

Both recover the conserved anchors — Ile7, Glu30, Phe31, Phe34, Val115. The difference is
the point: **Arg70 appears for MOT and not for LII**, which independently reproduces the
earlier finding that only the classical antifolate's charged glutamate tail reaches the
Arg70 subsite. That was derived in a separate codebase at a 4.5 Å cutoff with no residue
list carried over to this one, though both use a Biopython KD-tree over the same contact
criterion, so this is a reimplementation rather than a second method.

## Next

The screen this repository was designed around still has not been run. When a GPU is
available, the thing to run is the co-folding arm against **both** decoy sets, because the
property-only AUC above says the matched set is not yet a clean denominator — and a
structure-based method scored against a 0.841-separable decoy set would inherit exactly the
problem this campaign exists to avoid.

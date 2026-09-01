# Analysis

What was built, why it was built that way, and why the headline metric turned out to be the
wrong instrument.

## What could not be run

Boltz-2 co-folding with an affinity head needs a GPU this project has never had, and
REINVENT 4 generation the same. No structure-based enrichment is reported.

What ran instead is a ligand-based screen — max Tanimoto to the known actives,
leave-one-out — used as an **instrument** rather than reported as a result. One method
scored against two decoy sets isolates the denominator perfectly, which is the question the
repository is actually about.

## Design decisions, and the reasoning

**Exclusion from the decoy pool is on *any* recorded DHFR activity, not just potent
activity.** 1,234 compounds are barred though only 150 count as actives. A weak binder in a
decoy set is still not a non-binder, and scoring it as a false positive penalises the method
under test for being right.

**Decoys are assigned exclusively.** Reusing one candidate across several actives would make
the decoy set look larger and more diverse than it is — a quieter version of the bias the
module exists to remove.

**Leave-one-out is not optional.** Scoring an active against a reference set containing
itself gives similarity 1.0 and guarantees perfect enrichment. That is the purest form of
the circularity this repository is about, and an easy one to commit by accident.

**Ties break against the method.** Sorting by `(-score, label)` puts decoys first on ties.
An optimistic tie-break is a silent way to inflate early enrichment.

**The binding site is recomputed, not imported.** `prepare.binding_site` re-derives it from
the deposited coordinates with a KD-tree, additives stripped first, and is required to
reproduce `protein-ligand-interaction-pymol`. It does — including **Arg70 appearing for MOT
and not for LII**, which independently reproduces the earlier finding that only the
classical antifolate's charged glutamate tail reaches that subsite. Derived by a different
route, in a different repository, and not given to this one.

## The instrument failed, and that is reported

Enrichment factor could not distinguish the two arms. Both sit at **92–98% of a theoretical
ceiling of 43.93**, because DHFR antifolates share a 2,4-diaminopyrimidine head that a
similarity screen finds regardless of what the decoys look like. A metric pinned against its
maximum has no resolution left.

That is why `max_enrichment_factor` exists and why the ceiling prints beside every EF. An
enrichment factor quoted without its ceiling invites exactly the comparison that cannot be
made.

## The measurement that did work

Train a logistic regression on the six matched properties **alone** — no structure, no
fingerprint — and ask it to separate actives from decoys:

| Decoy set | Property-only AUC |
| --- | ---: |
| Unmatched | **0.913** |
| Property-matched | **0.841** |

**Property matching reduced decoy bias and did not remove it.** At 0.841 a model that never
sees a structure still separates actives from "matched" decoys most of the time. 0.5 is what
"matched" implies.

The mechanism is visible in the residual-mismatch table: every property is pressed against
its tolerance, and 27 of 150 actives could not be given a full 50 decoys from a
60,000-compound pool. Where the pool is thin the matcher takes what is *within tolerance*
rather than what is closest, and the offsets accumulate in one direction. Per-active
matching inside a box does not give distribution-level indistinguishability.

## What follows

Generating your own property-matched decoys beats downloading DUD-E and is not sufficient on
its own. A decoy set should be reported with its property-only AUC, the way a classifier is
reported with a baseline.

## What would change the conclusion

A GPU, and then the co-folding arm against **both** decoy sets. The 0.841 above says the
matched set is not yet a clean denominator, so a structure-based method scored against it
would inherit exactly the problem this campaign exists to avoid. Widening the pool or
tightening tolerances until property-only AUC approaches 0.5 is the prerequisite, not the
screen.

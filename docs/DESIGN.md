# dhfr-campaign — design notes

Written before the campaign was run, so the evaluation cannot be adjusted to fit the
results.

**The question:** how much enrichment is left when the decoys are matched to the actives on
every property except the one that should matter?

| | |
| --- | --- |
| **Target** | Human DHFR — 1HFR (MOT) and 1KMV (LII) |
| **Ground truth** | Binding-site residues from [`protein-ligand-interaction-pymol`](https://github.com/aposfys/protein-ligand-interaction-pymol), validated against deposited SITE records at 100% recall |
| **Screen** | Boltz-2 co-folding + affinity head; docking-guided prefilter |
| **Generation** | REINVENT 4 — scaffold decoration and R-group replacement around confirmed hits |
| **Decoys** | Property-matched, generated in-repo; **explicitly not DUD-E** |
| **Readout** | Enrichment factor and BEDROC at 1%, with a matched-decoy null |

## Why DHFR

Because the ground truth is already ours. `protein-ligand-interaction-pymol` established the
binding-site residues for both complexes with a Biopython KD-tree and validated them against
the deposited SITE records. A campaign whose target is already characterised in a previous
repo can be checked rather than believed — the screen has to rediscover contacts that were
derived independently, and failure to do so is a finding, not a footnote.

## Traps this pipeline is built to avoid

- **DUD-E and MUV decoys carry analog bias.** Actives resemble each other far more than they
  resemble the decoys, so a model can separate them on molecular weight and logP without
  learning anything about binding. This repo generates its own decoys, property-matched per
  active, and reports the property distributions alongside the enrichment so the match is
  auditable rather than asserted.
- **LIT-PCBA was the fix, and was itself audited in 2025** and found to carry leakage,
  duplication and analog redundancy. Using it uncritically as the "unbiased" option is no
  longer defensible; if it is used at all, it is used as a *third* comparison, labelled.
- **Scoring a generated molecule with the model that generated it is circular.** Generation
  (REINVENT 4) and scoring (co-folding) must not share a reward model, or the campaign
  measures agreement with itself. The scoring model is frozen before generation starts.
- **Crystallographic additives are not ligands.** The earlier DHFR work already hit this:
  a DMSO cryoprotectant counted as a binding-site residue. Structure preparation strips
  additives explicitly, from a named list, and the test suite asserts the list is applied.

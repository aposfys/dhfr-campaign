# dhfr-campaign
How much of a virtual screen's enrichment is manufactured by the choice of decoys?

[![CI](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

150 human DHFR actives against 6,440 decoys per arm, generated here rather than
downloaded and matched per active on the six properties that should be irrelevant
to binding.

```
make install
dhfrcamp prepare                             # structures, site definition
dhfrcamp campaign --catalog catalog.sqlite   # the decoy-bias experiment
dhfrcamp evaluate                            # print the table from an existing run
make test                                    # 28 tests
```

### Matching reduces decoy bias without removing it

Train a classifier on the six matched properties **alone** — no structure, no
fingerprint — and ask it to separate actives from decoys:

| Decoy set | Property-only AUC | 95% CI | EF 1% | Ceiling |
| --- | ---: | --- | ---: | ---: |
| Unmatched | **0.913** | [0.882, 0.944] | 43.27 | 43.93 |
| Property-matched | **0.841** | [0.801, 0.881] | 40.61 | 43.93 |

Both halves of the claim now carry a test. **Matching does help**: it removes 0.072 ± 0.026 of property-only AUC (z = 2.79, p = 0.005, and conservative because the arms share their actives). **Matching is nowhere near sufficient**: the matched arm sits 0.341 above chance with a lower interval bound of 0.801 (z = 16.8, p = 1.7 × 10⁻⁶³). A gap that size is not a sampling artefact.

At AUC 0.841 a model that never sees a molecule's structure still separates
actives from "matched" decoys most of the time. Per-active matching inside a
±25 Da, ±1 logP box does not produce distribution-level indistinguishability, and
27 of the 150 actives could not be given a full 50 decoys from a 60,000-compound
pool — where the pool is thin the matcher takes what is within tolerance rather
than what is closest, and the offsets accumulate one way.

The enrichment factors could not tell the arms apart at all: both sit at 92–98% of
the theoretical ceiling of 43.93, because DHFR antifolates share a
2,4-diaminopyrimidine head and a similarity screen finds them regardless of the
decoys. **An EF quoted without its ceiling invites a comparison that cannot be
made**, which is why `max_enrichment_factor` is printed beside every one.

Generating matched decoys beats downloading DUD-E and is not sufficient alone.
Report a decoy set with its property-only AUC, the way a classifier is reported
with a baseline.

### An independent reproduction

Binding-site residues recomputed here from deposited coordinates with a KD-tree
show **Arg70 for MOT and not for LII** — reproducing
[`protein-ligand-interaction-pymol`](https://github.com/aposfys/protein-ligand-interaction-pymol)'s
finding by a different route, in a different repository, and not given to this
one. Residue lists in [DESIGN.md](docs/DESIGN.md#binding-site-reproduction).

### Scope

The structure-based screen has not been run: Boltz-2 co-folding needs a GPU this
repo has never had, and `screen` and `generate` name that as the reason they are
unimplemented. What ran is a ligand-based similarity screen, used as an
*instrument* for the decoy question rather than reported as a screening result.

### Prior work

Decoy-set bias has its own literature — DUD-E and its critiques, the AVE-bias line of work,
and the general observation that a benchmark's decoys can be separated from its actives by
signals unrelated to binding. Generating property-matched decoys rather than downloading them
is a known recommendation, not a discovery here.

What this repository contributes is the pairing of a decoy set with a **property-only AUC
reported beside every enrichment factor**, and the demonstration that the recommendation is
insufficient rather than merely imperfect: per-active matching inside a ±25 Da, ±1 logP box
still leaves a classifier that never sees a structure separating actives from decoys at
0.841 [0.801, 0.881], which is 0.341 above chance at p = 1.7 × 10⁻⁶³. Matching helps —
0.072 ± 0.026, p = 0.005 — and does not come close to sufficing.

The companion observation, that an enrichment factor quoted without its ceiling invites a
comparison that cannot be made, is why `max_enrichment_factor` is printed beside every EF.

### More

- [Analysis](ANALYSIS.md) — what was done and why it was done that way
- [Results](results/RESULTS.md) — full results, including the residual mismatch table
- [Design](docs/DESIGN.md) — why DHFR, and the traps this pipeline avoids

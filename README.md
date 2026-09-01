# dhfr-campaign
How much of a virtual screen's enrichment is manufactured by the choice of decoys?

[![CI](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

Decoys are generated in this repository rather than downloaded, because the standard benchmark sets are the reason virtual screening papers report enrichment that does not survive contact with a real library. 150 human DHFR actives, 6,440 decoys per arm drawn from ChEMBL, matched per active on the six properties that should be irrelevant to binding.

> **The structure-based screen has not been run.** Boltz-2 co-folding needs a GPU this repo has never had. What ran is a ligand-based similarity screen, used as an instrument for the decoy question rather than reported as a result.

### The finding

Train a classifier on the six matched properties **alone** — no structure, no fingerprint — and ask it to tell actives from decoys:

| Decoy set | Property-only AUC | EF 1% | Ceiling |
| --- | ---: | ---: | ---: |
| Unmatched | **0.913** | 43.27 | 43.93 |
| Property-matched | **0.841** | 40.61 | 43.93 |

**Property matching reduced decoy bias and did not remove it.** At AUC 0.841 a model that never sees a molecule's structure still separates actives from "matched" decoys most of the time. Per-active matching inside a ±25 Da, ±1 logP box does not produce distribution-level indistinguishability, and 27 of the 150 actives could not be given a full 50 decoys from a 60,000-compound pool — where the pool is thin the matcher takes what is within tolerance rather than what is closest, and the offsets accumulate one way.

The enrichment factors, meanwhile, could not tell the arms apart at all: both sit at 92–98% of the theoretical ceiling of 43.93. DHFR antifolates share a 2,4-diaminopyrimidine head, so a similarity screen finds them regardless of the decoys. **An EF quoted without its ceiling invites a comparison that cannot be made**, which is why `max_enrichment_factor` is printed beside every one.

Generating your own matched decoys beats downloading DUD-E and is not sufficient on its own. Report a decoy set with its property-only AUC, the way a classifier is reported with a baseline.

### The binding site reproduces the earlier work

```
1HFR (MOT)  7 8 9 22 30 31 34 35 59 60 61 64 67 70 115 121 136
1KMV (LII)  7 8 9 22 30 31 34    56 59 60 61 64        115 121 136
```

Recomputed here from the deposited coordinates with a KD-tree, additives stripped first. **Arg70 appears for MOT and not for LII** — independently reproducing [`protein-ligand-interaction-pymol`](https://github.com/aposfys/protein-ligand-interaction-pymol)'s finding that only the classical antifolate's charged glutamate tail reaches that subsite. Derived by a different route, in a different repository, and not given to this one.

### Running it

```
make install
python3 -m dhfrcamp.cli prepare                    # structures, site definition
python3 -m dhfrcamp.cli campaign --catalog catalog.sqlite
make test
```

`catalog.sqlite` is the ChEMBL structure catalogue built by [`chem-explorer`](https://github.com/aposfys/chem-explorer)'s `tools/build_catalog.py`. The `screen` and `generate` subcommands remain unimplemented and say so — they are the GPU half.

### Layout

```
src/dhfrcamp/
  prepare.py    structure retrieval, additive stripping, KD-tree site definition
  decoys.py     property matching, the candidate pool, the unmatched comparison arm
  screen.py     the ligand-based screen, and why it is the instrument here
  campaign.py   one screen, two denominators
  evaluate.py   enrichment factor, its ceiling, BEDROC, property separability
  cli.py        `python -m dhfrcamp.cli`
```

21 tests.

### More

- [Full results, including the residual mismatch table](results/RESULTS.md)
- [Why DHFR, and the traps the pipeline is built to avoid](docs/DESIGN.md)

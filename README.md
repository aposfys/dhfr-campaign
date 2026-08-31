# dhfr-campaign
A retrospective structure-based virtual screen against human DHFR, with an honest denominator.

[![CI](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** Interfaces and tests are in place; no campaign has been run. Every number this repo will report is a target to earn, not a result.

Decoys are generated in this repository rather than downloaded from a benchmark set, because the standard benchmark sets are the reason virtual screening papers report enrichment that does not survive contact with a real library. The question: how much enrichment is left when the decoys are matched to the actives on every property except the one that should matter?

Binding-site ground truth comes from [`protein-ligand-interaction-pymol`](https://github.com/aposfys/protein-ligand-interaction-pymol), where it was derived independently and validated against the deposited SITE records — so the screen has to rediscover contacts it was not given, and failure to do so is a finding.

### Running it
```
make install
make data        # structures, known actives from ChEMBL, decoy generation
make screen      # the expensive step; resumable, cached per ligand
make analysis    # enrichment against the matched-decoy null
make test
```
The CLI subcommands are stubs at present and exit with "not implemented".

### Layout
```
src/dhfrcamp/
  prepare.py    structure retrieval, additive stripping, protonation, box definition
  decoys.py     property-matched decoy generation and the match report
  evaluate.py   enrichment factor, BEDROC, and the matched-decoy null
  cli.py        `python -m dhfrcamp.cli`
```
Planned: `screen.py` (co-folding + affinity scoring), `generate.py` (REINVENT 4 expansion).

Curation of the actives reuses `chembench` from [`chem-benchmark-audit`](https://github.com/aposfys/chem-benchmark-audit) rather than reimplementing standardisation.

### Design notes
[Why DHFR, the planned screen, and the traps the pipeline is built to avoid](docs/DESIGN.md)

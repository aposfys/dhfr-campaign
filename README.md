# dhfr-campaign
How much of a virtual screen's enrichment is manufactured by the choice of decoys?

[![CI](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/dhfr-campaign/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

150 human DHFR actives, a seed-0 sample of the 714 that qualify at pChEMBL 6.0,
against 6,440 decoys per arm, generated here rather than downloaded and matched per
active on the six properties that should be irrelevant to binding.

```
make install
dhfrcamp prepare        # structures, additives stripped, site definition
dhfrcamp evaluate       # the table below, from the committed run
make test               # 36 tests

# rerunning the experiment needs a ChEMBL structure catalogue, built locally
curl -O https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/chembl_36_chemreps.txt.gz
python3 tools/build_catalog.py chembl_36_chemreps.txt.gz data/catalog.sqlite
dhfrcamp campaign       # one screen, two decoy sets
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
actives from "matched" decoys most of the time. The matcher is first-fit rather
than nearest-neighbour, so it takes the first candidate inside a ±25 Da, ±1 logP
box and not the closest one, and 27 of the 150 actives could not be given a full
50 decoys from a 60,000-compound pool. So 0.841 bounds this protocol, not
per-active matching in general. Where the residual offsets sit is not measured
here, only their per-property maxima, and each of those lands on its tolerance.
Best-fit matching, tighter tolerances and a wider pool are the controls that would
separate the box from the matcher, and none of them has been run.

The enrichment factors could not tell the arms apart at all: both sit at 92–98% of
the theoretical ceiling of 43.93, because DHFR antifolates share a
2,4-diaminopyrimidine head and a similarity screen finds them regardless of the
decoys. **An EF quoted without its ceiling invites a comparison that cannot be
made**, which is why `max_enrichment_factor` is printed beside every one.

Generating matched decoys beats downloading DUD-E and is not sufficient alone.
Report a decoy set with its property-only AUC, the way a classifier is reported
with a baseline.

### Reproducing the earlier binding site

Binding-site residues recomputed here from deposited coordinates with a KD-tree
show **Arg70 for MOT and not for LII**, reproducing
[`protein-ligand-interaction-pymol`](https://github.com/aposfys/protein-ligand-interaction-pymol)'s
finding in a separate codebase at a 4.5 Å cutoff, with no residue list carried
over from it. Both use a Biopython KD-tree, so this is a reimplementation rather
than a second method. Residue lists in
[DESIGN.md](docs/DESIGN.md#binding-site-reproduction).

### Scope

The structure-based screen has not been run: Boltz-2 co-folding needs a GPU this
repo has never had, and `screen` and `generate` name that as the reason they are
unimplemented. What ran is a ligand-based similarity screen, used as an
*instrument* for the decoy question rather than reported as a screening result.

### Bundled data

`data/dhfr_activities.json` is a ChEMBL activity dump for `CHEMBL202`, fetched from the
ChEMBL web API on 2026-09-01 and redistributed here under CC BY-SA 3.0 (ChEMBL, EMBL-EBI).
`data/pdb/1hfr.pdb` and `data/pdb/1kmv.pdb` are RCSB PDB entries 1HFR and 1KMV.
`catalog.sqlite` is built locally from a ChEMBL chemreps dump and is not committed, and the
release it was built from was not recorded, so the decoy pool is not byte-reproducible.

### More

- [Analysis](ANALYSIS.md) — what was done, why it was done that way, and the prior work
- [Results](results/RESULTS.md) — full results, including the residual mismatch table
- [Design](docs/DESIGN.md) — why DHFR, and the traps this pipeline avoids

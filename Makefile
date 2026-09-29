.PHONY: install data catalog decoys campaign analysis screen test clean clean-data all

PYTHON ?= python3

## Print the result table from the committed run. Needs no catalogue and no GPU, so
## `make` works on a fresh clone.
all: analysis

## Install the package plus dev tooling.
install:
	$(PYTHON) -m pip install -e ".[dev]"

## Structures (1HFR, 1KMV), additives stripped, binding site defined
data:
	$(PYTHON) -m dhfrcamp.cli prepare

## The id -> structure catalogue the decoy pool is drawn from. CHEMREPS is a
## chembl_NN_chemreps.txt.gz dump from the EBI FTP site.
##   make catalog CHEMREPS=chembl_36_chemreps.txt.gz
catalog:
	$(PYTHON) tools/build_catalog.py $(CHEMREPS) data/catalog.sqlite

## Property-matched decoys, one set per active, plus the match report
decoys:
	$(PYTHON) -m dhfrcamp.cli decoys

## The decoy-bias experiment, one screen against two decoy sets. Needs
## data/catalog.sqlite from the catalog target, and rewrites results/findings.json.
campaign:
	$(PYTHON) -m dhfrcamp.cli campaign

## Enrichment factors with their ceiling, BEDROC, and property-only AUC with
## uncertainty, read from an existing run
analysis:
	$(PYTHON) -m dhfrcamp.cli evaluate

## Not implemented, and kept out of the default chain for that reason. Boltz-2
## co-folding needs a GPU this project has not had, so this target exits 1.
screen:
	$(PYTHON) -m dhfrcamp.cli screen

test:
	$(PYTHON) -m pytest -q

## Caches and generated models only. The committed findings.json and RESULTS.md stay.
clean:
	rm -rf results/models
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete the locally built catalogue and any cached structure downloads
clean-data: clean
	rm -f data/*.cif data/*.pdb data/*.csv data/*.smi data/*.sqlite

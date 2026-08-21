.PHONY: install data decoys screen analysis test clean clean-data all

PYTHON ?= python3

all: analysis

## Install the package plus dev tooling. Co-folding and generation backends are
## optional extras:  pip install -e ".[chem,folding,generate]"
install:
	$(PYTHON) -m pip install -e ".[dev]"

## Structures (1HFR, 1KMV), known DHFR actives from ChEMBL, and structure preparation
data:
	$(PYTHON) -m dhfrcamp.cli prepare

## Property-matched decoys, one set per active, plus the match report
decoys: data
	$(PYTHON) -m dhfrcamp.cli decoys

## The expensive step. Resumable: already-scored ligands are skipped, so a
## re-run after an interruption continues rather than restarting.
screen: decoys
	$(PYTHON) -m dhfrcamp.cli screen

## Enrichment factor and BEDROC against the matched-decoy null
analysis: screen
	$(PYTHON) -m dhfrcamp.cli evaluate

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/*
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete cached structures, actives and generated decoys
clean-data: clean
	rm -f data/*.cif data/*.pdb data/*.csv data/*.smi

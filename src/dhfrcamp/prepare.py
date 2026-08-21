"""Structure retrieval and preparation for the two DHFR complexes.

The binding-site definition is not recomputed from scratch. It is imported from the
earlier `protein-ligand-interaction-pymol` work, where it was derived with a Biopython
KD-tree and validated against the deposited SITE records at 100% recall, and this module
must reproduce it. Failing to reproduce it is a finding.
"""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

from dhfrcamp import CRYSTALLOGRAPHIC_ADDITIVES

#: The two complexes: PDB id -> ligand HET code.
COMPLEXES = {"1HFR": "MOT", "1KMV": "LII"}

#: Contact distance used to define the binding site, matching the earlier work.
CONTACT_CUTOFF_ANGSTROM = 4.5


def strip_additives(het_codes: Iterable[str]) -> list[str]:
    """Drop crystallographic additives, keeping genuine ligands.

    This exists because the earlier DHFR analysis counted a DMSO cryoprotectant as a
    binding-site residue. Stripping happens before the site is defined, never after.
    """
    return [code for code in het_codes if code.upper() not in CRYSTALLOGRAPHIC_ADDITIVES]


def fetch_structure(pdb_id: str, out_dir: Path) -> Path:
    """Download one structure to ``out_dir``, cached."""
    raise NotImplementedError("milestone 1: structure retrieval")


def binding_site(pdb_id: str, *, cutoff: float = CONTACT_CUTOFF_ANGSTROM) -> set[int]:
    """Residue numbers within ``cutoff`` of the ligand, additives already stripped."""
    raise NotImplementedError("milestone 1: KD-tree contact detection")

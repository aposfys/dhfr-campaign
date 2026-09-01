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


RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"


def fetch_structure(pdb_id: str, out_dir: Path) -> Path:
    """Download one structure to ``out_dir``, cached."""
    import urllib.request

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{pdb_id.lower()}.pdb"
    if path.exists():
        return path
    with urllib.request.urlopen(
        RCSB_URL.format(pdb_id=pdb_id.upper()), timeout=120
    ) as response:
        path.write_bytes(response.read())
    return path


def binding_site(
    pdb_id: str,
    *,
    cutoff: float = CONTACT_CUTOFF_ANGSTROM,
    data_dir: Path | None = None,
) -> set[int]:
    """Residue numbers within ``cutoff`` of the ligand, additives already stripped.

    Reproduces the definition from `protein-ligand-interaction-pymol`: heavy-atom contacts
    only, waters and the NADPH cofactor removed, and crystallographic additives dropped
    *before* the site is computed rather than filtered out of the answer afterwards.
    """
    from Bio.PDB import NeighborSearch, PDBParser

    ligand_code = COMPLEXES[pdb_id.upper()]
    path = fetch_structure(pdb_id, data_dir or Path("data/pdb"))

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, str(path))
    model = next(iter(structure))

    ligand_atoms = []
    environment_atoms = []
    for residue in model.get_residues():
        het_flag, _, _ = residue.get_id()
        name = residue.get_resname().strip().upper()
        atoms = [atom for atom in residue if atom.element != "H"]
        if name == ligand_code:
            ligand_atoms.extend(atoms)
        elif het_flag == " " and name not in CRYSTALLOGRAPHIC_ADDITIVES:
            # Standard residues only: the cofactor and every additive are excluded from
            # the environment, so neither can appear in the site definition.
            environment_atoms.extend(atoms)

    if not ligand_atoms:
        raise ValueError(f"{pdb_id}: ligand {ligand_code} not found")

    search = NeighborSearch(environment_atoms)
    residues: set[int] = set()
    for atom in ligand_atoms:
        for neighbour in search.search(atom.coord, cutoff, level="R"):
            residues.add(int(neighbour.get_id()[1]))
    return residues

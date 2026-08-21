"""Retrospective structure-based virtual screening against human DHFR.

The repository is organised so that the *denominator* of any enrichment claim is
constructed in this codebase and can be inspected: decoys are generated here, matched
here, and reported here. Screening results are meaningless without them.
"""

__version__ = "0.1.0"

#: Crystallographic additives that are not ligands. The earlier DHFR work counted a DMSO
#: cryoprotectant as a binding-site residue; this list exists so that cannot recur, and
#: the test suite asserts it is applied before any site definition.
CRYSTALLOGRAPHIC_ADDITIVES = frozenset(
    {
        "HOH",  # water
        "DMS",  # dimethyl sulfoxide
        "GOL",  # glycerol
        "EDO",  # ethylene glycol
        "PEG",
        "SO4",
        "PO4",
        "CL",
        "NA",
        "K",
        "MG",
        "ZN",
        "ACT",  # acetate
        "MPD",  # 2-methyl-2,4-pentanediol
        "TRS",  # tris buffer
    }
)

__all__ = ["CRYSTALLOGRAPHIC_ADDITIVES", "__version__"]

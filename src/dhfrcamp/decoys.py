"""Property-matched decoy generation, and the report that makes the match auditable.

The argument of this repo is that enrichment is only meaningful relative to a denominator
you can inspect. DUD-E and MUV decoys carry analog bias: actives resemble each other far
more than they resemble the decoys, so a classifier can separate them on molecular weight
and logP alone. Decoys here are matched per active on physicochemical properties that
should be irrelevant to binding, and the residual mismatch is reported rather than assumed
away.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass

#: Properties matched between an active and its decoys. Formal charge and rotatable bond
#: count are included because they are the two most common leaks in published decoy sets.
MATCHED_PROPERTIES = ("mw", "logp", "hbd", "hba", "rotatable_bonds", "formal_charge")

#: Default per-property tolerance for calling an active and a candidate "matched".
DEFAULT_TOLERANCES: dict[str, float] = {
    "mw": 25.0,
    "logp": 1.0,
    "hbd": 1.0,
    "hba": 2.0,
    "rotatable_bonds": 2.0,
    "formal_charge": 0.0,
}

Properties = Mapping[str, float]


@dataclass(frozen=True)
class MatchReport:
    """How well the decoy set actually matches the actives, per property."""

    n_actives: int
    n_decoys: int
    max_abs_difference: dict[str, float]
    unmatched_actives: list[str]

    @property
    def complete(self) -> bool:
        """True when every active received its full complement of decoys."""
        return not self.unmatched_actives


def is_match(
    active: Properties, candidate: Properties, tolerances: Mapping[str, float] | None = None
) -> bool:
    """Whether a candidate falls inside tolerance of an active on every matched property."""
    limits = dict(DEFAULT_TOLERANCES if tolerances is None else tolerances)
    for prop in MATCHED_PROPERTIES:
        if prop not in active or prop not in candidate:
            raise KeyError(f"property {prop!r} missing from active or candidate")
        if abs(active[prop] - candidate[prop]) > limits[prop]:
            return False
    return True


def match_decoys(
    actives: Mapping[str, Properties],
    pool: Mapping[str, Properties],
    *,
    n_per_active: int = 50,
    tolerances: Mapping[str, float] | None = None,
) -> tuple[dict[str, list[str]], MatchReport]:
    """Assign decoys from ``pool`` to each active, without reusing a decoy.

    Decoys are assigned exclusively: reusing one candidate across several actives would
    make the decoy set look larger and more diverse than it is, which is a quieter version
    of the same bias this module exists to remove.
    """
    assigned: dict[str, list[str]] = {}
    taken: set[str] = set()
    worst: dict[str, float] = dict.fromkeys(MATCHED_PROPERTIES, 0.0)

    for active_id, active_props in actives.items():
        chosen: list[str] = []
        for candidate_id, candidate_props in pool.items():
            if candidate_id in taken or len(chosen) >= n_per_active:
                continue
            if is_match(active_props, candidate_props, tolerances):
                chosen.append(candidate_id)
                taken.add(candidate_id)
                for prop in MATCHED_PROPERTIES:
                    delta = abs(active_props[prop] - candidate_props[prop])
                    worst[prop] = max(worst[prop], delta)
        assigned[active_id] = chosen

    unmatched = [key for key, decoys in assigned.items() if len(decoys) < n_per_active]
    report = MatchReport(
        n_actives=len(actives),
        n_decoys=sum(len(decoys) for decoys in assigned.values()),
        max_abs_difference=worst,
        unmatched_actives=unmatched,
    )
    return assigned, report


def generate_pool(actives: Sequence[str], *, size: int) -> dict[str, Properties]:
    """Draw a candidate pool from a purchasable library, excluding known DHFR binders.

    Excluding known binders matters: a decoy that is actually an unannotated active is
    scored as a false positive and silently penalises the method under test.
    """
    raise NotImplementedError("milestone 2: candidate pool from a purchasable library")

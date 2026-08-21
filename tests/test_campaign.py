"""Invariants that keep a screening result honest.

Each test guards a failure mode that would inflate the headline number rather than
crash the pipeline -- which is the only kind of bug that matters here.
"""

from __future__ import annotations

import pytest

from dhfrcamp import CRYSTALLOGRAPHIC_ADDITIVES, decoys, evaluate, prepare


def test_dmso_is_stripped_before_the_site_is_defined() -> None:
    """The exact trap the earlier DHFR analysis hit: DMSO counted as a site residue."""
    assert "DMS" in CRYSTALLOGRAPHIC_ADDITIVES
    assert prepare.strip_additives(["MOT", "DMS", "HOH", "GOL"]) == ["MOT"]


def test_genuine_ligands_survive_stripping() -> None:
    assert prepare.strip_additives(["MOT", "LII"]) == ["MOT", "LII"]


def test_perfect_ranking_beats_random_by_the_active_fraction() -> None:
    """100 compounds, 10 actives, all ranked first: EF at 10% is the maximum, 10x."""
    labels = [1] * 10 + [0] * 90
    assert evaluate.enrichment_factor(labels, fraction=0.10) == pytest.approx(10.0)


def test_random_ranking_scores_one() -> None:
    """Enrichment of 1.0 means no better than picking blind."""
    labels = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0] * 10
    assert evaluate.enrichment_factor(labels, fraction=0.10) == pytest.approx(1.0)


def test_bedroc_prefers_actives_at_the_top() -> None:
    """A method that finds the same actives later must score strictly lower."""
    early = [1] * 5 + [0] * 95
    late = [0] * 95 + [1] * 5
    assert evaluate.bedroc(early) > evaluate.bedroc(late)


def test_screen_without_actives_is_an_error_not_a_zero() -> None:
    with pytest.raises(ValueError):
        evaluate.enrichment_factor([0, 0, 0, 0])


def test_decoys_are_never_shared_between_actives() -> None:
    """Reusing a decoy inflates apparent diversity of the denominator."""
    props = {
        "mw": 300.0,
        "logp": 2.0,
        "hbd": 1.0,
        "hba": 3.0,
        "rotatable_bonds": 4.0,
        "formal_charge": 0.0,
    }
    actives = {"a1": props, "a2": props}
    pool = {f"d{i}": dict(props) for i in range(3)}
    assigned, report = decoys.match_decoys(actives, pool, n_per_active=2)
    all_assigned = [decoy for chosen in assigned.values() for decoy in chosen]
    assert len(all_assigned) == len(set(all_assigned))
    assert report.n_decoys == 3


def test_incomplete_match_is_reported_not_hidden() -> None:
    """An active that could not be matched must surface, not silently shrink the null."""
    props = {
        "mw": 300.0,
        "logp": 2.0,
        "hbd": 1.0,
        "hba": 3.0,
        "rotatable_bonds": 4.0,
        "formal_charge": 0.0,
    }
    far = dict(props, mw=900.0)
    assigned, report = decoys.match_decoys({"a1": props}, {"d1": far}, n_per_active=5)
    assert assigned["a1"] == []
    assert report.unmatched_actives == ["a1"]
    assert not report.complete

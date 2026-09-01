"""Metrics, the enrichment ceiling, and the direct measurement of decoy bias."""

from __future__ import annotations

import pytest

from dhfrcamp.decoys import MATCHED_PROPERTIES, is_match, match_decoys, sample_unmatched
from dhfrcamp.evaluate import (
    bedroc,
    enrichment_factor,
    max_enrichment_factor,
    property_separability,
)


def props(mw, logp, hbd=1.0, hba=2.0, rot=2.0, charge=0.0):
    return {
        "mw": mw,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot,
        "formal_charge": charge,
    }


def test_enrichment_of_a_perfect_ranking_hits_its_ceiling():
    labels = [1] * 10 + [0] * 990
    assert enrichment_factor(labels, 0.01) == pytest.approx(
        max_enrichment_factor(labels, 0.01)
    )


def test_enrichment_of_a_random_ranking_is_about_one():
    labels = ([1] + [0] * 9) * 100
    assert enrichment_factor(labels, 0.1) == pytest.approx(1.0, abs=0.15)


def test_the_ceiling_is_below_the_naive_reciprocal_when_actives_are_plentiful():
    """With more actives than top-slice slots, 1/active_fraction is unreachable."""
    labels = [1] * 50 + [0] * 50
    # 1% of 100 is 1 slot, but the active fraction is 0.5, so EF can reach at most 2.0.
    assert max_enrichment_factor(labels, 0.01) == pytest.approx(2.0)


def test_bedroc_rewards_early_recognition_over_late():
    early = [1] * 5 + [0] * 95
    late = [0] * 95 + [1] * 5
    assert bedroc(early) > 0.9
    assert bedroc(late) < 0.1


@pytest.mark.parametrize("bad", [[], [0, 0, 0], [1, 2, 0]])
def test_degenerate_rankings_are_refused(bad):
    with pytest.raises(ValueError):
        enrichment_factor(bad)


def test_property_separability_is_chance_for_identical_distributions():
    # Drawn from one distribution and split at random -- not sliced out of a sorted
    # sequence, which would be perfectly separable and would test nothing.
    import random

    rng = random.Random(0)
    drawn = [props(rng.gauss(300, 40), rng.gauss(2.0, 1.0)) for _ in range(200)]
    rng.shuffle(drawn)
    auc = property_separability(drawn[:100], drawn[100:], properties=MATCHED_PROPERTIES)
    assert 0.35 < auc < 0.65, f"identical distributions should be near chance, got {auc}"


def test_property_separability_detects_an_obviously_biased_decoy_set():
    """The failure DUD-E is criticised for: separable on properties alone."""
    actives = [props(500.0 + i, 4.0) for i in range(40)]
    decoys = [props(200.0 + i, 0.5) for i in range(40)]
    assert property_separability(actives, decoys, properties=MATCHED_PROPERTIES) > 0.95


def test_matching_refuses_a_candidate_outside_any_tolerance():
    active = props(300.0, 2.0)
    assert is_match(active, props(310.0, 2.5))
    assert not is_match(active, props(400.0, 2.0))  # mw
    assert not is_match(active, props(300.0, 5.0))  # logp
    assert not is_match(active, props(300.0, 2.0, charge=1.0))  # formal charge, zero tolerance


def test_a_decoy_is_never_assigned_to_two_actives():
    """Reuse would make the decoy set look larger and more diverse than it is."""
    actives = {f"a{i}": props(300.0, 2.0) for i in range(3)}
    pool = {f"d{i}": props(300.0, 2.0) for i in range(10)}
    assigned, report = match_decoys(actives, pool, n_per_active=5)
    everything = [key for decoys in assigned.values() for key in decoys]
    assert len(everything) == len(set(everything))
    assert report.n_decoys == len(everything)


def test_the_match_report_admits_when_actives_went_short():
    actives = {"a0": props(300.0, 2.0), "a1": props(900.0, 9.0)}
    pool = {f"d{i}": props(300.0, 2.0) for i in range(5)}
    _, report = match_decoys(actives, pool, n_per_active=5)
    assert not report.complete
    assert "a1" in report.unmatched_actives


def test_unmatched_sampling_excludes_what_it_is_told_to():
    pool = {f"d{i}": props(300.0, 2.0) for i in range(20)}
    picked = sample_unmatched(pool, n=5, exclude=[f"d{i}" for i in range(15)])
    assert len(picked) == 5
    assert all(int(key[1:]) >= 15 for key in picked)

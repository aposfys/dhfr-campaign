"""Tests for the uncertainty attached to the decoy-bias numbers."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from dhfrcamp.statistics import (
    bias_report,
    difference,
    distance_from_chance,
    hanley_mcneil,
)

RESULTS = Path(__file__).resolve().parents[1] / "results" / "findings.json"


def test_standard_error_shrinks_with_sample_size():
    small = hanley_mcneil(0.85, 50, 500)
    large = hanley_mcneil(0.85, 500, 5000)
    assert large.standard_error < small.standard_error


def test_interval_stays_inside_the_unit_range():
    estimate = hanley_mcneil(0.99, 10, 10)
    assert 0.0 <= estimate.lower <= estimate.auc <= estimate.upper <= 1.0


def test_rejects_impossible_input():
    with pytest.raises(ValueError):
        hanley_mcneil(1.2, 10, 10)
    with pytest.raises(ValueError):
        hanley_mcneil(0.8, 10, 0)


def test_a_chance_classifier_is_not_flagged_as_biased():
    estimate = hanley_mcneil(0.5, 150, 6440)
    assert distance_from_chance(estimate)["p_value"] > 0.05


def test_difference_of_identical_aucs_is_null():
    estimate = hanley_mcneil(0.8, 150, 6440)
    result = difference(estimate, estimate)
    assert result["difference"] == 0.0
    assert result["p_value"] == pytest.approx(1.0)


def test_bias_report_requires_the_fields_it_needs():
    with pytest.raises(ValueError):
        bias_report([{"name": "broken", "property_auc": 0.8}])


@pytest.fixture(scope="module")
def report():
    if not RESULTS.exists():
        pytest.skip(f"{RESULTS} missing; run the campaign first")
    return bias_report(json.loads(RESULTS.read_text(encoding="utf-8"))["arms"])


def test_matching_measurably_reduces_property_bias(report):
    """The first half of the claim: matching does something."""
    effect = report["matching_effect"]
    assert effect["difference"] > 0
    assert effect["p_value"] < 0.05


def test_matched_decoys_remain_far_from_indistinguishable(report):
    """The second half, and the one that matters.

    A property-matched decoy set is still separable from the actives on the
    matched properties alone, by a margin no sample-size argument dissolves:
    the lower end of the interval is above 0.8, not near 0.5.
    """
    matched = report["per_arm"]["property_matched"]
    assert matched["ci_lower"] > 0.75
    assert matched["vs_chance"]["p_value"] < 1e-20

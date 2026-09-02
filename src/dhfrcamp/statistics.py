"""Uncertainty for the decoy-bias numbers.

The campaign's central claim is comparative: property matching *reduces* decoy
bias but does not remove it. Both halves of that sentence are claims about a
difference, and neither is safe from point estimates alone -- the property-only
AUCs are cross-validated estimates from 150 actives, and a gap of 0.07 between
two such numbers could be nothing.

Two questions are answered here, from the counts the campaign already records,
so no molecule is re-fetched and no classifier is refitted:

1. **Did matching change anything?** The difference between the two arms'
   property-only AUCs, with a standard error.
2. **Is the matched arm still biased?** Its distance from chance, which is the
   question that decides whether "property-matched" means "unbiased".

AUC standard errors use the Hanley--McNeil closed form (Hanley & McNeil,
*Radiology* 1982), which needs only the AUC and the two class sizes.

The comparison in (1) treats the two arms as independent. They share their 150
actives, so the true standard error of the difference is smaller than the one
reported here and the test is **conservative** -- it understates the evidence
that matching helped, which is the safe direction for this claim.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

Z_95 = 1.959963984540054


@dataclass(frozen=True)
class AucEstimate:
    """An AUC with its Hanley--McNeil standard error and interval."""

    auc: float
    standard_error: float
    lower: float
    upper: float

    def as_dict(self) -> dict[str, float]:
        return {
            "auc": round(self.auc, 4),
            "standard_error": round(self.standard_error, 4),
            "ci_lower": round(self.lower, 4),
            "ci_upper": round(self.upper, 4),
        }


def _normal_two_sided_p(z: float) -> float:
    """Two-sided p-value for a standard normal deviate."""
    return math.erfc(abs(z) / math.sqrt(2.0))


def hanley_mcneil(auc: float, n_positive: int, n_negative: int) -> AucEstimate:
    """Standard error and 95% interval for a ROC AUC.

    Raises:
        ValueError: if the AUC is outside [0, 1] or either class is empty.
    """
    if not 0.0 <= auc <= 1.0:
        raise ValueError(f"auc {auc} outside [0, 1]")
    if n_positive <= 0 or n_negative <= 0:
        raise ValueError("both classes must be non-empty")

    q1 = auc / (2 - auc)
    q2 = 2 * auc**2 / (1 + auc)
    variance = (
        auc * (1 - auc) + (n_positive - 1) * (q1 - auc**2) + (n_negative - 1) * (q2 - auc**2)
    ) / (n_positive * n_negative)
    error = math.sqrt(max(variance, 0.0))
    return AucEstimate(
        auc=auc,
        standard_error=error,
        lower=max(0.0, auc - Z_95 * error),
        upper=min(1.0, auc + Z_95 * error),
    )


def difference(first: AucEstimate, second: AucEstimate) -> dict[str, float]:
    """Difference between two AUCs, treating them as independent.

    Conservative when the two arms share actives, as they do here.
    """
    delta = first.auc - second.auc
    error = math.sqrt(first.standard_error**2 + second.standard_error**2)
    z = delta / error if error else 0.0
    return {
        "difference": round(delta, 4),
        "standard_error": round(error, 4),
        "z": round(z, 3),
        "p_value": _normal_two_sided_p(z),
    }


def distance_from_chance(estimate: AucEstimate) -> dict[str, float]:
    """How far an AUC sits from 0.5.

    For a property-only classifier this is the bias measurement itself: at
    chance the decoys are indistinguishable from the actives on the matched
    properties, and anything above it is separable without seeing a structure.
    """
    z = (estimate.auc - 0.5) / estimate.standard_error if estimate.standard_error else 0.0
    return {
        "excess_over_chance": round(estimate.auc - 0.5, 4),
        "z": round(z, 3),
        "p_value": _normal_two_sided_p(z),
        "ci_lower": round(estimate.lower, 4),
    }


def bias_report(arms: list[dict[str, Any]]) -> dict[str, Any]:
    """Attach uncertainty to every arm's property-only AUC.

    Raises:
        ValueError: if an arm lacks the fields the estimate needs.
    """
    estimates: dict[str, AucEstimate] = {}
    for arm in arms:
        try:
            name = str(arm["name"])
            estimates[name] = hanley_mcneil(
                float(arm["property_auc"]),
                int(arm["n_actives"]),
                int(arm["n_decoys"]),
            )
        except KeyError as error:
            raise ValueError(f"arm is missing {error}") from error

    report: dict[str, Any] = {
        "method": {
            "auc_interval": "Hanley-McNeil, 95%",
            "caveat": (
                "The arms share their actives, so the between-arm test is "
                "conservative: the true standard error of the difference is "
                "smaller than the independent one used here."
            ),
        },
        "per_arm": {
            name: {
                **estimate.as_dict(),
                "vs_chance": distance_from_chance(estimate),
            }
            for name, estimate in estimates.items()
        },
    }

    if "unmatched" in estimates and "property_matched" in estimates:
        report["matching_effect"] = difference(
            estimates["unmatched"], estimates["property_matched"]
        )
    return report

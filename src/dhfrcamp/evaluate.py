"""Early-recognition metrics for a ranked screen.

A virtual screen is judged on the top of the list, because that is all anyone will ever
buy. Plain ROC AUC hides this: a method can score well while putting its actives at rank
5,000. Enrichment factor and BEDROC both weight the top, and both are reported here
against the matched-decoy null rather than against a downloaded benchmark set.
"""

from __future__ import annotations

import math
from collections.abc import Sequence


def enrichment_factor(labels_by_rank: Sequence[int], fraction: float = 0.01) -> float:
    """Enrichment factor at the top ``fraction`` of a ranked list.

    ``labels_by_rank`` is 1 for an active and 0 for a decoy, ordered best score first.
    A value of 1.0 means the method did exactly as well as picking at random.
    """
    _check_labels(labels_by_rank)
    if not 0.0 < fraction <= 1.0:
        raise ValueError(f"fraction must lie in (0, 1], got {fraction}")

    n_total = len(labels_by_rank)
    n_actives = sum(labels_by_rank)
    if n_actives == 0:
        raise ValueError("cannot compute enrichment without actives")

    n_top = max(1, round(n_total * fraction))
    hits_in_top = sum(labels_by_rank[:n_top])
    return (hits_in_top / n_top) / (n_actives / n_total)


def bedroc(labels_by_rank: Sequence[int], alpha: float = 20.0) -> float:
    """Boltzmann-enhanced discrimination of ROC (Truchon & Bayly, 2007).

    ``alpha`` sets how sharply the top of the list is weighted; 20.0 concentrates roughly
    80% of the weight in the top 8%. Returns a value in [0, 1], where the random
    expectation is not 0.5 but depends on the active fraction -- which is exactly why the
    matched-decoy null is reported alongside it.
    """
    _check_labels(labels_by_rank)
    n_total = len(labels_by_rank)
    n_actives = sum(labels_by_rank)
    if n_actives == 0:
        raise ValueError("cannot compute BEDROC without actives")
    if n_actives == n_total:
        return 1.0

    ratio = n_actives / n_total
    total = sum(
        math.exp(-alpha * (rank + 1) / n_total)
        for rank, label in enumerate(labels_by_rank)
        if label == 1
    )
    random_sum = (
        n_actives * (1 - math.exp(-alpha)) / (n_total * (math.exp(alpha / n_total) - 1))
    )
    scaled = total / random_sum
    factor = (
        ratio
        * math.sinh(alpha / 2)
        / (math.cosh(alpha / 2) - math.cosh(alpha / 2 - alpha * ratio))
    )
    offset = 1 / (1 - math.exp(alpha * (1 - ratio)))
    return scaled * factor + offset


def _check_labels(labels: Sequence[int]) -> None:
    if not labels:
        raise ValueError("cannot score an empty ranking")
    if any(label not in (0, 1) for label in labels):
        raise ValueError("labels must be 1 for actives and 0 for decoys")

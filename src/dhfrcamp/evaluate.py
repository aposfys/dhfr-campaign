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


def max_enrichment_factor(labels_by_rank: Sequence[int], fraction: float = 0.01) -> float:
    """The largest enrichment factor achievable at this fraction and active count.

    A perfect ranking cannot exceed ``1 / active_fraction``, and when the top slice is
    smaller than the number of actives the ceiling is lower still. Reporting EF without
    its ceiling invites comparing two numbers that were never on the same scale -- and a
    screen sitting at 95% of maximum has no room left to distinguish anything.
    """
    _check_labels(labels_by_rank)
    n_total = len(labels_by_rank)
    n_actives = sum(labels_by_rank)
    if n_actives == 0:
        raise ValueError("cannot compute enrichment without actives")
    n_top = max(1, round(n_total * fraction))
    best_hits = min(n_top, n_actives)
    return (best_hits / n_top) / (n_actives / n_total)


def property_separability(
    active_properties: Sequence[dict[str, float]],
    decoy_properties: Sequence[dict[str, float]],
    *,
    properties: Sequence[str],
    seed: int = 0,
) -> float:
    """Cross-validated ROC AUC for telling actives from decoys on properties alone.

    This is the direct measurement of decoy bias, and the one that matters. If a model can
    separate actives from decoys using only molecular weight, logP and hydrogen-bond
    counts, then any method evaluated on that decoy set can score well without learning
    anything about binding.

    0.5 is the target: it means the decoys are indistinguishable from the actives on
    everything that should be irrelevant. DUD-E's reported failure is that this number is
    far above 0.5.
    """
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    if not active_properties or not decoy_properties:
        raise ValueError("need both actives and decoys")

    features = np.array(
        [[row[name] for name in properties] for row in [*active_properties, *decoy_properties]]
    )
    labels = np.array([1] * len(active_properties) + [0] * len(decoy_properties))
    model = make_pipeline(
        StandardScaler(), LogisticRegression(max_iter=2000, random_state=seed)
    )
    scores = cross_val_score(model, features, labels, cv=5, scoring="roc_auc")
    return float(scores.mean())

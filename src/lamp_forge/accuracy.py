"""Diagnostic accuracy metrics for LAMP assay validation reporting.

Computes sensitivity, specificity, positive/negative predictive values,
likelihood ratios, and Youden's J statistic from the 2x2 contingency table
produced by a contrived-specimen or clinical validation study.  Wilson score
confidence intervals are used throughout because LAMP validation panels are
typically small (n = 20-100 specimens), where normal-approximation CIs are
unreliable.

Intended use
------------
After wet-lab validation of a designed primer set against a contrived-
specimen panel (positive specimens at known copy numbers + negative
specimens from the target off-target panel), feed the TP/FP/TN/FN counts
into this module to produce the accuracy table required for regulatory or
quality-assurance filings:

- ISO 15189 / CLSI EP12-A3 (User Protocol for Evaluation of Qualitative
  Test Performance)
- FDA 510(k) / EUA sensitivity and specificity table format
- BioVind commercial assay acceptance criteria (Sn >= 0.95, Sp >= 0.99,
  LR+ >= 20, Youden J >= 0.90)

Usage::

    from lamp_forge.accuracy import ContingencyTable, compute_accuracy

    table = ContingencyTable(tp=47, fp=2, tn=48, fn=3)
    metrics = compute_accuracy(table, prevalence=0.05)
    print(f"Sensitivity:  {metrics.sensitivity.point_estimate:.3f}")
    print(f"Specificity:  {metrics.specificity.point_estimate:.3f}")
    print(f"PPV @5% prev: {metrics.ppv:.3f}")

References:
    Altman DG & Bland JM (1994). Diagnostic tests 1: sensitivity and
        specificity. BMJ 308:1552.  doi:10.1136/bmj.308.6943.1552
    Wilson EB (1927). Probable inference, the law of succession, and
        statistical inference. J Am Stat Assoc 22:209-212.
        doi:10.1080/01621459.1927.10502953
    CLSI EP12-A3 (2023). User Protocol for Evaluation of Qualitative Test
        Performance. Clinical and Laboratory Standards Institute, Wayne PA.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

_NORMAL = NormalDist()


def _z_score(confidence: float) -> float:
    """Two-sided z-score for a given confidence level (0 < confidence < 1)."""
    if not (0.0 < confidence < 1.0):
        raise ValueError(f"confidence must be in (0, 1), got {confidence!r}")
    return _NORMAL.inv_cdf((1.0 + confidence) / 2.0)


@dataclass(frozen=True, slots=True)
class ContingencyTable:
    """2x2 contingency table from a contrived-specimen validation study.

    All counts must be non-negative integers.  At least one positive
    reference specimen (tp + fn > 0) and one negative reference specimen
    (tn + fp > 0) must be present to compute sensitivity and specificity.

    Attributes:
        tp: True positives (target-positive specimens that tested positive).
        fp: False positives (target-negative specimens that tested positive).
        tn: True negatives (target-negative specimens that tested negative).
        fn: False negatives (target-positive specimens that tested negative).
    """

    tp: int
    fp: int
    tn: int
    fn: int

    def __post_init__(self) -> None:  # noqa: D105
        for name, val in (("tp", self.tp), ("fp", self.fp), ("tn", self.tn), ("fn", self.fn)):
            if not isinstance(val, int) or val < 0:
                raise ValueError(f"{name} must be a non-negative integer, got {val!r}")
        if self.tp + self.fn == 0:
            raise ValueError(
                "At least one positive reference specimen required (tp + fn > 0). "
                "Cannot compute sensitivity with no positive reference samples."
            )
        if self.tn + self.fp == 0:
            raise ValueError(
                "At least one negative reference specimen required (tn + fp > 0). "
                "Cannot compute specificity with no negative reference samples."
            )

    @property
    def n_positive_ref(self) -> int:
        """Number of positive-reference specimens (tp + fn)."""
        return self.tp + self.fn

    @property
    def n_negative_ref(self) -> int:
        """Number of negative-reference specimens (tn + fp)."""
        return self.tn + self.fp

    @property
    def n_total(self) -> int:
        """Total specimen count."""
        return self.tp + self.fp + self.tn + self.fn


@dataclass(frozen=True, slots=True)
class WilsonCI:
    """Wilson score confidence interval for a proportion.

    Attributes:
        point_estimate: Observed proportion k/n.
        lower: Lower CI bound (clamped to [0, 1]).
        upper: Upper CI bound (clamped to [0, 1]).
        confidence: Nominal confidence level (e.g. 0.95 for 95% CI).
    """

    point_estimate: float
    lower: float
    upper: float
    confidence: float


def wilson_ci(k: int, n: int, confidence: float = 0.95) -> WilsonCI:
    """Compute the Wilson score confidence interval for a proportion k/n.

    The Wilson interval is preferred over Wald (normal approximation) for
    small n or extreme proportions (k near 0 or n).  It never exceeds [0, 1]
    and has better coverage properties at the boundaries.

    Args:
        k: Number of successes (e.g. TP count for sensitivity).
        n: Total trials (e.g. n_positive_ref for sensitivity).
        confidence: Two-sided confidence level, e.g. 0.95 for 95% CI.

    Returns:
        WilsonCI with point_estimate, lower, upper, and confidence.

    Raises:
        ValueError: If n == 0 or confidence is outside (0, 1).
    """
    if n == 0:
        raise ValueError("n must be > 0")
    z = _z_score(confidence)
    p = k / n
    z2_over_n = z**2 / n
    denominator = 1.0 + z2_over_n
    center = (p + z2_over_n / 2.0) / denominator
    half_width = z * math.sqrt(p * (1.0 - p) / n + z2_over_n / (4.0 * n)) / denominator
    lower = min(p, max(0.0, center - half_width))
    upper = max(p, min(1.0, center + half_width))
    return WilsonCI(
        point_estimate=p,
        lower=lower,
        upper=upper,
        confidence=confidence,
    )


@dataclass(frozen=True, slots=True)
class AccuracyMetrics:
    """Complete diagnostic accuracy metrics for one prevalence scenario.

    Attributes:
        table: The source 2x2 contingency table.
        sensitivity: Wilson CI for sensitivity (TP / (TP + FN)).
        specificity: Wilson CI for specificity (TN / (TN + FP)).
        lr_positive: Positive likelihood ratio = Sn / (1 - Sp).
        lr_negative: Negative likelihood ratio = (1 - Sn) / Sp.
        ppv: Positive predictive value at the given prevalence.
        npv: Negative predictive value at the given prevalence.
        prevalence: Disease prevalence used for PPV/NPV calculation.
        youden_j: Youden's J statistic (informedness) = Sn + Sp - 1.
    """

    table: ContingencyTable
    sensitivity: WilsonCI
    specificity: WilsonCI
    lr_positive: float
    lr_negative: float
    ppv: float
    npv: float
    prevalence: float
    youden_j: float


def compute_accuracy(
    table: ContingencyTable,
    prevalence: float = 0.10,
    confidence: float = 0.95,
) -> AccuracyMetrics:
    """Compute full diagnostic accuracy metrics from a 2x2 contingency table.

    PPV and NPV are computed via Bayes' theorem at the specified disease
    prevalence (which need not match the study's reference prevalence):

    ::

        PPV = (Sn * prev) / (Sn * prev + (1-Sp) * (1-prev))
        NPV = (Sp * (1-prev)) / (Sp * (1-prev) + (1-Sn) * prev)

    Likelihood ratios are prevalence-independent and can be used to update
    pre-test odds for any given patient or sample context:

    ::

        LR+ = Sn / (1 - Sp)     [typically want >= 10 for "strong" evidence]
        LR- = (1 - Sn) / Sp     [typically want <= 0.1 for "strong" exclusion]

    Args:
        table: 2x2 contingency table (ContingencyTable).
        prevalence: Disease/target prevalence for PPV/NPV calculation (0-1).
        confidence: CI confidence level (default 0.95 for 95% CI).

    Returns:
        AccuracyMetrics with all fields populated.

    Raises:
        ValueError: If prevalence is outside (0, 1) or confidence is invalid.
    """
    if not (0.0 < prevalence < 1.0):
        raise ValueError(f"prevalence must be in (0, 1), got {prevalence!r}")

    sn_ci = wilson_ci(table.tp, table.n_positive_ref, confidence)
    sp_ci = wilson_ci(table.tn, table.n_negative_ref, confidence)

    sn = sn_ci.point_estimate
    sp = sp_ci.point_estimate

    # Likelihood ratios (point estimates only; CI for LRs requires delta method)
    lr_pos = sn / (1.0 - sp) if sp < 1.0 else math.inf
    lr_neg = (1.0 - sn) / sp if sp > 0.0 else math.nan

    # PPV and NPV via Bayes' theorem
    denom_ppv = sn * prevalence + (1.0 - sp) * (1.0 - prevalence)
    ppv = (sn * prevalence) / denom_ppv if denom_ppv > 0.0 else 0.0

    denom_npv = sp * (1.0 - prevalence) + (1.0 - sn) * prevalence
    npv = (sp * (1.0 - prevalence)) / denom_npv if denom_npv > 0.0 else 0.0

    youden_j = sn + sp - 1.0

    return AccuracyMetrics(
        table=table,
        sensitivity=sn_ci,
        specificity=sp_ci,
        lr_positive=lr_pos,
        lr_negative=lr_neg,
        ppv=ppv,
        npv=npv,
        prevalence=prevalence,
        youden_j=youden_j,
    )


def accuracy_table(
    table: ContingencyTable,
    prevalences: tuple[float, ...] = (0.01, 0.05, 0.10, 0.20),
    confidence: float = 0.95,
) -> list[AccuracyMetrics]:
    """Compute accuracy metrics across a range of prevalences.

    Sensitivity, specificity, LR+, LR-, and Youden J are prevalence-
    independent (same for all rows).  PPV and NPV vary with prevalence.

    Args:
        table: 2x2 contingency table.
        prevalences: Tuple of prevalence values to evaluate.
        confidence: CI confidence level.

    Returns:
        List of AccuracyMetrics, one per prevalence value.
    """
    return [compute_accuracy(table, p, confidence) for p in prevalences]


def write_accuracy_csv(metrics_list: list[AccuracyMetrics], path: Path) -> None:
    """Write accuracy metrics to a CSV file.

    One row per prevalence scenario.  Sensitivity, specificity, and their
    Wilson CIs are the same across rows; PPV and NPV vary with prevalence.

    Args:
        metrics_list: List of AccuracyMetrics (typically from accuracy_table).
        path: Output CSV path (parent directory must exist).
    """
    fieldnames = [
        "prevalence",
        "tp",
        "fp",
        "tn",
        "fn",
        "n_positive_ref",
        "n_negative_ref",
        "sensitivity",
        "sensitivity_lower",
        "sensitivity_upper",
        "specificity",
        "specificity_lower",
        "specificity_upper",
        "lr_positive",
        "lr_negative",
        "ppv",
        "npv",
        "youden_j",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for m in metrics_list:
            writer.writerow(
                {
                    "prevalence": round(m.prevalence, 6),
                    "tp": m.table.tp,
                    "fp": m.table.fp,
                    "tn": m.table.tn,
                    "fn": m.table.fn,
                    "n_positive_ref": m.table.n_positive_ref,
                    "n_negative_ref": m.table.n_negative_ref,
                    "sensitivity": round(m.sensitivity.point_estimate, 6),
                    "sensitivity_lower": round(m.sensitivity.lower, 6),
                    "sensitivity_upper": round(m.sensitivity.upper, 6),
                    "specificity": round(m.specificity.point_estimate, 6),
                    "specificity_lower": round(m.specificity.lower, 6),
                    "specificity_upper": round(m.specificity.upper, 6),
                    "lr_positive": round(m.lr_positive, 4)
                    if math.isfinite(m.lr_positive)
                    else "inf",
                    "lr_negative": round(m.lr_negative, 4)
                    if math.isfinite(m.lr_negative)
                    else "nan",
                    "ppv": round(m.ppv, 6),
                    "npv": round(m.npv, 6),
                    "youden_j": round(m.youden_j, 6),
                }
            )

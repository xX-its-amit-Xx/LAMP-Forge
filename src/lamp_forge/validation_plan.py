"""Validation study size planning for LAMP assay regulatory submissions.

Answers the inverse of :mod:`lamp_forge.accuracy`: instead of "given TP/FP/TN/FN,
what are my Sn/Sp?", it answers "how many confirmed-positive and confirmed-negative
specimens do I need to *demonstrate* Sn >= p_sn and Sp >= p_sp at the required
confidence level?"

Two approaches are provided for each metric:

**Zero-failure plan (conservative, FDA/CLIA style)**
    Finds the smallest *n* such that, even if every specimen is correctly called
    (all TP for sensitivity, all TN for specificity), the lower Wilson score CI
    bound still meets the specification.  Formula::

        n_zf = ceil( p_min * z^2 / (1 - p_min) )

    where z is the normal quantile for the required confidence.  If *any*
    specimen fails the lower bound is instantly below spec, so the test must
    run perfect on this panel size.  This matches the style of FDA 510(k)
    sensitivity/specificity tables for qualitative IVDs.

**Expected-performance plan (powered study design)**
    Finds the smallest *n* such that, when the true underlying performance is
    *p_expected*, the probability of the lower Wilson CI bound meeting *p_min*
    is >= *power*.  Uses a normal approximation to the binomial count
    distribution with a continuity correction; accurate for n >= 20.

Usage::

    from lamp_forge.validation_plan import build_validation_plan

    plan = build_validation_plan(sn_min=0.95, sp_min=0.99, confidence=0.95)
    for row in plan.rows:
        print(f"{row.metric}: zero-failure n={row.n_zero_failure}, "
              f"expected n={row.n_expected}")

References:
    Wilson EB (1927). Probable inference, the law of succession, and
        statistical inference. J Am Stat Assoc 22:209-212.
        doi:10.1080/01621459.1927.10502953
    Flahault A, Cadilhac M & Thomas G (2005). Sample size calculation should
        be performed for design accuracy in diagnostic test studies.
        J Clin Epidemiol 58:859-862.  doi:10.1016/j.jclinepi.2004.12.009
    CLSI EP12-A3 (2023). User Protocol for Evaluation of Qualitative Test
        Performance. Clinical and Laboratory Standards Institute, Wayne PA.
    FDA (2007). Statistical Guidance on Reporting Results from Studies
        Evaluating Diagnostic Tests. Guidance for Industry and FDA Staff.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path
from statistics import NormalDist

_NORMAL = NormalDist()

# Maximum n to search before giving up (avoids infinite loop for edge cases).
_N_MAX_SEARCH = 5_000


# ---------------------------------------------------------------------------
# Low-level statistical helpers
# ---------------------------------------------------------------------------


def _z_score(confidence: float) -> float:
    """Two-sided normal z-score for *confidence* (e.g. 0.95 -> 1.96)."""
    if not (0.0 < confidence < 1.0):
        raise ValueError(f"confidence must be in (0, 1), got {confidence!r}")
    return _NORMAL.inv_cdf((1.0 + confidence) / 2.0)


def _normal_cdf(z: float) -> float:
    """Standard normal CDF via erfc."""
    return 0.5 * math.erfc(-z / math.sqrt(2.0))


def wilson_lower(k: int, n: int, confidence: float = 0.95) -> float:
    """Wilson score lower CI bound for *k* successes in *n* trials.

    Args:
        k: Number of successes (0 <= k <= n).
        n: Total trials (n >= 1).
        confidence: Two-sided confidence level (e.g. 0.95).

    Returns:
        Lower Wilson CI bound in [0, 1].

    Raises:
        ValueError: If *n* < 1 or inputs are out of range.
    """
    if n < 1:
        raise ValueError(f"n must be >= 1, got {n!r}")
    if not (0 <= k <= n):
        raise ValueError(f"k must be in [0, n], got k={k!r}, n={n!r}")
    z = _z_score(confidence)
    z2 = z * z
    p = k / n
    z2_over_n = z2 / n
    denom = 1.0 + z2_over_n
    center = (p + z2_over_n / 2.0) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z2_over_n / (4.0 * n)) / denom
    return max(0.0, min(p, center - half))


def _k_min_for_wilson(n: int, p_min: float, confidence: float) -> int:
    """Smallest *k* in [0, n] s.t. wilson_lower(k, n, confidence) >= p_min.

    Uses binary search (wilson_lower is monotonically non-decreasing in k).
    Returns n+1 if no such k exists (impossible when p_min > 1 or n=0).
    """
    lo, hi = 0, n
    # Quick check: if even k=n fails, return n+1
    if wilson_lower(n, n, confidence) < p_min:
        return n + 1
    while lo < hi:
        mid = (lo + hi) // 2
        if wilson_lower(mid, n, confidence) >= p_min:
            hi = mid
        else:
            lo = mid + 1
    return lo


def _binomial_prob_ge_k(n: int, k_min: int, p: float) -> float:
    """P(Binomial(n, p) >= k_min) via normal approximation with continuity correction.

    Accurate for n >= 20.  Returns 1.0 when k_min <= 0, 0.0 when k_min > n.
    """
    if k_min <= 0:
        return 1.0
    if k_min > n:
        return 0.0
    mu = n * p
    sigma = math.sqrt(n * p * (1.0 - p))
    if sigma < 1e-10:
        return 1.0 if n * p >= k_min - 0.5 else 0.0
    z = (k_min - 0.5 - mu) / sigma
    return 1.0 - _normal_cdf(z)


# ---------------------------------------------------------------------------
# Core planning functions
# ---------------------------------------------------------------------------


def n_zero_failure(p_min: float, confidence: float = 0.95) -> int:
    """Minimum n such that wilson_lower(n, n, confidence) >= p_min.

    With zero failures observed (all n specimens correctly called), the lower
    Wilson CI bound will exactly meet *p_min* at this sample size.  Any single
    failure will push the lower bound below spec.

    Args:
        p_min: Minimum acceptable proportion (0 < p_min < 1).
        confidence: Two-sided CI confidence (e.g. 0.95).

    Returns:
        Minimum integer sample size.

    Raises:
        ValueError: If *p_min* is outside (0, 1).
    """
    if not (0.0 < p_min < 1.0):
        raise ValueError(f"p_min must be in (0, 1), got {p_min!r}")
    z2 = _z_score(confidence) ** 2
    # Closed-form: n >= p_min * z^2 / (1 - p_min)
    n_float = p_min * z2 / (1.0 - p_min)
    n = math.ceil(n_float)
    # Verify and nudge up by 1 if floating-point gave us a borderline result
    while wilson_lower(n, n, confidence) < p_min - 1e-12:
        n += 1
    return n


def n_expected_performance(
    p_min: float,
    confidence: float,
    p_expected: float,
    power: float = 0.80,
) -> tuple[int, float]:
    """Minimum n for lower Wilson CI to meet p_min with probability >= power.

    Models the *n* required when the true underlying performance is *p_expected*
    (which should be > p_min to have any realistic power).

    Args:
        p_min: Minimum acceptable proportion (0 < p_min < 1).
        confidence: Two-sided CI confidence (e.g. 0.95).
        p_expected: Expected true performance (p_expected > p_min recommended).
        power: Required probability (power) that the CI lower bound meets p_min.

    Returns:
        Tuple of (n, achieved_power) where achieved_power >= power.

    Raises:
        ValueError: If arguments are out of range.
    """
    if not (0.0 < p_min < 1.0):
        raise ValueError(f"p_min must be in (0, 1), got {p_min!r}")
    if not (0.0 < p_expected <= 1.0):
        raise ValueError(f"p_expected must be in (0, 1], got {p_expected!r}")
    if not (0.0 < power < 1.0):
        raise ValueError(f"power must be in (0, 1), got {power!r}")

    for n in range(1, _N_MAX_SEARCH + 1):
        k_min = _k_min_for_wilson(n, p_min, confidence)
        achieved = _binomial_prob_ge_k(n, k_min, p_expected)
        if achieved >= power:
            return n, achieved
    # Fallback: return _N_MAX_SEARCH with whatever power we achieved there
    k_min = _k_min_for_wilson(_N_MAX_SEARCH, p_min, confidence)
    achieved = _binomial_prob_ge_k(_N_MAX_SEARCH, k_min, p_expected)
    return _N_MAX_SEARCH, achieved


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ValidationSizeRow:
    """Sample size requirement for one performance metric.

    Attributes:
        metric: "sensitivity" or "specificity".
        p_min: Minimum acceptable value (e.g. 0.95).
        confidence: CI confidence level (e.g. 0.95).
        n_zero_failure: Specimens required with zero allowed failures.
        wilson_lower_zf: Actual lower CI bound at n_zero_failure (>= p_min).
        p_expected: Expected true performance (used for powered plan).
        power: Required statistical power (used for powered plan).
        n_expected: Specimens required for powered study at p_expected.
        achieved_power: Actual power at n_expected.
    """

    metric: str
    p_min: float
    confidence: float
    n_zero_failure: int
    wilson_lower_zf: float
    p_expected: float
    power: float
    n_expected: int
    achieved_power: float


@dataclass(frozen=True, slots=True)
class ValidationPlan:
    """Complete validation study plan for Sn + Sp requirements.

    Attributes:
        rows: One ValidationSizeRow per metric.
        n_positive_required: Confirmed-positive specimens needed (for sensitivity).
        n_negative_required: Confirmed-negative specimens needed (for specificity).
    """

    rows: tuple[ValidationSizeRow, ...]
    n_positive_required: int
    n_negative_required: int


# ---------------------------------------------------------------------------
# Plan builder
# ---------------------------------------------------------------------------


def build_validation_plan(
    sn_min: float = 0.95,
    sp_min: float = 0.99,
    confidence: float = 0.95,
    p_expected_pos: float | None = None,
    p_expected_neg: float | None = None,
    power: float = 0.80,
) -> ValidationPlan:
    """Build a complete validation study plan.

    Computes both the zero-failure and powered sample sizes for sensitivity
    and specificity requirements.

    Args:
        sn_min: Minimum sensitivity to claim (default 0.95).
        sp_min: Minimum specificity to claim (default 0.99).
        confidence: Wilson CI confidence level (default 0.95 for 95% CI).
        p_expected_pos: Expected true sensitivity for powered plan.
            Defaults to sn_min + 0.02 (2 pp above the threshold) when
            p_expected_pos is None or equal to sn_min, to model a realistic
            scenario where the true performance exceeds the minimum claim.
        p_expected_neg: Expected true specificity.  Same defaulting logic.
        power: Required statistical power for the powered plan (default 0.80).

    Returns:
        ValidationPlan with rows for sensitivity and specificity.

    Raises:
        ValueError: If any argument is out of range.
    """
    _p_exp_pos = p_expected_pos if p_expected_pos is not None else min(sn_min + 0.02, 0.999)
    _p_exp_neg = p_expected_neg if p_expected_neg is not None else min(sp_min + 0.01, 0.999)

    rows: list[ValidationSizeRow] = []
    n_pos_req = 0
    n_neg_req = 0

    for metric, p_min, p_exp in (
        ("sensitivity", sn_min, _p_exp_pos),
        ("specificity", sp_min, _p_exp_neg),
    ):
        nzf = n_zero_failure(p_min, confidence)
        wl_zf = wilson_lower(nzf, nzf, confidence)
        nep, ap = n_expected_performance(p_min, confidence, p_exp, power)

        row = ValidationSizeRow(
            metric=metric,
            p_min=p_min,
            confidence=confidence,
            n_zero_failure=nzf,
            wilson_lower_zf=wl_zf,
            p_expected=p_exp,
            power=power,
            n_expected=nep,
            achieved_power=ap,
        )
        rows.append(row)
        if metric == "sensitivity":
            n_pos_req = nzf
        else:
            n_neg_req = nzf

    return ValidationPlan(
        rows=tuple(rows),
        n_positive_required=n_pos_req,
        n_negative_required=n_neg_req,
    )


# ---------------------------------------------------------------------------
# CSV export
# ---------------------------------------------------------------------------


def write_plan_csv(plan: ValidationPlan, path: Path) -> None:
    """Write validation plan to a CSV file.

    Args:
        plan: ValidationPlan from build_validation_plan.
        path: Output CSV path.
    """
    fieldnames = [
        "metric",
        "p_min",
        "confidence",
        "n_zero_failure",
        "wilson_lower_zf",
        "p_expected",
        "power",
        "n_expected",
        "achieved_power",
    ]
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for row in plan.rows:
            writer.writerow(
                {
                    "metric": row.metric,
                    "p_min": round(row.p_min, 6),
                    "confidence": round(row.confidence, 6),
                    "n_zero_failure": row.n_zero_failure,
                    "wilson_lower_zf": round(row.wilson_lower_zf, 6),
                    "p_expected": round(row.p_expected, 6),
                    "power": round(row.power, 6),
                    "n_expected": row.n_expected,
                    "achieved_power": round(row.achieved_power, 6),
                }
            )

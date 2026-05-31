"""Surveillance sample-size planner for portable LAMP panel deployments.

Answers: "How many wells / animals must I test to be confident I would
detect the target pathogen if it is present at or above a minimum
prevalence?"

This is the *discovery phase* counterpart to the assay LOD tool (lod.py).
lod.py asks "how concentrated must the target be for one reaction to
detect it?"; this module asks "how many sampling units must be tested to
find at least one positive unit in the field population?"

Model (infinite population)
----------------------------
The standard "freedom-from-disease" surveillance sample-size formula:

    n = ceil[ log(1 - confidence) / log(1 - prevalence * sensitivity) ]

where:

  prevalence   - minimum design prevalence (fraction of units that are
                 positive, e.g. 0.05 for 5%).
  sensitivity  - assay analytical sensitivity (probability that the LAMP
                 reaction returns a positive result given the target is
                 truly present above the assay LOD).  Default 0.95.
  confidence   - required probability of detecting at least one positive
                 unit if the true prevalence equals the design prevalence.
                 Default 0.95 (95%).

Finite-population correction
-----------------------------
When the population size N is known and n approaches N, the
hypergeometric correction avoids over-sampling:

    n_finite = ceil[ n_infinite * N / (n_infinite + N - 1) ]

When n_infinite >= N the entire population must be tested; n_finite = N.

References:
----------
Cannon RM & Roe RT (1982). Livestock Disease Surveys: A Field Manual
for Veterinarians. Aust Bureau of Animal Health, Canberra.

Cameron AR & Baldock FC (1998). A new probability formula for surveys to
substantiate freedom from disease. *Prev Vet Med* 34:1-17.
doi:10.1016/S0167-5877(97)00081-0
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SamplingParams:
    """Parameters for a surveillance sample-size calculation.

    Attributes:
        target_prevalence: Minimum prevalence to detect, as a fraction in
            (0, 1).  E.g. 0.05 means "I want to detect the condition if
            >= 5% of units are positive."
        confidence: Required detection confidence, as a fraction in (0, 1).
            E.g. 0.95 means "95% probability of finding at least one
            positive sample if the true prevalence >= target_prevalence."
        test_sensitivity: Probability that the LAMP assay returns a
            positive result when the target pathogen is truly present above
            the assay LOD.  Must be in (0, 1].  Default 0.95.
        population_size: Total number of units in the population (e.g. wells
            in a field, animals in a herd).  ``None`` (default) means the
            population is effectively infinite; supply a positive integer to
            apply the finite-population hypergeometric correction.
    """

    target_prevalence: float
    confidence: float
    test_sensitivity: float = 0.95
    population_size: int | None = None

    def __post_init__(self) -> None:  # noqa: D105
        if not (0.0 < self.target_prevalence < 1.0):
            raise ValueError(f"target_prevalence must be in (0, 1), got {self.target_prevalence!r}")
        if not (0.0 < self.confidence < 1.0):
            raise ValueError(f"confidence must be in (0, 1), got {self.confidence!r}")
        if not (0.0 < self.test_sensitivity <= 1.0):
            raise ValueError(f"test_sensitivity must be in (0, 1], got {self.test_sensitivity!r}")
        if self.population_size is not None and self.population_size < 1:
            raise ValueError(
                f"population_size must be a positive integer, got {self.population_size!r}"
            )


@dataclass(frozen=True, slots=True)
class SamplingEstimate:
    """Sample-size estimate for one parameter combination.

    Attributes:
        target_prevalence: Design prevalence used in the calculation.
        confidence: Required detection confidence used in the calculation.
        test_sensitivity: Assay sensitivity used in the calculation.
        population_size: Population size if finite-population correction was
            applied; ``None`` for infinite-population calculation.
        n_samples: Minimum number of samples required.
        achieved_confidence: Back-calculated probability of detection with
            exactly ``n_samples`` tests (>= requested confidence due to
            ceiling rounding).
        finite_population_corrected: ``True`` when the hypergeometric
            correction reduced the sample size below the infinite-population
            estimate.
    """

    target_prevalence: float
    confidence: float
    test_sensitivity: float
    population_size: int | None
    n_samples: int
    achieved_confidence: float
    finite_population_corrected: bool


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def _infinite_n(prevalence: float, confidence: float, sensitivity: float) -> int:
    """Minimum sample size for an infinite population.

    Args:
        prevalence: Design prevalence in (0, 1).
        confidence: Required confidence in (0, 1).
        sensitivity: Assay sensitivity in (0, 1].

    Returns:
        Smallest integer n such that P(detect) >= confidence.
    """
    effective = prevalence * sensitivity
    if effective <= 0.0:
        return 0
    if effective >= 1.0:
        return 1
    n_exact = math.log(1.0 - confidence) / math.log(1.0 - effective)
    return math.ceil(n_exact)


def _achieved_confidence_infinite(n: int, prevalence: float, sensitivity: float) -> float:
    """Back-calculate P(at least one positive in n tests) for an infinite population.

    Uses the binomial complement: 1 - (1 - prevalence * sensitivity)^n.

    Args:
        n: Number of samples.
        prevalence: True prevalence in (0, 1).
        sensitivity: Assay sensitivity in (0, 1].

    Returns:
        Detection probability in [0, 1].
    """
    effective = prevalence * sensitivity
    if n <= 0 or effective <= 0.0:
        return 0.0
    if effective >= 1.0:
        return 1.0
    return 1.0 - (1.0 - effective) ** n


def _achieved_confidence_finite(n: int, N: int, prevalence: float, sensitivity: float) -> float:
    """P(at least one positive in n tests) for a finite population size N.

    Uses the hypergeometric complement probability.

    Args:
        n: Number of samples (n <= N).
        N: Population size.
        prevalence: True prevalence in (0, 1).
        sensitivity: Assay sensitivity in (0, 1].

    Returns:
        Detection probability in [0, 1].
    """
    D = round(prevalence * N)
    if D <= 0:
        return 0.0
    if n >= N:
        return 1.0

    # P(zero positives in n draws from N units with D positives)
    # = C(N-D, n) / C(N, n) via log-space computation to avoid overflow
    neg_count = N - D
    if n > neg_count:
        return 1.0  # impossible to draw n negatives when fewer than n are available

    # log C(neg_count, n) - log C(N, n)
    log_p_zero = 0.0
    for i in range(n):
        log_p_zero += math.log(neg_count - i) - math.log(N - i)

    p_zero = math.exp(log_p_zero)
    # Apply sensitivity: each truly-positive unit is detected with probability Se
    # For the conservative lower bound, compute as if sensitivity applies per unit
    # P(zero detected in n) = P(zero positive in n) + sum over k>0 of
    #   P(k positive in n) * (1-Se)^k.  For simplicity, use the
    # effective prevalence = D*Se / N approximation (equivalent to binomial
    # blend), which is conservative and widely used in practice:
    effective_prevalence = (D * sensitivity) / N
    if effective_prevalence >= 1.0:
        return 1.0
    if effective_prevalence <= 0.0:
        return 0.0
    # Re-compute with effective prevalence using hypergeometric complement
    # (approximation: treat D_eff = round(N * effective_prevalence))
    D_eff = max(1, round(N * effective_prevalence))
    neg_eff = N - D_eff
    if n > neg_eff:
        return 1.0
    log_p_zero_eff = 0.0
    for i in range(n):
        log_p_zero_eff += math.log(max(1, neg_eff - i)) - math.log(N - i)
    _ = p_zero  # original unused after effective-prevalence adjustment
    return 1.0 - math.exp(log_p_zero_eff)


def sample_size(params: SamplingParams) -> SamplingEstimate:
    """Calculate minimum sample size for surveillance detection.

    For an infinite population this uses the exact binomial complement
    inverse.  For a finite population the hypergeometric correction is
    applied, which can substantially reduce the required n when prevalence
    is high relative to population size.

    Args:
        params: Sampling-plan parameters.

    Returns:
        :class:`SamplingEstimate` with the minimum n and achieved confidence.

    Raises:
        ValueError: Propagated from :class:`SamplingParams` validation if
            any parameter is out of range.
    """
    n_inf = _infinite_n(
        params.target_prevalence,
        params.confidence,
        params.test_sensitivity,
    )

    fpc_applied = False
    N = params.population_size
    if N is not None:
        if n_inf >= N:
            n_final = N
        else:
            n_raw = n_inf * N / (n_inf + N - 1)
            n_final = math.ceil(n_raw)
            fpc_applied = n_final < n_inf
    else:
        n_final = n_inf

    if N is not None:
        achieved = _achieved_confidence_finite(
            n_final, N, params.target_prevalence, params.test_sensitivity
        )
        # The FPC approximation can slightly undershoot; bump n until confident.
        while achieved < params.confidence and n_final < N:
            n_final += 1
            achieved = _achieved_confidence_finite(
                n_final, N, params.target_prevalence, params.test_sensitivity
            )
        fpc_applied = n_final < n_inf
    else:
        achieved = _achieved_confidence_infinite(
            n_final, params.target_prevalence, params.test_sensitivity
        )

    return SamplingEstimate(
        target_prevalence=params.target_prevalence,
        confidence=params.confidence,
        test_sensitivity=params.test_sensitivity,
        population_size=N,
        n_samples=n_final,
        achieved_confidence=achieved,
        finite_population_corrected=fpc_applied,
    )


def sampling_plan(
    prevalences: tuple[float, ...],
    confidence: float = 0.95,
    sensitivity: float = 0.95,
    population_size: int | None = None,
) -> list[SamplingEstimate]:
    """Build a sampling-plan table across multiple design prevalences.

    A convenience wrapper around :func:`sample_size` that produces the
    table most useful for pre-campaign planning: how does the required
    sample size change as I lower the target prevalence?

    Args:
        prevalences: Design prevalences to evaluate (fractions in (0, 1)).
        confidence: Required detection confidence.  Default 0.95.
        sensitivity: Assay analytical sensitivity.  Default 0.95.
        population_size: Finite population size; ``None`` for infinite.

    Returns:
        List of :class:`SamplingEstimate`, one per prevalence, in input order.
    """
    return [
        sample_size(
            SamplingParams(
                target_prevalence=p,
                confidence=confidence,
                test_sensitivity=sensitivity,
                population_size=population_size,
            )
        )
        for p in prevalences
    ]


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_sampling_csv(estimates: list[SamplingEstimate], path: Path) -> None:
    """Write a sampling-plan table to a CSV file.

    Args:
        estimates: List of :class:`SamplingEstimate` to serialise.
        path: Output file path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "target_prevalence_pct",
                "confidence_pct",
                "test_sensitivity_pct",
                "population_size",
                "n_samples",
                "achieved_confidence_pct",
                "finite_population_corrected",
            ]
        )
        for e in estimates:
            writer.writerow(
                [
                    f"{e.target_prevalence * 100:.2f}",
                    f"{e.confidence * 100:.1f}",
                    f"{e.test_sensitivity * 100:.1f}",
                    "" if e.population_size is None else str(e.population_size),
                    str(e.n_samples),
                    f"{e.achieved_confidence * 100:.2f}",
                    "yes" if e.finite_population_corrected else "no",
                ]
            )

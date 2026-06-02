r"""Swine respiratory / PRDC outbreak trend analysis for on-farm surveillance.

In BioVind-style farm biosecurity, a single LAMP panel result tells you whether
a pathogen is present right now.  Trend analysis across consecutive sampling
intervals reveals whether the respiratory situation is improving, stable, or
deteriorating -- the key signal for escalating or relaxing PRDC biosecurity
and treatment decisions.

This module converts a chronological sequence of three-pathogen swine
respiratory LAMP panel results into a :class:`SwineRespTrend` that captures:

* Overall biosecurity **trajectory** (EMERGING / STABLE_CLEAR /
  STABLE_ENDEMIC / RESOLVING / INSUFFICIENT_DATA).
* **PRRSV first-detection and clearance events** -- the movement-restriction
  inflection points that require immediate response before respiratory
  disease escalates across the barn or to adjacent premises.
* Worst alert level observed across the time series.

Typical workflow::

    lamp-forge swine-respiratory-risk --prrsv \
        --out-json results/barn-A/2026-01.json
    lamp-forge swine-respiratory-risk --prrsv --pcv2 \
        --out-json results/barn-A/2026-02.json
    lamp-forge swine-respiratory-risk --pcv2 --mhp \
        --out-json results/barn-A/2026-03.json

    lamp-forge swine-respiratory-trend \
        --resp-result results/barn-A/2026-01.json \
        --resp-result results/barn-A/2026-02.json \
        --resp-result results/barn-A/2026-03.json \
        --out-json results/barn-A/resp_trend_Q1.json

CSV format (header required)::

    sample_id,date,prrsv,pcv2,mhp
    BarnA-2026-01,2026-01-15,1,0,0
    BarnA-2026-02,2026-02-15,1,1,0
    BarnA-2026-03,2026-03-15,0,1,1

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    Holtkamp DJ et al. (2013) Assessment of the economic impact of PRRSV
        on United States pork producers. J Swine Health Prod 21(2):72-84.
    Segalés J, Allan GM, Domingo M (2005) Porcine circovirus diseases.
        Anim Health Res Rev 6(2):119-142. doi:10.1079/AHR2005106
    Maes D et al. (2018) Swine enzootic pneumonia: a review. Vet Rec Open
        5(1):e000228. doi:10.1136/vetreco-2017-000228
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.swine_respiratory_risk import (
    PrdcAlertLevel,
    PrdcFlags,
    assess_prdc_risk,
)

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

# Alert level ordering for worst-level computation.
_LEVEL_ORDER: list[PrdcAlertLevel] = [
    PrdcAlertLevel.NEGATIVE,
    PrdcAlertLevel.LOW,
    PrdcAlertLevel.MODERATE,
    PrdcAlertLevel.HIGH,
    PrdcAlertLevel.CRITICAL,
]


@dataclass(frozen=True, slots=True)
class SwineRespRecord:
    """One chronological swine respiratory LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"BarnA-2026-01"``).
        flags: Three-pathogen positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the sample
            collection date.  Used only for display and export; ordering is
            based on sequence position in the input list.
    """

    sample_id: str
    flags: PrdcFlags
    date: str | None = None


class SwineRespTrendDirection(StrEnum):
    """Swine respiratory / PRDC trend direction derived from a time series.

    Attributes:
        EMERGING: PRRSV newly detected or overall pathogen burden is
            increasing.  Barn biosecurity is failing or under acute threat.
        STABLE_CLEAR: All panels consistently negative; no active infection
            detected across the monitoring period.
        STABLE_ENDEMIC: Pathogens (typically Mhp or PCV2) consistently
            detected without a clear upward or downward trend; managed endemic
            state.
        RESOLVING: PRRSV previously detected but now cleared, or overall
            pathogen burden is decreasing; biosecurity response appears
            effective.
        INSUFFICIENT_DATA: Fewer than two samples; a trend cannot be computed.
    """

    EMERGING = "EMERGING"
    STABLE_CLEAR = "STABLE_CLEAR"
    STABLE_ENDEMIC = "STABLE_ENDEMIC"
    RESOLVING = "RESOLVING"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class SwineRespTrend:
    """Swine respiratory / PRDC trajectory assessment for a time series.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        prrsv_newly_detected: True if PRRSV is positive in the most recent
            sample after being absent in all prior samples -- the
            movement-restriction inflection point requiring immediate
            veterinary assessment.
        prrsv_cleared: True if PRRSV was positive in at least one prior
            sample but is negative in the most recent sample.
        worst_alert_level: Highest :class:`~lamp_forge.swine_respiratory_risk.PrdcAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: SwineRespTrendDirection
    prrsv_newly_detected: bool
    prrsv_cleared: bool
    worst_alert_level: PrdcAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[SwineRespRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample pathogen flags in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "prrsv_newly_detected": self.prrsv_newly_detected,
            "prrsv_cleared": self.prrsv_cleared,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "prrsv": r.flags.prrsv,
                    "pcv2": r.flags.pcv2,
                    "mhp": r.flags.mhp,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _burden(flags: PrdcFlags) -> int:
    """Count the number of positive pathogens in a single panel result."""
    return int(flags.prrsv) + int(flags.pcv2) + int(flags.mhp)


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[SwineRespRecord, ...],
) -> tuple[SwineRespTrendDirection, bool, bool]:
    """Compute trend direction and PRRSV event flags.

    PRRSV first-detection supersedes all other signals and returns EMERGING
    immediately.  Otherwise the direction is determined by changes in overall
    pathogen burden across the time series.

    Args:
        records: At least two :class:`SwineRespRecord` in chronological order.

    Returns:
        Tuple of (direction, prrsv_newly_detected, prrsv_cleared).
    """
    n = len(records)
    last_flags = records[-1].flags

    prrsv_in_prior = any(r.flags.prrsv for r in records[:-1])
    prrsv_in_last = last_flags.prrsv

    prrsv_newly_detected = prrsv_in_last and not prrsv_in_prior
    prrsv_cleared = prrsv_in_prior and not prrsv_in_last

    if prrsv_newly_detected:
        return SwineRespTrendDirection.EMERGING, prrsv_newly_detected, prrsv_cleared

    burdens = [_burden(r.flags) for r in records]

    if all(b == 0 for b in burdens):
        return SwineRespTrendDirection.STABLE_CLEAR, prrsv_newly_detected, prrsv_cleared

    if burdens[-1] == 0 and burdens[0] > 0:
        return SwineRespTrendDirection.RESOLVING, prrsv_newly_detected, prrsv_cleared

    mid = max(1, n // 2)
    early_mean = sum(burdens[:mid]) / mid
    recent_mean = sum(burdens[mid:]) / (n - mid)

    if recent_mean - early_mean >= 0.5:
        return SwineRespTrendDirection.EMERGING, prrsv_newly_detected, prrsv_cleared
    if early_mean - recent_mean >= 0.5:
        return SwineRespTrendDirection.RESOLVING, prrsv_newly_detected, prrsv_cleared

    return SwineRespTrendDirection.STABLE_ENDEMIC, prrsv_newly_detected, prrsv_cleared


def _build_trend_text(
    direction: SwineRespTrendDirection,
    n_samples: int,
    prrsv_newly_detected: bool,
    prrsv_cleared: bool,
    worst_level: PrdcAlertLevel,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a swine respiratory trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        prrsv_newly_detected: True if PRRSV first appeared in the last sample.
        prrsv_cleared: True if PRRSV was present before but is now absent.
        worst_level: Highest alert level observed across all records.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is SwineRespTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the "
            "swine respiratory / PRDC trajectory.",
            "Collect the next scheduled monitoring sample and re-run "
            "swine-respiratory-trend after the second sample is available.",
        )

    if direction is SwineRespTrendDirection.EMERGING:
        if prrsv_newly_detected:
            interp = (
                "PRRSV newly detected in the most recent sample after being absent in "
                "all prior intervals. The respiratory biosecurity situation is "
                "escalating -- movement restriction and immediate veterinary assessment "
                "are required before PRRSV spreads to adjacent barns or premises. "
                "PRRSV spreads by aerosol at distances up to 3 km under favourable "
                "wind conditions and can devastate a PRRSV-naive grow-finish herd "
                "within 4-6 weeks."
            )
            action = (
                "CRITICAL: Restrict all pig movement from the affected unit to "
                "PRRSV-negative premises. Notify the site veterinarian and the farm "
                "biosecurity coordinator immediately. Assess PRRSV vaccination status "
                "of the affected cohort and the breeding herd. Consider whole-herd "
                "modified-live vaccine (MLV) exposure protocol under veterinary "
                "supervision. Increase monitoring frequency across adjacent barns. "
                "Implement air-filtration and nose-to-nose contact prevention if "
                "available. Submit additional samples for genotyping to identify the "
                "PRRSV lineage and inform the vaccination strategy."
            )
        else:
            interp = (
                f"Swine respiratory pathogen burden is increasing across the "
                f"{n_samples}-sample time series. Barn biosecurity may be deteriorating "
                f"or new pathogen introductions are occurring (worst alert: {worst_level.value})."
            )
            action = (
                "Increase monitoring frequency. Notify the site veterinarian. "
                "Review biosecurity protocols (vehicle decontamination, dedicated "
                "barn PPE, batch management). Assess vaccination compliance for PRRSV "
                "and PCV2. If clinical signs worsen (increased coughing, dyspnoea, "
                "mortality), implement barn-level movement restrictions immediately."
            )
        return interp, action

    if direction is SwineRespTrendDirection.RESOLVING:
        if prrsv_cleared:
            interp = (
                f"PRRSV previously detected but now negative in the most recent "
                f"sample across {n_samples} intervals. Biosecurity response appears "
                "effective, but confirmatory clearance testing is recommended before "
                "relaxing movement restrictions."
            )
            action = (
                "Maintain movement restrictions and enhanced biosecurity until at "
                "least three consecutive PRRSV-negative results are obtained. Confirm "
                "clearance with RT-PCR at an accredited laboratory before restoring "
                "normal pig movement. Continue monitoring PCV2 and Mhp; endemic "
                "pathogens may persist after PRRSV clears. Review the PRRSV "
                "elimination or stabilisation protocol with the herd health veterinarian."
            )
        else:
            interp = (
                f"Swine respiratory pathogen burden is decreasing across the "
                f"{n_samples}-sample time series. Biosecurity and treatment measures "
                "appear to be working."
            )
            action = (
                "Maintain the current biosecurity and treatment programme. Continue "
                "regular surveillance until at least three consecutive negative or "
                "low-burden results are achieved before relaxing enhanced controls."
            )
        return interp, action

    if direction is SwineRespTrendDirection.STABLE_CLEAR:
        interp = (
            f"All swine respiratory / PRDC targets negative across {n_samples} "
            "monitoring interval(s). No active respiratory pathogen detected."
        )
        action = (
            "Maintain the routine surveillance schedule. Confirm an internal "
            "extraction control (e.g. 16S rRNA) is positive to confirm extraction "
            "quality. Review biosecurity protocols before any high-risk events "
            "(introduction of pigs from external sources, transport, shows)."
        )
        return interp, action

    # STABLE_ENDEMIC
    if worst_level in (PrdcAlertLevel.CRITICAL, PrdcAlertLevel.HIGH):
        interp = (
            f"PRRSV or multi-target PRDC detected in multiple of the {n_samples} "
            f"monitoring intervals without a clear downward trend "
            f"(worst alert: {worst_level.value}). Barn biosecurity or treatment "
            "response is insufficient to resolve the outbreak."
        )
        action = (
            "Escalate to the herd health veterinarian for a formal PRDC management "
            "plan. Investigate persistent PRRSV circulation (consider whole-herd MLV "
            "exposure or a closed-herd PRRSV elimination protocol). Review PCV2 "
            "vaccination compliance. Assess all-in/all-out batch management and "
            "air-space separation between age groups."
        )
    else:
        interp = (
            f"Endemic respiratory pathogens (PCV2 and/or Mhp) detected in one or "
            f"more of the {n_samples} monitoring intervals without a clear trend "
            f"(worst alert: {worst_level.value}). This is a managed-endemic state "
            "typical of many commercial herds."
        )
        action = (
            "Review the metaphylaxis and vaccination programme with the site "
            "veterinarian. Ensure Mhp metaphylaxis timing aligns with the expected "
            "colonisation window. Confirm PCV2 vaccination is current in all age "
            "groups. Monitor for PRRSV emergence; a PRRSV introduction into this "
            "background would rapidly escalate the clinical situation."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_swine_resp_trend(records: Sequence[SwineRespRecord]) -> SwineRespTrend:
    """Compute a swine respiratory / PRDC trajectory assessment from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`SwineRespRecord`
            (oldest first).  Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`SwineRespTrend` with trajectory direction, key event flags,
        and a field-actionable recommendation.

    Example::

        from lamp_forge.swine_respiratory_risk import PrdcFlags
        from lamp_forge.swine_respiratory_trend import (
            SwineRespRecord, analyse_swine_resp_trend
        )

        recs = [
            SwineRespRecord("BarnA-Jan", PrdcFlags(prrsv=False, pcv2=True, mhp=True)),
            SwineRespRecord("BarnA-Feb", PrdcFlags(prrsv=True, pcv2=True, mhp=True)),
        ]
        trend = analyse_swine_resp_trend(recs)
        print(trend.direction)              # SwineRespTrendDirection.EMERGING
        print(trend.prrsv_newly_detected)   # True
    """
    if not records:
        raise ValueError("records must contain at least one SwineRespRecord")

    recs = tuple(records)
    n = len(recs)

    worst: PrdcAlertLevel = PrdcAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_prdc_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    if n < 2:
        interp, action = _build_trend_text(
            SwineRespTrendDirection.INSUFFICIENT_DATA, n, False, False, worst
        )
        return SwineRespTrend(
            n_samples=n,
            direction=SwineRespTrendDirection.INSUFFICIENT_DATA,
            prrsv_newly_detected=False,
            prrsv_cleared=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, prrsv_newly_detected, prrsv_cleared = _compute_direction(recs)
    interp, action = _build_trend_text(direction, n, prrsv_newly_detected, prrsv_cleared, worst)

    return SwineRespTrend(
        n_samples=n,
        direction=direction,
        prrsv_newly_detected=prrsv_newly_detected,
        prrsv_cleared=prrsv_cleared,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[SwineRespRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`SwineRespRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier for the interval (required).
        ``date``      -- ISO date string (optional; leave blank or omit column).
        ``prrsv``, ``pcv2``, ``mhp`` -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`SwineRespRecord` in file order.
    """
    required = {"sample_id", "prrsv", "pcv2", "mhp"}
    result: list[SwineRespRecord] = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"Empty CSV file: {path}")
        missing = required - {f.strip().lower() for f in reader.fieldnames}
        if missing:
            raise ValueError(f"CSV missing required column(s): {', '.join(sorted(missing))}")
        for i, row in enumerate(reader, start=2):
            try:
                sample_id = row["sample_id"].strip()
                if not sample_id:
                    raise ValueError("sample_id must not be empty")
                date_raw = row.get("date", "").strip() or None
                flags = PrdcFlags(
                    prrsv=_is_truthy(row["prrsv"]),
                    pcv2=_is_truthy(row["pcv2"]),
                    mhp=_is_truthy(row["mhp"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(SwineRespRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_resp_json(paths: Sequence[Path]) -> list[SwineRespRecord]:
    """Load a sequence of swine-respiratory-risk JSON outputs as :class:`SwineRespRecord`.

    Each JSON file (produced by ``lamp-forge swine-respiratory-risk --out-json``)
    supplies the ``panel_flags`` dict for one sample.  The ``sample_id``
    defaults to the file stem; an optional top-level ``"date"`` key is used
    if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key.

    Returns:
        List of :class:`SwineRespRecord` in the supplied order.
    """
    result: list[SwineRespRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        pf = data.get("panel_flags")
        if not isinstance(pf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'panel_flags' dict. "
                "Generate it with: lamp-forge swine-respiratory-risk --out-json <path>"
            )
        flags = PrdcFlags(
            prrsv=bool(pf.get("prrsv", False)),
            pcv2=bool(pf.get("pcv2", False)),
            mhp=bool(pf.get("mhp", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(SwineRespRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: SwineRespTrend, path: Path) -> None:
    """Write the swine respiratory trend assessment to a JSON file.

    Args:
        trend: :class:`SwineRespTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: SwineRespTrend, path: Path) -> None:
    """Write the chronological timeline of pathogen flags to a CSV file.

    Each row represents one monitoring interval with per-pathogen boolean
    columns and the per-sample alert level from
    :func:`~lamp_forge.swine_respiratory_risk.assess_prdc_risk`.

    Args:
        trend: :class:`SwineRespTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "prrsv", "pcv2", "mhp", "alert_level"])
        for r in trend.records:
            lvl = assess_prdc_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.prrsv),
                    int(r.flags.pcv2),
                    int(r.flags.mhp),
                    lvl,
                ]
            )

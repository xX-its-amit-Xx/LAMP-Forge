r"""Swine neonatal enteric outbreak trend analysis for on-farm surveillance.

In BioVind-style farm biosecurity, a single LAMP panel result tells you whether
a pathogen is present right now.  Trend analysis across consecutive sampling
intervals reveals whether the enteric situation is improving, stable, or
deteriorating -- the key signal for escalating or relaxing barn biosecurity
and treatment decisions.

This module converts a chronological sequence of three-pathogen swine neonatal
enteric LAMP panel results into a :class:`SwineEntericTrend` that captures:

* Overall biosecurity **trajectory** (EMERGING / STABLE_CLEAR /
  STABLE_ENDEMIC / RESOLVING / INSUFFICIENT_DATA).
* **PEDV first-detection and clearance events** -- the barn-lockdown inflection
  points that require immediate response before piglet mortality escalates.
* Worst alert level observed across the time series.

Typical workflow::

    lamp-forge swine-enteric-risk --pedv \
        --out-json results/barn-A/2026-01.json
    lamp-forge swine-enteric-risk --pedv --rota-a \
        --out-json results/barn-A/2026-02.json
    lamp-forge swine-enteric-risk --rota-a \
        --out-json results/barn-A/2026-03.json

    lamp-forge swine-enteric-trend \
        --enteric-result results/barn-A/2026-01.json \
        --enteric-result results/barn-A/2026-02.json \
        --enteric-result results/barn-A/2026-03.json \
        --out-json results/barn-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,pedv,pdcov,rota_a
    BarnA-2026-01,2026-01-15,1,0,0
    BarnA-2026-02,2026-02-15,1,0,1
    BarnA-2026-03,2026-03-15,0,0,1

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    Stevenson GW et al. (2013) Emergence of PEDV in the United States.
        J Vet Diagn Invest 25(5):649-654. doi:10.1177/1040638713501675
    Chen Q et al. (2015) Porcine deltacoronavirus (PDCoV) in the United States.
        Vet Microbiol 179(3-4):285-293. doi:10.1016/j.vetmic.2015.04.022
    He X et al. (2019) Development and evaluation of a RT-LAMP assay for
        rapid detection of PEDV. Transbound Emerg Dis 66(1):431-439.
        doi:10.1111/tbed.13044
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.swine_enteric_risk import (
    SwineEntericAlertLevel,
    SwineEntericFlags,
    assess_swine_enteric_risk,
)

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

# Alert level ordering for worst-level computation.
_LEVEL_ORDER: list[SwineEntericAlertLevel] = [
    SwineEntericAlertLevel.NEGATIVE,
    SwineEntericAlertLevel.LOW,
    SwineEntericAlertLevel.MODERATE,
    SwineEntericAlertLevel.HIGH,
    SwineEntericAlertLevel.CRITICAL,
]


@dataclass(frozen=True, slots=True)
class SwineEntericRecord:
    """One chronological swine neonatal enteric LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"BarnA-2026-01"``).
        flags: Three-pathogen positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the sample
            collection date.  Used only for display and export; ordering is
            based on sequence position in the input list.
    """

    sample_id: str
    flags: SwineEntericFlags
    date: str | None = None


class SwineEntericTrendDirection(StrEnum):
    """Swine neonatal enteric trend direction derived from a time series.

    Attributes:
        EMERGING: PEDV newly detected or overall pathogen burden is
            increasing.  Barn biosecurity is failing or under acute threat.
        STABLE_CLEAR: All panels consistently negative; no active infection
            detected across the monitoring period.
        STABLE_ENDEMIC: Pathogens (typically PDCoV or Rotavirus A) consistently
            detected without a clear upward or downward trend; managed endemic
            state.
        RESOLVING: PEDV previously detected but now cleared, or overall
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
class SwineEntericTrend:
    """Swine neonatal enteric trajectory assessment for a time series.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        pedv_newly_detected: True if PEDV is positive in the most recent
            sample after being absent in all prior samples -- the barn-lockdown
            inflection point requiring immediate cessation of pig movement.
        pedv_cleared: True if PEDV was positive in at least one prior sample
            but is negative in the most recent sample.
        worst_alert_level: Highest :class:`~lamp_forge.swine_enteric_risk.SwineEntericAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: SwineEntericTrendDirection
    pedv_newly_detected: bool
    pedv_cleared: bool
    worst_alert_level: SwineEntericAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[SwineEntericRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample pathogen flags in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "pedv_newly_detected": self.pedv_newly_detected,
            "pedv_cleared": self.pedv_cleared,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "pedv": r.flags.pedv,
                    "pdcov": r.flags.pdcov,
                    "rota_a": r.flags.rota_a,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _burden(flags: SwineEntericFlags) -> int:
    """Count the number of positive pathogens in a single panel result."""
    return int(flags.pedv) + int(flags.pdcov) + int(flags.rota_a)


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[SwineEntericRecord, ...],
) -> tuple[SwineEntericTrendDirection, bool, bool]:
    """Compute trend direction and PEDV event flags.

    PEDV first-detection supersedes all other signals and returns EMERGING
    immediately.  Otherwise the direction is determined by changes in overall
    pathogen burden across the time series.

    Args:
        records: At least two :class:`SwineEntericRecord` in chronological
            order.

    Returns:
        Tuple of (direction, pedv_newly_detected, pedv_cleared).
    """
    n = len(records)
    last_flags = records[-1].flags

    # PEDV event flags.
    pedv_in_prior = any(r.flags.pedv for r in records[:-1])
    pedv_in_last = last_flags.pedv

    pedv_newly_detected = pedv_in_last and not pedv_in_prior
    pedv_cleared = pedv_in_prior and not pedv_in_last

    # A new PEDV detection always triggers EMERGING -- it is the highest-
    # severity event in this panel and requires immediate action before any
    # other classification matters.
    if pedv_newly_detected:
        return SwineEntericTrendDirection.EMERGING, pedv_newly_detected, pedv_cleared

    burdens = [_burden(r.flags) for r in records]

    # All zero throughout -> STABLE_CLEAR.
    if all(b == 0 for b in burdens):
        return SwineEntericTrendDirection.STABLE_CLEAR, pedv_newly_detected, pedv_cleared

    # Last burden is zero and first was positive -> RESOLVING.
    if burdens[-1] == 0 and burdens[0] > 0:
        return SwineEntericTrendDirection.RESOLVING, pedv_newly_detected, pedv_cleared

    # Compare first-half vs second-half burden means.
    mid = max(1, n // 2)
    early_mean = sum(burdens[:mid]) / mid
    recent_mean = sum(burdens[mid:]) / (n - mid)

    if recent_mean - early_mean >= 0.5:
        return SwineEntericTrendDirection.EMERGING, pedv_newly_detected, pedv_cleared
    if early_mean - recent_mean >= 0.5:
        return SwineEntericTrendDirection.RESOLVING, pedv_newly_detected, pedv_cleared

    # Stable with some positives -> STABLE_ENDEMIC.
    return SwineEntericTrendDirection.STABLE_ENDEMIC, pedv_newly_detected, pedv_cleared


def _build_trend_text(
    direction: SwineEntericTrendDirection,
    n_samples: int,
    pedv_newly_detected: bool,
    pedv_cleared: bool,
    worst_level: SwineEntericAlertLevel,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a swine enteric trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        pedv_newly_detected: True if PEDV first appeared in the last sample.
        pedv_cleared: True if PEDV was present before but is now absent.
        worst_level: Highest alert level observed across all records.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is SwineEntericTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the "
            "enteric disease trajectory.",
            "Collect the next scheduled monitoring sample and re-run "
            "swine-enteric-trend after the second sample is available.",
        )

    if direction is SwineEntericTrendDirection.EMERGING:
        if pedv_newly_detected:
            interp = (
                "PEDV newly detected in the most recent sample after being absent in "
                "all prior intervals. The barn biosecurity situation is escalating -- "
                "immediate lockdown and pen quarantine are required before piglet "
                "mortality escalates. PEDV causes near-100% neonatal mortality and "
                "spreads via fomites and aerosolised faeces within hours of the first "
                "clinical signs."
            )
            action = (
                "CRITICAL: Immediately quarantine all affected pens and cease all pig "
                "movement. Enforce barn lockdown (dedicated PPE, boot dips, one-way "
                "traffic). Notify the site veterinarian and farm biosecurity coordinator. "
                "Submit intestinal content samples from acutely affected piglets to an "
                "accredited laboratory for confirmatory RT-PCR. Review gilt/sow PEDV "
                "vaccination history. Increase monitoring frequency on adjacent barns."
            )
        else:
            interp = (
                f"Enteric pathogen burden is increasing across the "
                f"{n_samples}-sample time series. Barn biosecurity may be deteriorating."
            )
            action = (
                "Increase monitoring frequency. Enforce enhanced biosecurity protocols "
                "(vehicle decontamination, dedicated barn PPE, strict one-way traffic "
                "flow). Review vaccination and all-in/all-out management. If clinical "
                "signs worsen or PEDV is suspected, initiate barn lockdown and notify "
                "the site veterinarian immediately."
            )
        if pedv_cleared:
            interp += (
                " Note: PEDV had cleared in the most recent sample, though overall "
                "pathogen burden is still rising."
            )
        return interp, action

    if direction is SwineEntericTrendDirection.RESOLVING:
        if pedv_cleared:
            interp = (
                f"PEDV previously detected but now negative in the most recent sample "
                f"across {n_samples} intervals. Biosecurity response appears effective, "
                "but confirmatory clearance testing is recommended before relaxing "
                "controls."
            )
            action = (
                "Maintain enhanced biosecurity and monitoring. Continue barn lockdown "
                "and enhanced PPE until at least three consecutive PEDV-negative results "
                "are obtained. Confirm clearance with RT-PCR at an accredited laboratory "
                "before restoring normal pig movement. Monitor for reintroduction from "
                "adjacent premises."
            )
        else:
            interp = (
                f"Enteric pathogen burden is decreasing across the {n_samples}-sample "
                "time series. Biosecurity measures appear to be working."
            )
            action = (
                "Maintain the current biosecurity programme. Continue regular "
                "surveillance until at least three consecutive negative results are "
                "achieved before relaxing enhanced controls."
            )
        return interp, action

    if direction is SwineEntericTrendDirection.STABLE_CLEAR:
        interp = (
            f"All swine enteric targets negative across {n_samples} monitoring "
            "interval(s). No active enteric infection detected."
        )
        action = (
            "Maintain the routine surveillance schedule. Ensure an internal extraction "
            "control (e.g. 16S rRNA) is positive to confirm extraction quality. "
            "Review biosecurity protocols before any high-risk events (introduction "
            "of new animals, transport from positive-status farms)."
        )
        return interp, action

    # STABLE_ENDEMIC
    if worst_level in (SwineEntericAlertLevel.CRITICAL, SwineEntericAlertLevel.HIGH):
        interp = (
            f"Coronavirus (PEDV or PDCoV) detected in multiple of the {n_samples} "
            f"monitoring intervals without a clear downward trend "
            f"(worst alert: {worst_level.value}). Barn biosecurity is insufficient."
        )
        action = (
            "Investigate the persistent coronavirus detection. Review the all-in/all-out "
            "farrowing management and barn sanitation between batches. Consult the site "
            "veterinarian regarding vaccination and colostrum management. Consider "
            "feedback exposure protocol under veterinary supervision."
        )
    else:
        interp = (
            f"Enteric pathogens detected in one or more of the {n_samples} monitoring "
            f"intervals without a clear trend (worst alert: {worst_level.value}). "
            "Rotavirus A is the most likely managed-endemic pathogen at this alert level."
        )
        action = (
            "Review piglet vaccination and colostrum management. Ensure adequate "
            "sanitation between batches. Consult the site veterinarian. Consider "
            "increasing monitoring frequency if clinical signs in neonatal litters "
            "are observed."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_swine_enteric_trend(records: Sequence[SwineEntericRecord]) -> SwineEntericTrend:
    """Compute a swine neonatal enteric trajectory assessment from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`SwineEntericRecord`
            (oldest first).  Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`SwineEntericTrend` with trajectory direction, key event
        flags, and a field-actionable recommendation.

    Example::

        from lamp_forge.swine_enteric_risk import SwineEntericFlags
        from lamp_forge.swine_enteric_trend import (
            SwineEntericRecord, analyse_swine_enteric_trend
        )

        recs = [
            SwineEntericRecord("Barn-Jan", SwineEntericFlags(
                pedv=False, pdcov=True, rota_a=False)),
            SwineEntericRecord("Barn-Feb", SwineEntericFlags(
                pedv=True, pdcov=True, rota_a=False)),
        ]
        trend = analyse_swine_enteric_trend(recs)
        print(trend.direction)          # SwineEntericTrendDirection.EMERGING
        print(trend.pedv_newly_detected)  # True
    """
    if not records:
        raise ValueError("records must contain at least one SwineEntericRecord")

    recs = tuple(records)
    n = len(recs)

    # Worst alert level across all records.
    worst: SwineEntericAlertLevel = SwineEntericAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_swine_enteric_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    if n < 2:
        interp, action = _build_trend_text(
            SwineEntericTrendDirection.INSUFFICIENT_DATA, n, False, False, worst
        )
        return SwineEntericTrend(
            n_samples=n,
            direction=SwineEntericTrendDirection.INSUFFICIENT_DATA,
            pedv_newly_detected=False,
            pedv_cleared=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, pedv_newly_detected, pedv_cleared = _compute_direction(recs)
    interp, action = _build_trend_text(direction, n, pedv_newly_detected, pedv_cleared, worst)

    return SwineEntericTrend(
        n_samples=n,
        direction=direction,
        pedv_newly_detected=pedv_newly_detected,
        pedv_cleared=pedv_cleared,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[SwineEntericRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`SwineEntericRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier for the interval (required).
        ``date``      -- ISO date string (optional; leave blank or omit column).
        ``pedv``, ``pdcov``, ``rota_a`` -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`SwineEntericRecord` in file order.
    """
    required = {"sample_id", "pedv", "pdcov", "rota_a"}
    result: list[SwineEntericRecord] = []
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
                flags = SwineEntericFlags(
                    pedv=_is_truthy(row["pedv"]),
                    pdcov=_is_truthy(row["pdcov"]),
                    rota_a=_is_truthy(row["rota_a"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(SwineEntericRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_enteric_json(paths: Sequence[Path]) -> list[SwineEntericRecord]:
    """Load a sequence of swine-enteric-risk JSON outputs as :class:`SwineEntericRecord`.

    Each JSON file (produced by ``lamp-forge swine-enteric-risk --out-json``)
    supplies the ``panel_flags`` dict for one sample.  The ``sample_id``
    defaults to the file stem; an optional top-level ``"date"`` key is used
    if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key.

    Returns:
        List of :class:`SwineEntericRecord` in the supplied order.
    """
    result: list[SwineEntericRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        pf = data.get("panel_flags")
        if not isinstance(pf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'panel_flags' dict. "
                "Generate it with: lamp-forge swine-enteric-risk --out-json <path>"
            )
        flags = SwineEntericFlags(
            pedv=bool(pf.get("pedv", False)),
            pdcov=bool(pf.get("pdcov", False)),
            rota_a=bool(pf.get("rota_a", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(SwineEntericRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: SwineEntericTrend, path: Path) -> None:
    """Write the swine enteric trend assessment to a JSON file.

    Args:
        trend: :class:`SwineEntericTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: SwineEntericTrend, path: Path) -> None:
    """Write the chronological timeline of pathogen flags to a CSV file.

    Each row represents one monitoring interval with per-pathogen boolean
    columns and the per-sample alert level from
    :func:`~lamp_forge.swine_enteric_risk.assess_swine_enteric_risk`.

    Args:
        trend: :class:`SwineEntericTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "pedv", "pdcov", "rota_a", "alert_level"])
        for r in trend.records:
            lvl = assess_swine_enteric_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.pedv),
                    int(r.flags.pdcov),
                    int(r.flags.rota_a),
                    lvl,
                ]
            )

r"""Poultry-biosecurity outbreak trend analysis for on-farm surveillance programs.

In BioVind-style poultry biosecurity, a single LAMP panel result tells you
whether a pathogen is present right now.  Trend analysis across consecutive
sampling intervals reveals whether the biosecurity situation is improving,
stable, or deteriorating -- the key signal for quarantine escalation,
depopulation decisions, and flock-movement restrictions.

This module converts a chronological sequence of four-pathogen poultry LAMP
panel results into a :class:`PoultryTrend` that captures:

* Overall biosecurity **trajectory** (EMERGING / STABLE_CLEAR /
  STABLE_ENDEMIC / RESOLVING / INSUFFICIENT_DATA).
* **WOAH-notifiable** first-detection and clearance events -- the inflection
  points that require immediate national veterinary authority notification.
* **IBV endemic** tracking -- distinguishing a managed endemic IBV situation
  (consistently positive, common in commercial poultry) from a new-onset
  respiratory outbreak.
* Worst alert level observed across the time series.

Typical workflow::

    lamp-forge poultry-risk --ibv \
        --out-json results/flock-A/2026-01.json
    lamp-forge poultry-risk --ibv --ibdv \
        --out-json results/flock-A/2026-02.json
    lamp-forge poultry-risk --aiv \
        --out-json results/flock-A/2026-03.json

    lamp-forge poultry-trend \
        --poultry-result results/flock-A/2026-01.json \
        --poultry-result results/flock-A/2026-02.json \
        --poultry-result results/flock-A/2026-03.json \
        --out-json results/flock-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,aiv,ndv,ibdv,ibv
    FlockA-2026-01,2026-01-15,0,0,0,1
    FlockA-2026-02,2026-02-15,0,0,1,1
    FlockA-2026-03,2026-03-15,1,0,0,0

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    WOAH (2024) List of OIE-notifiable diseases, infections and infestations.
    Swayne DE (2012) Impact of vaccines and vaccination on global control of
        avian influenza. Avian Dis 56(4 suppl):818-828.
        doi:10.1637/10138-091511-Review.1
    Alexander DJ (2000) Newcastle disease and other avian paramyxoviruses.
        Rev Sci Tech OIE 19(2):443-462. doi:10.20506/rst.19.2.1231
    Cavanagh D (2007) Coronavirus avian infectious bronchitis virus.
        Vet Res 38(2):281-297. doi:10.1051/vetres:2006055
    Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
        Nucleic Acids Res 28(12):e63. doi:10.1093/nar/28.12.e63
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.poultry_risk import PoultryAlertLevel, PoultryFlags, assess_poultry_risk

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

_LEVEL_ORDER: list[PoultryAlertLevel] = [
    PoultryAlertLevel.NEGATIVE,
    PoultryAlertLevel.LOW,
    PoultryAlertLevel.MODERATE,
    PoultryAlertLevel.HIGH,
    PoultryAlertLevel.CRITICAL,
]


@dataclass(frozen=True, slots=True)
class PoultryRecord:
    """One chronological poultry biosecurity LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"FlockA-2026-01"``).
        flags: Four-pathogen positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the
            sample collection date.  Used only for display and export;
            ordering is based on sequence position in the input list.
    """

    sample_id: str
    flags: PoultryFlags
    date: str | None = None


class PoultryTrendDirection(StrEnum):
    """Poultry biosecurity trend direction derived from a time series of panel results.

    Attributes:
        EMERGING: New pathogens are being detected or the overall pathogen
            burden is increasing.  Biosecurity may be failing; a new WOAH-
            notifiable detection always triggers this direction.
        STABLE_CLEAR: All panels consistently negative; clean biosecurity
            status maintained across the monitoring period.
        STABLE_ENDEMIC: Pathogens consistently detected without a clear
            upward or downward trend; typical for IBV-positive flocks in a
            managed endemic state.
        RESOLVING: Pathogen burden is decreasing or a previously detected
            pathogen cleared in the most recent interval; biosecurity
            measures appear effective.
        INSUFFICIENT_DATA: Fewer than two samples; a trend cannot be
            computed.
    """

    EMERGING = "EMERGING"
    STABLE_CLEAR = "STABLE_CLEAR"
    STABLE_ENDEMIC = "STABLE_ENDEMIC"
    RESOLVING = "RESOLVING"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class PoultryTrend:
    """Poultry biosecurity trajectory assessment for a time series of panel results.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        notifiable_newly_detected: WOAH-notifiable pathogens (AIV/NDV) first
            detected in the most recent sample after being absent in all prior
            samples -- a critical inflection point requiring immediate national
            veterinary authority notification.
        notifiable_cleared: WOAH-notifiable pathogens detected in earlier
            samples that are now negative in the most recent sample.
        ibv_endemic: True if IBV was positive in every sample in the series
            -- indicates the flock is in a managed endemic IBV state.
        worst_alert_level: Highest :class:`~lamp_forge.poultry_risk.PoultryAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: PoultryTrendDirection
    notifiable_newly_detected: tuple[str, ...]
    notifiable_cleared: tuple[str, ...]
    ibv_endemic: bool
    worst_alert_level: PoultryAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[PoultryRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample pathogen flags in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "notifiable_newly_detected": list(self.notifiable_newly_detected),
            "notifiable_cleared": list(self.notifiable_cleared),
            "ibv_endemic": self.ibv_endemic,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "aiv": r.flags.aiv,
                    "ndv": r.flags.ndv,
                    "ibdv": r.flags.ibdv,
                    "ibv": r.flags.ibv,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _burden(flags: PoultryFlags) -> int:
    """Count the number of positive pathogens in a single panel result."""
    return int(flags.aiv) + int(flags.ndv) + int(flags.ibdv) + int(flags.ibv)


def _notifiable_set(flags: PoultryFlags) -> frozenset[str]:
    """Return the set of detected WOAH-notifiable pathogen names."""
    result: set[str] = set()
    if flags.aiv:
        result.add("AIV")
    if flags.ndv:
        result.add("NDV")
    return frozenset(result)


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[PoultryRecord, ...],
) -> tuple[PoultryTrendDirection, tuple[str, ...], tuple[str, ...]]:
    """Compute trend direction and notifiable pathogen event flags.

    Args:
        records: At least two :class:`PoultryRecord` in chronological order.

    Returns:
        Tuple of (direction, notifiable_newly_detected, notifiable_cleared).
    """
    n = len(records)
    last_flags = records[-1].flags

    prior_notifiable: set[str] = set()
    for r in records[:-1]:
        prior_notifiable |= _notifiable_set(r.flags)

    last_notifiable = _notifiable_set(last_flags)

    newly_detected = tuple(sorted(last_notifiable - prior_notifiable))
    cleared = tuple(sorted(prior_notifiable - last_notifiable))

    if newly_detected:
        return PoultryTrendDirection.EMERGING, newly_detected, cleared

    burdens = [_burden(r.flags) for r in records]

    if all(b == 0 for b in burdens):
        return PoultryTrendDirection.STABLE_CLEAR, newly_detected, cleared

    if burdens[-1] == 0 and burdens[0] > 0:
        return PoultryTrendDirection.RESOLVING, newly_detected, cleared

    mid = max(1, n // 2)
    early_mean = sum(burdens[:mid]) / mid
    recent_mean = sum(burdens[mid:]) / (n - mid)

    if recent_mean - early_mean >= 0.5:
        return PoultryTrendDirection.EMERGING, newly_detected, cleared
    if early_mean - recent_mean >= 0.5:
        return PoultryTrendDirection.RESOLVING, newly_detected, cleared

    return PoultryTrendDirection.STABLE_ENDEMIC, newly_detected, cleared


def _build_trend_text(
    direction: PoultryTrendDirection,
    n_samples: int,
    newly_detected: tuple[str, ...],
    cleared: tuple[str, ...],
    ibv_endemic: bool,
    worst_level: PoultryAlertLevel,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a poultry trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        newly_detected: WOAH-notifiable pathogens first seen in the last sample.
        cleared: WOAH-notifiable pathogens present before but absent in last sample.
        ibv_endemic: True if IBV detected in every sample.
        worst_level: Highest alert level observed across all records.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is PoultryTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the biosecurity "
            "trajectory.",
            "Collect the next scheduled monitoring sample and re-run poultry-trend "
            "after the second sample is available.",
        )

    if direction is PoultryTrendDirection.EMERGING:
        if newly_detected:
            targets = ", ".join(newly_detected)
            interp = (
                f"WOAH-notifiable pathogen(s) newly detected in the most recent sample: "
                f"{targets}. Biosecurity situation is escalating -- immediate response "
                "required."
            )
            action = (
                f"Notify the national veterinary authority immediately regarding {targets} "
                "detection. Enforce flock movement restrictions and quarantine the affected "
                "house. Deploy PPE (N95, eye protection, gloves) for all personnel entering "
                "the affected house. Submit samples to a WOAH Reference Laboratory for "
                "confirmatory testing and subtyping."
            )
        else:
            interp = (
                f"Pathogen burden is increasing across the {n_samples}-sample time series. "
                "Biosecurity measures may be insufficient -- outbreak risk is rising."
            )
            action = (
                "Increase monitoring frequency. Review and enforce flock biosecurity "
                "protocols (personnel flow, equipment and vehicle decontamination, "
                "visitor restriction). Review the vaccination programme and consult the "
                "flock health veterinarian. If a WOAH-notifiable pathogen is confirmed, "
                "notify the national veterinary authority immediately."
            )
        if cleared:
            interp += (
                f" Note: {', '.join(cleared)} cleared since earlier samples, "
                "though other pathogens are now detected."
            )
        return interp, action

    if direction is PoultryTrendDirection.RESOLVING:
        if cleared:
            targets = ", ".join(cleared)
            interp = (
                f"WOAH-notifiable pathogen(s) {targets} previously detected but now negative "
                f"in the most recent sample across {n_samples} intervals. "
                "Biosecurity response appears effective, but confirmatory clearance testing "
                "is recommended before lifting restrictions."
            )
            action = (
                "Maintain enhanced biosecurity and monitoring. Obtain confirmatory negative "
                "results from a WOAH Reference Laboratory before lifting movement "
                "restrictions or restocking. Continue surveillance for at least three "
                "consecutive negative intervals before returning to routine monitoring."
            )
        else:
            interp = (
                f"Pathogen burden is decreasing across the {n_samples}-sample time series. "
                "Biosecurity measures appear to be working."
            )
            action = (
                "Maintain the current biosecurity programme. Continue regular surveillance "
                "until at least three consecutive negative results are achieved before "
                "relaxing controls."
            )
        return interp, action

    if direction is PoultryTrendDirection.STABLE_CLEAR:
        interp = (
            f"All pathogen targets negative across {n_samples} monitoring interval(s). "
            "No active infection detected."
        )
        action = (
            "Maintain the routine surveillance schedule. Ensure sample adequacy controls "
            "(e.g. 16S rRNA or avian beta-actin control) are positive to confirm extraction "
            "quality. Review biosecurity protocols before high-risk events (introduction of "
            "new birds, live-bird market attendance, or flock movements)."
        )
        return interp, action

    # STABLE_ENDEMIC
    ibv_note = (
        " IBV consistently detected: the flock appears to be in a managed endemic IBV "
        "state; consider reviewing the IBV vaccination programme and serotype match."
        if ibv_endemic
        else ""
    )
    if worst_level in (PoultryAlertLevel.CRITICAL, PoultryAlertLevel.HIGH):
        interp = (
            f"WOAH-notifiable pathogen(s) detected across multiple of the {n_samples} "
            f"monitoring intervals without a clear downward trend "
            f"(worst alert: {worst_level.value})." + ibv_note
        )
        action = (
            "Investigate the persistent notifiable-pathogen detection. Confirm results at "
            "a WOAH Reference Laboratory and maintain quarantine until cleared. Review the "
            "biosecurity programme with the flock health veterinarian and national authority."
        )
    else:
        interp = (
            f"Pathogens detected in one or more of the {n_samples} monitoring intervals "
            f"without a clear trend (worst alert: {worst_level.value})." + ibv_note
        )
        action = (
            "Review flock biosecurity controls and the vaccination programme. Consult the "
            "flock health veterinarian; consider increased monitoring frequency if the "
            "situation does not resolve within the next two sampling intervals."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_poultry_trend(records: Sequence[PoultryRecord]) -> PoultryTrend:
    """Compute a poultry biosecurity trajectory from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`PoultryRecord` (oldest first).
            Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`PoultryTrend` with trajectory direction, key event flags, and
        a field-actionable recommendation.

    Example::

        from lamp_forge.poultry_risk import PoultryFlags
        from lamp_forge.poultry_trend import PoultryRecord, analyse_poultry_trend

        recs = [
            PoultryRecord("FlockA-Jan", PoultryFlags(
                aiv=False, ndv=False, ibdv=False, ibv=True)),
            PoultryRecord("FlockA-Feb", PoultryFlags(
                aiv=True, ndv=False, ibdv=False, ibv=False)),
        ]
        trend = analyse_poultry_trend(recs)
        print(trend.direction)                 # PoultryTrendDirection.EMERGING
        print(trend.notifiable_newly_detected) # ('AIV',)
    """
    if not records:
        raise ValueError("records must contain at least one PoultryRecord")

    recs = tuple(records)
    n = len(recs)

    worst: PoultryAlertLevel = PoultryAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_poultry_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    ibv_endemic = all(r.flags.ibv for r in recs)

    if n < 2:
        interp, action = _build_trend_text(
            PoultryTrendDirection.INSUFFICIENT_DATA, n, (), (), ibv_endemic, worst
        )
        return PoultryTrend(
            n_samples=n,
            direction=PoultryTrendDirection.INSUFFICIENT_DATA,
            notifiable_newly_detected=(),
            notifiable_cleared=(),
            ibv_endemic=ibv_endemic,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, newly_detected, cleared = _compute_direction(recs)
    interp, action = _build_trend_text(direction, n, newly_detected, cleared, ibv_endemic, worst)

    return PoultryTrend(
        n_samples=n,
        direction=direction,
        notifiable_newly_detected=newly_detected,
        notifiable_cleared=cleared,
        ibv_endemic=ibv_endemic,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[PoultryRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`PoultryRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier for the interval (required).
        ``date``      -- ISO date string (optional; leave blank or omit column).
        ``aiv``, ``ndv``, ``ibdv``, ``ibv`` -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`PoultryRecord` in file order.
    """
    required = {"sample_id", "aiv", "ndv", "ibdv", "ibv"}
    result: list[PoultryRecord] = []
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
                flags = PoultryFlags(
                    aiv=_is_truthy(row["aiv"]),
                    ndv=_is_truthy(row["ndv"]),
                    ibdv=_is_truthy(row["ibdv"]),
                    ibv=_is_truthy(row["ibv"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(PoultryRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_poultry_json(paths: Sequence[Path]) -> list[PoultryRecord]:
    """Load a sequence of poultry-risk JSON outputs as :class:`PoultryRecord` list.

    Each JSON file (produced by ``lamp-forge poultry-risk --out-json``) supplies
    the ``panel_flags`` dict for one sample.  The ``sample_id`` defaults to the
    file stem; an optional top-level ``"date"`` key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key.

    Returns:
        List of :class:`PoultryRecord` in the supplied order.
    """
    result: list[PoultryRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        pf = data.get("panel_flags")
        if not isinstance(pf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'panel_flags' dict. "
                "Generate it with: lamp-forge poultry-risk --out-json <path>"
            )
        flags = PoultryFlags(
            aiv=bool(pf.get("aiv", False)),
            ndv=bool(pf.get("ndv", False)),
            ibdv=bool(pf.get("ibdv", False)),
            ibv=bool(pf.get("ibv", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(PoultryRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: PoultryTrend, path: Path) -> None:
    """Write the poultry biosecurity trend assessment to a JSON file.

    Args:
        trend: :class:`PoultryTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: PoultryTrend, path: Path) -> None:
    """Write the chronological timeline of pathogen flags to a CSV file.

    Each row represents one monitoring interval with per-pathogen boolean
    columns and the per-sample alert level from
    :func:`~lamp_forge.poultry_risk.assess_poultry_risk`.

    Args:
        trend: :class:`PoultryTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "aiv", "ndv", "ibdv", "ibv", "alert_level"])
        for r in trend.records:
            lvl = assess_poultry_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.aiv),
                    int(r.flags.ndv),
                    int(r.flags.ibdv),
                    int(r.flags.ibv),
                    lvl,
                ]
            )

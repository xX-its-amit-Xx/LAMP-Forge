r"""Farm-biosecurity outbreak trend analysis for on-farm surveillance programs.

In BioVind-style farm biosecurity, a single LAMP panel result tells you whether
a pathogen is present right now.  Trend analysis across consecutive sampling
intervals reveals whether the biosecurity situation is improving, stable, or
deteriorating -- the key signal for quarantine escalation, treatment decisions,
and livestock-movement restrictions.

This module converts a chronological sequence of five-pathogen LAMP panel
results into a :class:`FarmTrend` that captures:

* Overall biosecurity **trajectory** (EMERGING / STABLE_CLEAR /
  STABLE_ENDEMIC / RESOLVING / INSUFFICIENT_DATA).
* **Notifiable-pathogen** first-detection and clearance events -- the inflection
  points that require immediate national veterinary authority notification.
* **PRRSV endemic** tracking -- distinguishing a managed endemic PRRSV situation
  (consistently positive without escalation) from a new-onset outbreak.
* Worst alert level observed across the time series.

Typical workflow::

    lamp-forge farm-risk --prrsv \
        --out-json results/farm-A/2026-01.json
    lamp-forge farm-risk --prrsv --aiv \
        --out-json results/farm-A/2026-02.json
    lamp-forge farm-risk --aiv \
        --out-json results/farm-A/2026-03.json

    lamp-forge farm-trend \
        --farm-result results/farm-A/2026-01.json \
        --farm-result results/farm-A/2026-02.json \
        --farm-result results/farm-A/2026-03.json \
        --out-json results/farm-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,asfv,fmdv,aiv,ndv,prrsv
    FarmA-2026-01,2026-01-15,0,0,0,0,1
    FarmA-2026-02,2026-02-15,0,0,1,0,1
    FarmA-2026-03,2026-03-15,0,0,1,0,0

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    WOAH (2024) List of OIE-notifiable diseases, infections and infestations.
    Holtkamp DJ et al. (2013) Assessment of PRRSV economic impact.
        J Swine Health Prod 21:72-84.
    Grubman MJ & Baxt B (2004) Foot-and-mouth disease.
        Clin Microbiol Rev 17(2):465-493.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.farm_risk import FarmAlertLevel, FarmPanelFlags, assess_farm_risk

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

# Alert level ordering for worst-level computation.
_LEVEL_ORDER: list[FarmAlertLevel] = [
    FarmAlertLevel.NEGATIVE,
    FarmAlertLevel.LOW,
    FarmAlertLevel.MODERATE,
    FarmAlertLevel.HIGH,
    FarmAlertLevel.CRITICAL,
]


@dataclass(frozen=True, slots=True)
class FarmRecord:
    """One chronological farm-biosecurity LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"FarmA-2026-01"``).
        flags: Five-pathogen positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the
            sample collection date.  Used only for display and export;
            ordering is based on sequence position in the input list.
    """

    sample_id: str
    flags: FarmPanelFlags
    date: str | None = None


class FarmTrendDirection(StrEnum):
    """Farm-biosecurity trend direction derived from a time series of panel results.

    Attributes:
        EMERGING: New pathogens are being detected or the overall pathogen
            burden is increasing.  Biosecurity may be failing.
        STABLE_CLEAR: All panels consistently negative; clean biosecurity
            status maintained across the monitoring period.
        STABLE_ENDEMIC: Pathogens consistently detected without a clear
            upward or downward trend; typical for PRRSV-positive herds in
            a managed endemic state.
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
class FarmTrend:
    """Farm-biosecurity trajectory assessment for a time series of panel results.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        notifiable_newly_detected: WOAH-notifiable pathogens (ASFV/FMDV/AIV/NDV)
            first detected in the most recent sample after being absent in all
            prior samples -- a critical inflection point requiring immediate
            national veterinary authority notification.
        notifiable_cleared: WOAH-notifiable pathogens detected in earlier
            samples that are now negative in the most recent sample.
        prrsv_endemic: True if PRRSV was positive in every sample in the
            series -- indicates the herd is in a managed endemic PRRSV state.
        worst_alert_level: Highest :class:`~lamp_forge.farm_risk.FarmAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: FarmTrendDirection
    notifiable_newly_detected: tuple[str, ...]
    notifiable_cleared: tuple[str, ...]
    prrsv_endemic: bool
    worst_alert_level: FarmAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[FarmRecord, ...]

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
            "prrsv_endemic": self.prrsv_endemic,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "asfv": r.flags.asfv,
                    "fmdv": r.flags.fmdv,
                    "aiv": r.flags.aiv,
                    "ndv": r.flags.ndv,
                    "prrsv": r.flags.prrsv,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _burden(flags: FarmPanelFlags) -> int:
    """Count the number of positive pathogens in a single panel result."""
    return int(flags.asfv) + int(flags.fmdv) + int(flags.aiv) + int(flags.ndv) + int(flags.prrsv)


def _notifiable_set(flags: FarmPanelFlags) -> frozenset[str]:
    """Return the set of detected WOAH-notifiable pathogen names."""
    result: set[str] = set()
    if flags.asfv:
        result.add("ASFV")
    if flags.fmdv:
        result.add("FMDV")
    if flags.aiv:
        result.add("AIV")
    if flags.ndv:
        result.add("NDV")
    return frozenset(result)


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[FarmRecord, ...],
) -> tuple[FarmTrendDirection, tuple[str, ...], tuple[str, ...]]:
    """Compute trend direction and notifiable pathogen event flags.

    The primary signal is the overall pathogen burden (count of positive targets)
    and whether any WOAH-notifiable pathogen first appeared or cleared in the most
    recent sample.  A new notifiable detection always supersedes a burden-based
    classification and returns EMERGING immediately.

    Args:
        records: At least two :class:`FarmRecord` in chronological order.

    Returns:
        Tuple of (direction, notifiable_newly_detected, notifiable_cleared).
    """
    n = len(records)
    last_flags = records[-1].flags

    # Notifiable pathogens detected at any point before the last sample.
    prior_notifiable: set[str] = set()
    for r in records[:-1]:
        prior_notifiable |= _notifiable_set(r.flags)

    last_notifiable = _notifiable_set(last_flags)

    # Newly detected: present in the last sample but absent in ALL prior samples.
    newly_detected = tuple(sorted(last_notifiable - prior_notifiable))
    # Cleared: present in at least one prior sample but absent in the last sample.
    cleared = tuple(sorted(prior_notifiable - last_notifiable))

    # A new notifiable detection always triggers EMERGING.
    if newly_detected:
        return FarmTrendDirection.EMERGING, newly_detected, cleared

    # Pathogen burden per interval.
    burdens = [_burden(r.flags) for r in records]

    # All zero throughout -> STABLE_CLEAR.
    if all(b == 0 for b in burdens):
        return FarmTrendDirection.STABLE_CLEAR, newly_detected, cleared

    # Last burden is zero and first was positive -> RESOLVING.
    if burdens[-1] == 0 and burdens[0] > 0:
        return FarmTrendDirection.RESOLVING, newly_detected, cleared

    # Compare first-half vs second-half burden means.
    mid = max(1, n // 2)
    early_mean = sum(burdens[:mid]) / mid
    recent_mean = sum(burdens[mid:]) / (n - mid)

    if recent_mean - early_mean >= 0.5:
        return FarmTrendDirection.EMERGING, newly_detected, cleared
    if early_mean - recent_mean >= 0.5:
        return FarmTrendDirection.RESOLVING, newly_detected, cleared

    # Stable with some positives -> STABLE_ENDEMIC.
    return FarmTrendDirection.STABLE_ENDEMIC, newly_detected, cleared


def _build_trend_text(
    direction: FarmTrendDirection,
    n_samples: int,
    newly_detected: tuple[str, ...],
    cleared: tuple[str, ...],
    prrsv_endemic: bool,
    worst_level: FarmAlertLevel,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a farm trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        newly_detected: WOAH-notifiable pathogens first seen in the last sample.
        cleared: WOAH-notifiable pathogens present before but absent in last sample.
        prrsv_endemic: True if PRRSV detected in every sample.
        worst_level: Highest alert level observed across all records.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is FarmTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the biosecurity trajectory.",
            "Collect the next scheduled monitoring sample and re-run farm-trend "
            "after the second sample is available.",
        )

    if direction is FarmTrendDirection.EMERGING:
        if newly_detected:
            targets = ", ".join(newly_detected)
            interp = (
                f"WOAH-notifiable pathogen(s) newly detected in the most recent sample: "
                f"{targets}. Biosecurity situation is escalating -- immediate response required."
            )
            action = (
                f"Notify the national veterinary authority immediately regarding {targets} "
                "detection. Enforce movement restrictions and quarantine the affected premises. "
                "Submit samples to a WOAH Reference Laboratory for confirmatory testing. "
                "Increase monitoring frequency on adjacent premises."
            )
        else:
            interp = (
                f"Pathogen burden is increasing across the {n_samples}-sample time series. "
                "Biosecurity measures may be insufficient -- outbreak risk is rising."
            )
            action = (
                "Increase monitoring frequency. Enforce strict biosecurity protocols "
                "(personnel, equipment, and vehicle decontamination). Review the vaccination "
                "programme and consult the site veterinarian. If a notifiable pathogen is "
                "confirmed, notify the national veterinary authority."
            )
        if cleared:
            interp += (
                f" Note: {', '.join(cleared)} cleared since earlier samples, "
                "though other pathogens are now detected."
            )
        return interp, action

    if direction is FarmTrendDirection.RESOLVING:
        if cleared:
            targets = ", ".join(cleared)
            interp = (
                f"Notifiable pathogen(s) {targets} previously detected but now negative in the "
                f"most recent sample across {n_samples} intervals. "
                "Biosecurity response appears effective, but confirmatory clearance testing "
                "is recommended."
            )
            action = (
                "Maintain enhanced biosecurity and monitoring. Obtain confirmatory negative "
                "results from a reference laboratory before lifting movement restrictions. "
                "Continue regular surveillance for at least three consecutive negative intervals."
            )
        else:
            interp = (
                f"Pathogen burden is decreasing across the {n_samples}-sample time series. "
                "Biosecurity measures appear to be working."
            )
            action = (
                "Maintain the current biosecurity programme. Continue regular surveillance "
                "until at least three consecutive negative results are achieved before relaxing "
                "controls."
            )
        return interp, action

    if direction is FarmTrendDirection.STABLE_CLEAR:
        interp = (
            f"All pathogen targets negative across {n_samples} monitoring interval(s). "
            "No active infection detected."
        )
        action = (
            "Maintain the routine surveillance schedule. Ensure sample adequacy controls "
            "(e.g. 16S rRNA) are positive to confirm extraction quality. Review biosecurity "
            "protocols before any high-risk events (introduction of new animals, movements "
            "to shows or markets)."
        )
        return interp, action

    # STABLE_ENDEMIC
    prrsv_note = (
        " PRRSV consistently detected: the herd appears to be in a managed endemic state."
        if prrsv_endemic
        else ""
    )
    if worst_level in (FarmAlertLevel.CRITICAL, FarmAlertLevel.HIGH):
        interp = (
            f"Notifiable pathogen(s) detected across multiple of the {n_samples} monitoring "
            f"intervals without a clear downward trend (worst alert: {worst_level.value})."
            + prrsv_note
        )
        action = (
            "Investigate the persistent notifiable-pathogen detection. Confirm results at a "
            "reference laboratory and maintain quarantine until cleared. Review the biosecurity "
            "programme with the site veterinarian."
        )
    else:
        interp = (
            f"Pathogens detected in one or more of the {n_samples} monitoring intervals "
            f"without a clear trend (worst alert: {worst_level.value})." + prrsv_note
        )
        action = (
            "Review site biosecurity controls. Ensure the vaccination programme is current. "
            "Consult the site veterinarian; consider increased monitoring frequency if the "
            "situation does not resolve."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_farm_trend(records: Sequence[FarmRecord]) -> FarmTrend:
    """Compute a farm-biosecurity trajectory assessment from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`FarmRecord` (oldest first).
            Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`FarmTrend` with trajectory direction, key event flags, and
        a field-actionable recommendation.

    Example::

        from lamp_forge.farm_risk import FarmPanelFlags
        from lamp_forge.farm_trend import FarmRecord, analyse_farm_trend

        recs = [
            FarmRecord("FarmA-Jan", FarmPanelFlags(
                asfv=False, fmdv=False, aiv=False, ndv=False, prrsv=True)),
            FarmRecord("FarmA-Feb", FarmPanelFlags(
                asfv=False, fmdv=False, aiv=True, ndv=False, prrsv=True)),
        ]
        trend = analyse_farm_trend(recs)
        print(trend.direction)                 # FarmTrendDirection.EMERGING
        print(trend.notifiable_newly_detected) # ('AIV',)
    """
    if not records:
        raise ValueError("records must contain at least one FarmRecord")

    recs = tuple(records)
    n = len(recs)

    # Worst alert level across all records.
    worst: FarmAlertLevel = FarmAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_farm_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    # PRRSV endemic: positive in every sample.
    prrsv_endemic = all(r.flags.prrsv for r in recs)

    if n < 2:
        interp, action = _build_trend_text(
            FarmTrendDirection.INSUFFICIENT_DATA, n, (), (), prrsv_endemic, worst
        )
        return FarmTrend(
            n_samples=n,
            direction=FarmTrendDirection.INSUFFICIENT_DATA,
            notifiable_newly_detected=(),
            notifiable_cleared=(),
            prrsv_endemic=prrsv_endemic,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, newly_detected, cleared = _compute_direction(recs)
    interp, action = _build_trend_text(direction, n, newly_detected, cleared, prrsv_endemic, worst)

    return FarmTrend(
        n_samples=n,
        direction=direction,
        notifiable_newly_detected=newly_detected,
        notifiable_cleared=cleared,
        prrsv_endemic=prrsv_endemic,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[FarmRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`FarmRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier for the interval (required).
        ``date``      -- ISO date string (optional; leave blank or omit column).
        ``asfv``, ``fmdv``, ``aiv``, ``ndv``, ``prrsv`` -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`FarmRecord` in file order.
    """
    required = {"sample_id", "asfv", "fmdv", "aiv", "ndv", "prrsv"}
    result: list[FarmRecord] = []
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
                flags = FarmPanelFlags(
                    asfv=_is_truthy(row["asfv"]),
                    fmdv=_is_truthy(row["fmdv"]),
                    aiv=_is_truthy(row["aiv"]),
                    ndv=_is_truthy(row["ndv"]),
                    prrsv=_is_truthy(row["prrsv"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(FarmRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_farm_json(paths: Sequence[Path]) -> list[FarmRecord]:
    """Load a sequence of farm-risk JSON outputs as :class:`FarmRecord` list.

    Each JSON file (produced by ``lamp-forge farm-risk --out-json``) supplies
    the ``panel_flags`` dict for one sample.  The ``sample_id`` defaults to the
    file stem; an optional top-level ``"date"`` key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key.

    Returns:
        List of :class:`FarmRecord` in the supplied order.
    """
    result: list[FarmRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        pf = data.get("panel_flags")
        if not isinstance(pf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'panel_flags' dict. "
                "Generate it with: lamp-forge farm-risk --out-json <path>"
            )
        flags = FarmPanelFlags(
            asfv=bool(pf.get("asfv", False)),
            fmdv=bool(pf.get("fmdv", False)),
            aiv=bool(pf.get("aiv", False)),
            ndv=bool(pf.get("ndv", False)),
            prrsv=bool(pf.get("prrsv", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(FarmRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: FarmTrend, path: Path) -> None:
    """Write the farm-biosecurity trend assessment to a JSON file.

    Args:
        trend: :class:`FarmTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: FarmTrend, path: Path) -> None:
    """Write the chronological timeline of pathogen flags to a CSV file.

    Each row represents one monitoring interval with per-pathogen boolean
    columns and the per-sample alert level from
    :func:`~lamp_forge.farm_risk.assess_farm_risk`.

    Args:
        trend: :class:`FarmTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "asfv", "fmdv", "aiv", "ndv", "prrsv", "alert_level"])
        for r in trend.records:
            lvl = assess_farm_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.asfv),
                    int(r.flags.fmdv),
                    int(r.flags.aiv),
                    int(r.flags.ndv),
                    int(r.flags.prrsv),
                    lvl,
                ]
            )

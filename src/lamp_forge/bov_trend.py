r"""Bovine respiratory disease (BRD) LAMP panel trend analysis.

Converts a chronological sequence of five-pathogen bovine respiratory LAMP
panel results into a :class:`BovTrend` that captures the herd-level BRD
trajectory (EMERGING / STABLE_CLEAR / STABLE_ENDEMIC / RESOLVING /
INSUFFICIENT_DATA), IBR detection and clearance events, BVDV endemic status,
and the worst alert level -- the key signals for treatment escalation,
movement restriction, and regulatory notification decisions.

BRD pathogens tracked
---------------------
    BRSV   - Bovine respiratory syncytial virus (N gene)    -- primary trigger
    BCOV   - Bovine coronavirus (N gene)                    -- respiratory + enteric
    BVDV   - Bovine viral diarrhea virus (5-UTR, types 1+2) -- immunosuppressive
    IBR    - Infect. bovine rhinotracheitis / BoHV-1 (gB)   -- mandatory control EU/UK
    MHAE   - Mannheimia haemolytica (lktA leukotoxin A)     -- primary bacterial cause

Special events tracked
----------------------
* **IBR newly detected** -- BoHV-1 mandatory eradication / control programmes in
  many EU member states, the UK, and Scandinavia require rapid regulatory notification.
  First detection in a surveillance series is a critical inflection point.
* **IBR cleared** -- relevant when monitoring herds enrolled in eradication programmes
  (must be confirmed by serology before official clearance).
* **BVDV endemic** -- BVDV positive in every sample interval suggests a persistently
  infected (PI) animal is shedding continuously in the herd.  Removing PI animals
  is the primary control measure.
* **Bacterial co-infection count** -- the number of intervals in which M. haemolytica
  was co-detected with a viral trigger (the highest-mortality BRD pattern); rising
  counts indicate inadequate antimicrobial cover.

Typical workflow::

    lamp-forge bov-risk --brsv --mhae \
        --out-json results/herd-A/2026-01.json
    lamp-forge bov-risk --brsv --bcov --mhae \
        --out-json results/herd-A/2026-02.json
    lamp-forge bov-risk --bcov \
        --out-json results/herd-A/2026-03.json

    lamp-forge bov-trend \
        --bov-result results/herd-A/2026-01.json \
        --bov-result results/herd-A/2026-02.json \
        --bov-result results/herd-A/2026-03.json \
        --out-json results/herd-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,brsv,bcov,bvdv,ibr,mhae
    HerdA-2026-01,2026-01-15,1,0,0,0,1
    HerdA-2026-02,2026-02-15,1,1,0,0,1
    HerdA-2026-03,2026-03-15,0,1,0,0,0

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    Snowder GD et al. (2006) Bovine respiratory disease in feedlot cattle: economic
        impact and preventive measures. J Anim Sci 84(4):860-870.
    Grissett GP et al. (2015) Structured literature review of responses of cattle
        to viral and bacterial pathogens causing BRD. J Vet Intern Med 29:771-794.
    Ackermann M & Engels M (2006) Pro and contra IBR eradication.
        Vet Microbiol 113(3-4):293-302.
    Baker JC (1995) Clinical manifestations of bovine viral diarrhea infection.
        Vet Clin North Am Food Anim Pract 11(3):425-445.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.bov_risk import BovAlertLevel, BovPanelFlags, assess_bov_risk

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

_LEVEL_ORDER: list[BovAlertLevel] = [
    BovAlertLevel.NEGATIVE,
    BovAlertLevel.LOW,
    BovAlertLevel.MODERATE,
    BovAlertLevel.HIGH,
    BovAlertLevel.CRITICAL,
]


@dataclass(frozen=True, slots=True)
class BovRecord:
    """One chronological bovine respiratory LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"HerdA-2026-01"``).
        flags: Five-pathogen positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the
            sample collection date.  Used only for display and export;
            ordering is based on sequence position in the input list.
    """

    sample_id: str
    flags: BovPanelFlags
    date: str | None = None


class BovTrendDirection(StrEnum):
    """BRD trend direction derived from a time series of panel results.

    Attributes:
        EMERGING: Pathogen burden increasing or a new high-weight pathogen
            (IBR, BVDV) first detected -- outbreak risk rising.
        STABLE_CLEAR: All panels consistently negative.
        STABLE_ENDEMIC: Pathogens consistently detected without a clear
            trend; typical for BVDV-positive herds with a PI animal.
        RESOLVING: Pathogen burden decreasing or a key pathogen cleared in
            the most recent interval; control measures appear effective.
        INSUFFICIENT_DATA: Fewer than two samples; trend cannot be computed.
    """

    EMERGING = "EMERGING"
    STABLE_CLEAR = "STABLE_CLEAR"
    STABLE_ENDEMIC = "STABLE_ENDEMIC"
    RESOLVING = "RESOLVING"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class BovTrend:
    """BRD trajectory assessment for a time series of bovine respiratory panel results.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        ibr_newly_detected: True if IBR (BoHV-1) was absent in all prior
            samples but positive in the most recent one -- a regulatory
            notification trigger in EU/UK/Scandinavian eradication programmes.
        ibr_cleared: True if IBR was positive in at least one earlier sample
            but is negative in the most recent sample.
        bvdv_endemic: True if BVDV was positive in every sample -- suggests a
            persistently infected (PI) animal continuously shedding in the herd.
        bacterial_coinfection_count: Number of intervals in which M. haemolytica
            (MHAE) was co-detected with at least one viral target.
        worst_alert_level: Highest :class:`~lamp_forge.bov_risk.BovAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: BovTrendDirection
    ibr_newly_detected: bool
    ibr_cleared: bool
    bvdv_endemic: bool
    bacterial_coinfection_count: int
    worst_alert_level: BovAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[BovRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample pathogen flags in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "ibr_newly_detected": self.ibr_newly_detected,
            "ibr_cleared": self.ibr_cleared,
            "bvdv_endemic": self.bvdv_endemic,
            "bacterial_coinfection_count": self.bacterial_coinfection_count,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "brsv": r.flags.brsv,
                    "bcov": r.flags.bcov,
                    "bvdv": r.flags.bvdv,
                    "ibr": r.flags.ibr,
                    "mhae": r.flags.mhae,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _burden(flags: BovPanelFlags) -> int:
    """Count the number of positive targets in a single panel result."""
    return int(flags.brsv) + int(flags.bcov) + int(flags.bvdv) + int(flags.ibr) + int(flags.mhae)


def _has_bacterial_coinfection(flags: BovPanelFlags) -> bool:
    """True when M. haemolytica is co-detected with any viral target."""
    any_viral = flags.brsv or flags.bcov or flags.bvdv or flags.ibr
    return flags.mhae and any_viral


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[BovRecord, ...],
) -> tuple[BovTrendDirection, bool, bool]:
    """Compute trend direction and IBR event flags.

    IBR (BoHV-1) is treated as the notifiable-pathogen signal for BRD trend
    analysis (analogous to WOAH-notifiable detection in farm_trend).  A new
    IBR detection always returns EMERGING; IBR clearance while other burdens
    are zero returns RESOLVING.

    Args:
        records: At least two :class:`BovRecord` in chronological order.

    Returns:
        Tuple of (direction, ibr_newly_detected, ibr_cleared).
    """
    n = len(records)
    last_flags = records[-1].flags

    ibr_in_prior = any(r.flags.ibr for r in records[:-1])
    ibr_newly_detected = last_flags.ibr and not ibr_in_prior
    ibr_cleared = ibr_in_prior and not last_flags.ibr

    if ibr_newly_detected:
        return BovTrendDirection.EMERGING, ibr_newly_detected, ibr_cleared

    burdens = [_burden(r.flags) for r in records]

    if all(b == 0 for b in burdens):
        return BovTrendDirection.STABLE_CLEAR, ibr_newly_detected, ibr_cleared

    if burdens[-1] == 0 and burdens[0] > 0:
        return BovTrendDirection.RESOLVING, ibr_newly_detected, ibr_cleared

    mid = max(1, n // 2)
    early_mean = sum(burdens[:mid]) / mid
    recent_mean = sum(burdens[mid:]) / (n - mid)

    if recent_mean - early_mean >= 0.5:
        return BovTrendDirection.EMERGING, ibr_newly_detected, ibr_cleared
    if early_mean - recent_mean >= 0.5:
        return BovTrendDirection.RESOLVING, ibr_newly_detected, ibr_cleared

    return BovTrendDirection.STABLE_ENDEMIC, ibr_newly_detected, ibr_cleared


def _build_trend_text(
    direction: BovTrendDirection,
    n_samples: int,
    ibr_newly_detected: bool,
    ibr_cleared: bool,
    bvdv_endemic: bool,
    bacterial_coinfection_count: int,
    worst_level: BovAlertLevel,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a BRD trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        ibr_newly_detected: True if IBR first appeared in the last sample.
        ibr_cleared: True if IBR was positive before but is now negative.
        bvdv_endemic: True if BVDV detected in every sample.
        bacterial_coinfection_count: Number of intervals with MHAE + viral co-detection.
        worst_level: Highest alert level across all records.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is BovTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the BRD trajectory.",
            "Collect the next scheduled monitoring sample and re-run bov-trend "
            "after the second sample is available.",
        )

    bvdv_note = (
        " BVDV detected in every interval: a persistently infected (PI) animal "
        "is the likely source -- individual antigen-ELISA screening of the herd is indicated."
        if bvdv_endemic
        else ""
    )
    coinfection_note = (
        f" Bacterial co-infection (M. haemolytica + viral trigger) detected in "
        f"{bacterial_coinfection_count} interval(s) -- the highest-mortality BRD pattern."
        if bacterial_coinfection_count > 0
        else ""
    )

    if direction is BovTrendDirection.EMERGING:
        if ibr_newly_detected:
            interp = (
                "IBR (BoHV-1) newly detected in the most recent sample. "
                "IBR is subject to mandatory eradication or control programmes in the EU, "
                "United Kingdom, and Scandinavia -- regulatory notification may be required."
                + bvdv_note
                + coinfection_note
            )
            action = (
                "Notify the site veterinarian and check local regulatory requirements for IBR "
                "reporting (mandatory in many EU/UK jurisdictions). Restrict animal movement "
                "and verify vaccination status. If the mandatory control programme applies, "
                "notify the national competent authority. Confirm with serology at a reference "
                "laboratory. Increase monitoring frequency on adjacent groups."
            )
        else:
            interp = (
                f"BRD pathogen burden is increasing across the {n_samples}-sample time series. "
                "Control measures may be insufficient -- outbreak risk is rising."
                + bvdv_note
                + coinfection_note
            )
            action = (
                "Increase monitoring frequency. Review the herd vaccination programme "
                "(BRSV/BCoV/BVDV/IBR polyvalent vaccines). Reduce stress triggers "
                "(stocking density, feed changes, mixing of groups). "
                "Consult the site veterinarian for treatment protocols, particularly "
                "if M. haemolytica co-infection is detected."
            )
        if ibr_cleared:
            interp += " Note: IBR cleared since earlier samples, though overall burden is rising."
        return interp, action

    if direction is BovTrendDirection.RESOLVING:
        if ibr_cleared:
            interp = (
                "IBR (BoHV-1) was detected in earlier samples but is now negative in the "
                f"most recent sample across {n_samples} intervals. "
                "Control measures appear effective, but sero-conversion may persist."
                + bvdv_note
                + coinfection_note
            )
            action = (
                "Maintain enhanced biosecurity and movement restrictions. Confirm IBR clearance "
                "with individual serology (gE-ELISA DIVA test) before formally lifting control "
                "programme obligations. Continue surveillance for at least three consecutive "
                "negative intervals. Notify the competent authority of the clearance if "
                "the mandatory control programme requires it."
            )
        else:
            interp = (
                f"BRD pathogen burden is decreasing across the {n_samples}-sample time series. "
                "Control measures appear to be working." + bvdv_note + coinfection_note
            )
            action = (
                "Maintain the current vaccination, biosecurity, and treatment programme. "
                "Continue surveillance until at least three consecutive negative results "
                "are achieved before relaxing controls."
            )
        return interp, action

    if direction is BovTrendDirection.STABLE_CLEAR:
        interp = (
            f"All BRD targets negative across {n_samples} monitoring interval(s). "
            "No active BRD detected."
        )
        action = (
            "Maintain the routine surveillance schedule. Ensure sample-adequacy controls "
            "(e.g. 16S rRNA) are positive to confirm extraction quality. "
            "Review the vaccination programme before high-risk events "
            "(weaning, shipping, introduction of new animals)."
        )
        return interp, action

    # STABLE_ENDEMIC
    if worst_level in (BovAlertLevel.CRITICAL, BovAlertLevel.HIGH):
        interp = (
            f"High-alert BRD pathogens detected across multiple of the {n_samples} monitoring "
            f"intervals without a clear downward trend (worst alert: {worst_level.value})."
            + bvdv_note
            + coinfection_note
        )
        action = (
            "Reassess the BRD control programme with the site veterinarian. "
            "Review antimicrobial therapy choices and susceptibility data for M. haemolytica. "
            "Ensure the polyvalent BRD vaccination schedule is current. "
            "If BVDV is consistently detected, initiate PI-animal testing "
            "(individual antigen ELISA on the full cohort)."
        )
    else:
        interp = (
            f"BRD pathogens detected in one or more of the {n_samples} monitoring intervals "
            f"without a clear trend (worst alert: {worst_level.value})."
            + bvdv_note
            + coinfection_note
        )
        action = (
            "Maintain current biosecurity and vaccination protocols. Consult the site "
            "veterinarian; consider increased monitoring frequency if clinical signs "
            "of respiratory disease are observed. Review biosecurity at re-stocking events."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_bov_trend(records: Sequence[BovRecord]) -> BovTrend:
    """Compute a BRD trajectory assessment from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`BovRecord` (oldest first).
            Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`BovTrend` with trajectory direction, IBR event flags, BVDV
        endemic status, bacterial co-infection count, and a field-actionable
        recommendation.

    Example::

        from lamp_forge.bov_risk import BovPanelFlags
        from lamp_forge.bov_trend import BovRecord, analyse_bov_trend

        recs = [
            BovRecord("HerdA-Jan", BovPanelFlags(
                brsv=True, bcov=False, bvdv=False, ibr=False, mhae=True)),
            BovRecord("HerdA-Feb", BovPanelFlags(
                brsv=True, bcov=False, bvdv=False, ibr=True, mhae=True)),
        ]
        trend = analyse_bov_trend(recs)
        print(trend.direction)          # BovTrendDirection.EMERGING
        print(trend.ibr_newly_detected) # True
    """
    if not records:
        raise ValueError("records must contain at least one BovRecord")

    recs = tuple(records)
    n = len(recs)

    worst: BovAlertLevel = BovAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_bov_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    bvdv_endemic = all(r.flags.bvdv for r in recs)
    bacterial_coinfection_count = sum(1 for r in recs if _has_bacterial_coinfection(r.flags))

    if n < 2:
        interp, action = _build_trend_text(
            BovTrendDirection.INSUFFICIENT_DATA,
            n,
            False,
            False,
            bvdv_endemic,
            bacterial_coinfection_count,
            worst,
        )
        return BovTrend(
            n_samples=n,
            direction=BovTrendDirection.INSUFFICIENT_DATA,
            ibr_newly_detected=False,
            ibr_cleared=False,
            bvdv_endemic=bvdv_endemic,
            bacterial_coinfection_count=bacterial_coinfection_count,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, ibr_newly_detected, ibr_cleared = _compute_direction(recs)
    interp, action = _build_trend_text(
        direction,
        n,
        ibr_newly_detected,
        ibr_cleared,
        bvdv_endemic,
        bacterial_coinfection_count,
        worst,
    )

    return BovTrend(
        n_samples=n,
        direction=direction,
        ibr_newly_detected=ibr_newly_detected,
        ibr_cleared=ibr_cleared,
        bvdv_endemic=bvdv_endemic,
        bacterial_coinfection_count=bacterial_coinfection_count,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[BovRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`BovRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier for the interval (required).
        ``date``      -- ISO date string (optional; leave blank or omit column).
        ``brsv``, ``bcov``, ``bvdv``, ``ibr``, ``mhae`` -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`BovRecord` in file order.
    """
    required = {"sample_id", "brsv", "bcov", "bvdv", "ibr", "mhae"}
    result: list[BovRecord] = []
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
                flags = BovPanelFlags(
                    brsv=_is_truthy(row["brsv"]),
                    bcov=_is_truthy(row["bcov"]),
                    bvdv=_is_truthy(row["bvdv"]),
                    ibr=_is_truthy(row["ibr"]),
                    mhae=_is_truthy(row["mhae"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(BovRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_bov_json(paths: Sequence[Path]) -> list[BovRecord]:
    """Load a sequence of bov-risk JSON outputs as :class:`BovRecord` list.

    Each JSON file (produced by ``lamp-forge bov-risk --out-json``) supplies
    the ``panel_flags`` dict for one sample.  The ``sample_id`` defaults to the
    file stem; an optional top-level ``"date"`` key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key.

    Returns:
        List of :class:`BovRecord` in the supplied order.
    """
    result: list[BovRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        pf = data.get("panel_flags")
        if not isinstance(pf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'panel_flags' dict. "
                "Generate it with: lamp-forge bov-risk --out-json <path>"
            )
        flags = BovPanelFlags(
            brsv=bool(pf.get("brsv", False)),
            bcov=bool(pf.get("bcov", False)),
            bvdv=bool(pf.get("bvdv", False)),
            ibr=bool(pf.get("ibr", False)),
            mhae=bool(pf.get("mhae", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(BovRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: BovTrend, path: Path) -> None:
    """Write the BRD trend assessment to a JSON file.

    Args:
        trend: :class:`BovTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: BovTrend, path: Path) -> None:
    """Write the chronological timeline of pathogen flags to a CSV file.

    Each row represents one monitoring interval with per-pathogen boolean
    columns and the per-sample alert level from
    :func:`~lamp_forge.bov_risk.assess_bov_risk`.

    Args:
        trend: :class:`BovTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "brsv", "bcov", "bvdv", "ibr", "mhae", "alert_level"])
        for r in trend.records:
            lvl = assess_bov_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.brsv),
                    int(r.flags.bcov),
                    int(r.flags.bvdv),
                    int(r.flags.ibr),
                    int(r.flags.mhae),
                    lvl,
                ]
            )

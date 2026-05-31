r"""Human point-of-care (POC) LAMP panel outbreak trend analysis.

Converts a chronological sequence of five-pathogen human POC LAMP panel
results into a :class:`PocTrend` that captures the clinic-level surveillance
trajectory (EMERGING / RESOLVING / STABLE_CLEAR / STABLE_ENDEMIC /
INSUFFICIENT_DATA), TB detection and clearance events, and the respiratory
wave activity level -- the key signals for clinic capacity planning,
public-health notification, and patient-flow management.

This module mirrors :mod:`lamp_forge.farm_trend` and
:mod:`lamp_forge.souring_trend` for the human POC vertical: a rural clinic,
urgent-care centre, or telemedicine hub accumulates POC results over time
and needs a single-command assessment of the evolving outbreak landscape.

Pathogens tracked
-----------------
    MTB    - Mycobacterium tuberculosis (rpoB) - mandatory notifiable
    CDIFF  - Clostridioides difficile (tcdB)   - HAI cluster alert
    SARS2  - SARS-CoV-2 (N gene)               - respiratory wave marker
    FLU_A  - Influenza A (M gene, pan-IAV)     - respiratory wave marker
    GAS    - Group A Streptococcus (speB)       - community pharyngitis

Typical workflow::

    lamp-forge poc-risk --mtb --out-json results/clinic-A/2026-01.json
    lamp-forge poc-risk --flu-a --sars2 \
        --out-json results/clinic-A/2026-02.json
    lamp-forge poc-risk --flu-a --out-json results/clinic-A/2026-03.json

    lamp-forge poc-trend \
        --poc-result results/clinic-A/2026-01.json \
        --poc-result results/clinic-A/2026-02.json \
        --poc-result results/clinic-A/2026-03.json \
        --out-json results/clinic-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,mtb,cdiff,sars2,flu_a,gas
    Clinic-2026-01,2026-01-15,0,0,0,0,0
    Clinic-2026-02,2026-02-15,0,0,1,1,0
    Clinic-2026-03,2026-03-15,0,0,0,1,0

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    WHO (2024) Global tuberculosis report. WHO/HTM/TB/2024.
    Shulman ST et al. (2012) IDSA GAS pharyngitis guidelines.
        Clin Infect Dis 55(10):e86-e102. doi:10.1093/cid/cis629
    McDonald LC et al. (2018) CDI clinical practice guidelines, 2017 update.
        Clin Infect Dis 66(7):e1-e48. doi:10.1093/cid/cix1085
    Notomi T et al. (2000) Loop-mediated isothermal amplification of DNA.
        Nucleic Acids Res 28(12):e63. doi:10.1093/nar/28.12.e63
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from lamp_forge.poc_risk import PocAlertLevel, PocPanelFlags, assess_poc_risk

# ---------------------------------------------------------------------------
# Public enumerations
# ---------------------------------------------------------------------------


class PocTrendDirection(StrEnum):
    """Clinic-level surveillance trajectory for the human POC panel.

    EMERGING: High-urgency pathogen newly detected, or overall pathogen
        burden rising across the time series.  Immediate action required.
    RESOLVING: Pathogen(s) present in early samples, absent in the most-
        recent sample.  Outbreak appears to be waning.
    STABLE_CLEAR: All samples consistently negative.  Clinic is not in an
        active outbreak period.
    STABLE_ENDEMIC: Consistent low-level pathogen positivity without
        escalation -- typical of an active respiratory virus season or
        endemic community GAS circulation.
    INSUFFICIENT_DATA: Fewer than two samples; no trend can be inferred.
    """

    EMERGING = "EMERGING"
    RESOLVING = "RESOLVING"
    STABLE_CLEAR = "STABLE_CLEAR"
    STABLE_ENDEMIC = "STABLE_ENDEMIC"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PocRecord:
    """One time-point in a clinic POC monitoring series.

    Attributes:
        sample_id: Identifier for this sampling event (e.g. date string,
            batch number, or clinic pseudonym for aggregate results).
        flags: Five-target positivity flags from the POC panel.
        date: Optional ISO-8601 date string (YYYY-MM-DD).
    """

    sample_id: str
    flags: PocPanelFlags
    date: str | None = None


@dataclass(frozen=True, slots=True)
class PocTrend:
    """Clinic-level surveillance trend derived from a series of POC results.

    Attributes:
        n_samples: Number of records in the time series.
        direction: Overall surveillance trajectory.
        records: Ordered tuple of records used to compute the trend.
        tb_newly_detected: True when MTB was absent in all preceding samples
            and positive in the most-recent sample -- a TB onset event
            requiring immediate public-health notification.
        tb_clearing: True when MTB was positive in the first sample and is
            negative in the most-recent sample, suggesting treatment response.
        respiratory_wave_active: True when SARS-CoV-2 and/or Influenza A
            were detected in the majority of the most-recent half of samples.
        worst_alert_level: Highest :class:`~lamp_forge.poc_risk.PocAlertLevel`
            observed across the entire time series.
        interpretation: Human-readable summary of the trend pattern.
        recommended_action: Actionable clinic guidance for the trajectory.
    """

    n_samples: int
    direction: PocTrendDirection
    records: tuple[PocRecord, ...]
    tb_newly_detected: bool
    tb_clearing: bool
    respiratory_wave_active: bool
    worst_alert_level: PocAlertLevel
    interpretation: str
    recommended_action: str

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a JSON-safe plain dict.

        Returns:
            Dict with all trend fields; per-sample timeline is a list of
            dicts each containing sample_id, date, five flag booleans, and
            the computed alert_level for that sample.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "tb_newly_detected": self.tb_newly_detected,
            "tb_clearing": self.tb_clearing,
            "respiratory_wave_active": self.respiratory_wave_active,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "mtb": r.flags.mtb,
                    "cdiff": r.flags.cdiff,
                    "sars2": r.flags.sars2,
                    "flu_a": r.flags.flu_a,
                    "gas": r.flags.gas,
                    "alert_level": assess_poc_risk(r.flags).alert_level.value,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_LEVEL_ORDER: dict[PocAlertLevel, int] = {
    PocAlertLevel.NEGATIVE: 0,
    PocAlertLevel.LOW: 1,
    PocAlertLevel.MODERATE: 2,
    PocAlertLevel.HIGH: 3,
    PocAlertLevel.CRITICAL: 4,
}


def _parse_bool(value: str) -> bool:
    return value.strip().lower() in ("1", "true", "yes")


def _worst_level(records: list[PocRecord]) -> PocAlertLevel:
    best = PocAlertLevel.NEGATIVE
    for r in records:
        lvl = assess_poc_risk(r.flags).alert_level
        if _LEVEL_ORDER[lvl] > _LEVEL_ORDER[best]:
            best = lvl
    return best


def _burden(flags: PocPanelFlags) -> int:
    return sum([flags.mtb, flags.cdiff, flags.sars2, flags.flu_a, flags.gas])


def _any_positive(flags: PocPanelFlags) -> bool:
    return any([flags.mtb, flags.cdiff, flags.sars2, flags.flu_a, flags.gas])


def _build_trend_text(
    direction: PocTrendDirection,
    tb_newly_detected: bool,
    tb_clearing: bool,
    respiratory_wave_active: bool,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for the trend.

    Args:
        direction: Computed surveillance trajectory.
        tb_newly_detected: MTB first seen in the last sample.
        tb_clearing: MTB was positive, now negative.
        respiratory_wave_active: Respiratory virus in majority of recent half.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    if tb_newly_detected:
        return (
            "Mycobacterium tuberculosis detected in the most-recent POC result "
            "with no TB signal in prior samples. This is a new TB onset event "
            "requiring immediate mandatory public-health notification. Confirmatory "
            "diagnostics (sputum smear, culture, Xpert MTB/RIF) must be initiated "
            "before the patient leaves the clinic.",
            "IMMEDIATE: Place patient in airborne isolation (negative-pressure room "
            "or equivalent; N95/FFP2 respirators for HCW). Notify the national TB "
            "programme / public health authority per local mandatory notification "
            "protocol. Initiate contact tracing for clinic exposures. Do not delay "
            "confirmatory diagnostics.",
        )

    if direction is PocTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data: at least two samples are required to determine "
            "a surveillance trend.",
            "Continue sampling and rerun trend analysis once a second result is available.",
        )

    if direction is PocTrendDirection.STABLE_CLEAR:
        return (
            "All POC panel samples are negative. No active pathogen circulation "
            "detected at this clinic.",
            "Maintain routine surveillance schedule. Confirm each negative run "
            "with a sample-adequacy control (16S LAMP). No outbreak response required.",
        )

    if direction is PocTrendDirection.RESOLVING:
        extra = (
            " MTB was detected previously and is now negative, suggesting "
            "treatment response; continue confirmatory follow-up per local "
            "TB programme end-of-treatment protocol."
            if tb_clearing
            else ""
        )
        return (
            "Pathogen positivity is declining across the sampling series. "
            "The clinic outbreak situation appears to be resolving." + extra,
            "Maintain current infection-control measures and continue surveillance "
            "until two consecutive fully-negative results are achieved. "
            "If TB clearance is involved, follow the TB programme protocol for "
            "end-of-treatment confirmation.",
        )

    if direction is PocTrendDirection.STABLE_ENDEMIC:
        resp_note = (
            " An active respiratory virus wave (SARS-CoV-2 and/or Influenza A) "
            "is the dominant driver, consistent with a seasonal outbreak."
            if respiratory_wave_active
            else ""
        )
        return (
            "Consistent low-level pathogen positivity without clear escalation. "
            "The pattern suggests stable endemic circulation or an ongoing "
            "seasonal outbreak without worsening." + resp_note,
            "Maintain heightened infection-control precautions (cohorting of "
            "symptomatic patients, droplet/contact precautions, staff PPE). "
            "Review antiviral prophylaxis eligibility for high-risk patients. "
            "Increase sampling frequency if positivity rate rises.",
        )

    # EMERGING
    resp_note = " Respiratory wave activity is elevated." if respiratory_wave_active else ""
    return (
        "Pathogen burden is increasing across the sampling series. "
        "A new or escalating outbreak is indicated." + resp_note,
        "Activate clinic outbreak response protocol. Implement cohort isolation "
        "for all pathogen-positive patients. Notify the local public health "
        "authority. Increase testing frequency and consider contact tracing "
        "for notifiable detections.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_poc_trend(records: list[PocRecord]) -> PocTrend:
    """Compute a clinic surveillance trend from an ordered series of POC records.

    Uses a split-half burden comparison to classify the trajectory.  The
    TB-onset override takes precedence over all other signals.

    Args:
        records: Chronologically ordered list of :class:`PocRecord` objects.
            Must contain at least one record.

    Returns:
        :class:`PocTrend` with direction, special-event flags, worst alert
        level, and human-readable guidance.

    Raises:
        ValueError: If ``records`` is empty.

    Example::

        from lamp_forge.poc_risk import PocPanelFlags
        from lamp_forge.poc_trend import PocRecord, analyse_poc_trend

        recs = [
            PocRecord("Jan", PocPanelFlags(mtb=False, cdiff=False,
                                           sars2=True, flu_a=True, gas=False)),
            PocRecord("Feb", PocPanelFlags(mtb=False, cdiff=False,
                                           sars2=False, flu_a=False, gas=False)),
        ]
        trend = analyse_poc_trend(recs)
        print(trend.direction)   # PocTrendDirection.RESOLVING
    """
    if not records:
        raise ValueError("analyse_poc_trend requires at least one record.")

    n = len(records)
    worst = _worst_level(records)

    if n < 2:
        return PocTrend(
            n_samples=n,
            direction=PocTrendDirection.INSUFFICIENT_DATA,
            records=tuple(records),
            tb_newly_detected=False,
            tb_clearing=False,
            respiratory_wave_active=False,
            worst_alert_level=worst,
            interpretation=(
                "Insufficient data: at least two samples are required to determine "
                "a surveillance trend."
            ),
            recommended_action=(
                "Continue sampling and rerun trend analysis once a second result is available."
            ),
        )

    first, last = records[0], records[-1]

    # TB-specific event flags
    previous_any_mtb = any(r.flags.mtb for r in records[:-1])
    tb_newly_detected = last.flags.mtb and not previous_any_mtb
    tb_clearing = first.flags.mtb and not last.flags.mtb

    # Respiratory wave: SARS2 or FLU_A in majority of recent half
    mid = max(n // 2, 1)
    recent_half = records[mid:]
    resp_positives = sum(1 for r in recent_half if r.flags.sars2 or r.flags.flu_a)
    respiratory_wave_active = resp_positives > len(recent_half) / 2

    # All-negative shortcut
    if all(not _any_positive(r.flags) for r in records):
        interp, action = _build_trend_text(PocTrendDirection.STABLE_CLEAR, False, False, False)
        return PocTrend(
            n_samples=n,
            direction=PocTrendDirection.STABLE_CLEAR,
            records=tuple(records),
            tb_newly_detected=False,
            tb_clearing=False,
            respiratory_wave_active=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
        )

    # TB onset overrides all other direction signals
    if tb_newly_detected:
        interp, action = _build_trend_text(
            PocTrendDirection.EMERGING, True, False, respiratory_wave_active
        )
        return PocTrend(
            n_samples=n,
            direction=PocTrendDirection.EMERGING,
            records=tuple(records),
            tb_newly_detected=True,
            tb_clearing=False,
            respiratory_wave_active=respiratory_wave_active,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
        )

    # Split-half burden comparison
    early_half = records[:mid]
    early_burden = sum(_burden(r.flags) for r in early_half) / len(early_half)
    recent_burden = sum(_burden(r.flags) for r in recent_half) / len(recent_half)
    delta = recent_burden - early_burden

    if delta >= 0.5:
        direction = PocTrendDirection.EMERGING
    elif (not _any_positive(last.flags) and early_burden > 0) or delta <= -0.5:
        direction = PocTrendDirection.RESOLVING
    else:
        direction = PocTrendDirection.STABLE_ENDEMIC

    interp, action = _build_trend_text(
        direction, tb_newly_detected, tb_clearing, respiratory_wave_active
    )
    return PocTrend(
        n_samples=n,
        direction=direction,
        records=tuple(records),
        tb_newly_detected=tb_newly_detected,
        tb_clearing=tb_clearing,
        respiratory_wave_active=respiratory_wave_active,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
    )


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

_REQUIRED_CSV_COLUMNS = frozenset({"sample_id", "mtb", "cdiff", "sars2", "flu_a", "gas"})


def records_from_csv(path: Path) -> list[PocRecord]:
    """Parse a monitoring CSV into an ordered list of :class:`PocRecord`.

    Required header columns: ``sample_id``, ``mtb``, ``cdiff``, ``sars2``,
    ``flu_a``, ``gas``.  An optional ``date`` column (ISO-8601) is preserved.
    Boolean columns accept ``0``/``1`` or ``true``/``false`` (case-insensitive).

    Args:
        path: Path to the CSV file.

    Returns:
        List of :class:`PocRecord` in file order.

    Raises:
        ValueError: If a required column is absent from the header.
    """
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            return []
        missing = _REQUIRED_CSV_COLUMNS - set(reader.fieldnames)
        if missing:
            raise ValueError(f"CSV is missing required columns: {', '.join(sorted(missing))}")
        records: list[PocRecord] = []
        for row in reader:
            flags = PocPanelFlags(
                mtb=_parse_bool(row["mtb"]),
                cdiff=_parse_bool(row["cdiff"]),
                sars2=_parse_bool(row["sars2"]),
                flu_a=_parse_bool(row["flu_a"]),
                gas=_parse_bool(row["gas"]),
            )
            date = row.get("date") or None
            records.append(PocRecord(sample_id=row["sample_id"], flags=flags, date=date))
    return records


def records_from_poc_json(paths: list[Path]) -> list[PocRecord]:
    """Parse ``lamp-forge poc-risk --out-json`` files into :class:`PocRecord` objects.

    Each file must contain a ``panel_flags`` dict with keys ``mtb``, ``cdiff``,
    ``sars2``, ``flu_a``, ``gas``.  The file stem is used as the sample ID;
    an optional top-level ``date`` string is preserved.

    Args:
        paths: Ordered list of JSON paths (chronological order).

    Returns:
        List of :class:`PocRecord` in paths order.

    Raises:
        ValueError: If a file is missing the required ``panel_flags`` key.
    """
    records: list[PocRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, Any] = json.load(fh)
        if "panel_flags" not in data:
            raise ValueError(f"{p}: missing 'panel_flags' key.")
        pf = data["panel_flags"]
        flags = PocPanelFlags(
            mtb=bool(pf.get("mtb", False)),
            cdiff=bool(pf.get("cdiff", False)),
            sars2=bool(pf.get("sars2", False)),
            flu_a=bool(pf.get("flu_a", False)),
            gas=bool(pf.get("gas", False)),
        )
        date = data.get("date") or None
        records.append(PocRecord(sample_id=p.stem, flags=flags, date=date))
    return records


def write_trend_json(trend: PocTrend, path: Path) -> None:
    """Write a :class:`PocTrend` to a JSON file.

    Args:
        trend: Trend object to serialise.
        path: Output path; parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: PocTrend, path: Path) -> None:
    """Write the per-sample timeline of a :class:`PocTrend` to a CSV file.

    Args:
        trend: Trend object whose records will be written.
        path: Output path; parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["sample_id", "date", "mtb", "cdiff", "sars2", "flu_a", "gas", "alert_level"]
        )
        for r in trend.records:
            lvl = assess_poc_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.mtb),
                    int(r.flags.cdiff),
                    int(r.flags.sars2),
                    int(r.flags.flu_a),
                    int(r.flags.gas),
                    lvl,
                ]
            )

r"""Respiratory RT-LAMP panel outbreak trend analysis for point-of-care settings.

Converts a chronological sequence of three-target respiratory LAMP panel
results (hRSV, Influenza A, SARS-CoV-2) into a :class:`PocRespTrend` that
captures the surveillance trajectory (EMERGING / RESOLVING / STABLE_CLEAR /
STABLE_ENDEMIC / INSUFFICIENT_DATA), dual-wave co-circulation, and
season-onset events -- the key signals for clinic capacity planning and
antiviral stockpile decisions.

This module mirrors :mod:`lamp_forge.poc_trend` for the focused respiratory
sub-panel (Recipe 34 in the cookbook): a portable BioVind-style device at a
rural clinic, urgent-care centre, or telemedicine hub accumulates three-target
results over a respiratory season and needs a single-command trend summary.

Targets tracked
---------------
    hRSV  - Human RSV (N gene, pan-RSV-A/B)          -- supportive / nirsevimab
    IAV   - Influenza A (M gene, pan-IAV)              -- oseltamivir within 48 h
    SARS2 - SARS-CoV-2 (N gene)                        -- nirmatrelvir / remdesivir

Key events detected
-------------------
    ``iav_newly_detected``:
        IAV absent in all prior samples but positive in the most recent one.
        Signals the start of an influenza season at this site; antiviral
        pre-positioning and post-exposure prophylaxis planning apply.

    ``dual_wave_active``:
        Both IAV and SARS-CoV-2 detected in the majority of the most-recent
        half of samples.  A co-epidemic ('twindemic') increases demand for
        two antiviral classes simultaneously and may produce clinically
        ambiguous presentations requiring both tests.

Typical workflow::

    lamp-forge poc-respiratory-risk --iav \
        --out-json results/clinic-A/2026-01.json
    lamp-forge poc-respiratory-risk --iav --sars2 \
        --out-json results/clinic-A/2026-02.json
    lamp-forge poc-respiratory-risk --sars2 \
        --out-json results/clinic-A/2026-03.json

    lamp-forge poc-respiratory-trend \
        --resp-result results/clinic-A/2026-01.json \
        --resp-result results/clinic-A/2026-02.json \
        --resp-result results/clinic-A/2026-03.json \
        --out-json results/clinic-A/resp_trend_Q1.json

CSV format (header required)::

    sample_id,date,hrsv,iav,sars2
    Clinic-2026-01,2026-01-15,0,1,0
    Clinic-2026-02,2026-02-15,0,1,1
    Clinic-2026-03,2026-03-15,0,0,1

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    Shi T et al. (2017) RSV global burden. Lancet 390:946-958.
        doi:10.1016/S0140-6736(17)30938-8
    Guo Z et al. (2021) Co-circulation of SARS-CoV-2 and influenza.
        J Infect 82(4):e1-e3. doi:10.1016/j.jinf.2021.01.028
    Nie K et al. (2017) RT-LAMP for hRSV. Arch Virol 162(4):1077-1083.
        doi:10.1007/s00705-016-3182-3
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

from lamp_forge.poc_respiratory_risk import PocRespAlertLevel, PocRespFlags, assess_poc_resp_risk

# ---------------------------------------------------------------------------
# Public enumerations
# ---------------------------------------------------------------------------


class PocRespTrendDirection(StrEnum):
    """Clinic-level surveillance trajectory for the respiratory RT-LAMP panel.

    Attributes:
        EMERGING: Pathogen burden increasing, or IAV detected for the first
            time (season onset).  Antiviral pre-positioning and enhanced
            surveillance recommended.
        RESOLVING: Pathogen(s) detected in early samples but absent in the
            most-recent sample.  Respiratory season appears to be waning.
        STABLE_CLEAR: All samples consistently negative.  Not in an active
            respiratory wave.
        STABLE_ENDEMIC: Consistent low-level pathogen positivity without
            escalation -- typical of an ongoing respiratory season.
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
class PocRespRecord:
    """One time-point in a respiratory RT-LAMP monitoring series.

    Attributes:
        sample_id: Identifier for this sampling event (e.g. a date string,
            batch number, or clinic pseudonym for aggregate results).
        flags: Three-target positivity flags from the respiratory panel.
        date: Optional ISO-8601 date string (YYYY-MM-DD).  Used for display
            and export only; chronological ordering is determined by list
            position.
    """

    sample_id: str
    flags: PocRespFlags
    date: str | None = None


@dataclass(frozen=True, slots=True)
class PocRespTrend:
    """Clinic surveillance trend derived from a series of respiratory LAMP results.

    Attributes:
        n_samples: Number of records in the time series.
        direction: Overall surveillance trajectory.
        records: Chronological tuple of records used to compute the trend.
        iav_newly_detected: True when IAV was absent in all preceding samples
            and positive in the most-recent sample -- a season-onset event
            signalling the start of an influenza wave.
        dual_wave_active: True when both IAV and SARS-CoV-2 were detected in
            the majority of the most-recent half of samples -- a co-epidemic
            ('twindemic') pattern.
        worst_alert_level: Highest :class:`~lamp_forge.poc_respiratory_risk\
.PocRespAlertLevel` observed across the entire time series.
        interpretation: Human-readable summary of the trend pattern.
        recommended_action: Actionable clinic guidance for this trajectory.
    """

    n_samples: int
    direction: PocRespTrendDirection
    records: tuple[PocRespRecord, ...]
    iav_newly_detected: bool
    dual_wave_active: bool
    worst_alert_level: PocRespAlertLevel
    interpretation: str
    recommended_action: str

    def to_dict(self) -> dict[str, Any]:
        """Serialise to a JSON-safe plain dict.

        Returns:
            Dict with all trend fields; per-sample timeline is a list of
            dicts each containing sample_id, date, three flag booleans, and
            the computed alert_level for that sample.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "iav_newly_detected": self.iav_newly_detected,
            "dual_wave_active": self.dual_wave_active,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "hrsv": r.flags.hrsv,
                    "iav": r.flags.iav,
                    "sars2": r.flags.sars2,
                    "alert_level": assess_poc_resp_risk(r.flags).alert_level.value,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

_LEVEL_ORDER: dict[PocRespAlertLevel, int] = {
    PocRespAlertLevel.NEGATIVE: 0,
    PocRespAlertLevel.LOW: 1,
    PocRespAlertLevel.MODERATE: 2,
    PocRespAlertLevel.HIGH: 3,
}


def _parse_bool(value: str) -> bool:
    return value.strip().lower() in ("1", "true", "yes")


def _worst_level(records: list[PocRespRecord]) -> PocRespAlertLevel:
    best = PocRespAlertLevel.NEGATIVE
    for r in records:
        lvl = assess_poc_resp_risk(r.flags).alert_level
        if _LEVEL_ORDER[lvl] > _LEVEL_ORDER[best]:
            best = lvl
    return best


def _burden(flags: PocRespFlags) -> int:
    return sum([flags.hrsv, flags.iav, flags.sars2])


def _any_positive(flags: PocRespFlags) -> bool:
    return flags.hrsv or flags.iav or flags.sars2


def _build_trend_text(
    direction: PocRespTrendDirection,
    iav_newly_detected: bool,
    dual_wave_active: bool,
) -> tuple[str, str]:
    """Construct (interpretation, recommended_action) for the detected trajectory.

    Args:
        direction: Computed surveillance trajectory.
        iav_newly_detected: IAV first seen in the last sample.
        dual_wave_active: IAV and SARS-CoV-2 co-dominant in recent samples.

    Returns:
        Two strings: (interpretation, recommended_action).
    """
    if direction is PocRespTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data: at least two samples are required to determine "
            "a respiratory surveillance trend.",
            "Continue sampling and rerun trend analysis once a second result is available.",
        )

    if direction is PocRespTrendDirection.STABLE_CLEAR:
        return (
            "All respiratory differential panel samples are negative "
            "(hRSV, IAV, SARS-CoV-2). No active respiratory viral circulation "
            "detected at this clinic.",
            "Maintain routine surveillance schedule. Confirm each negative run "
            "with a sample-adequacy control. No outbreak response required.",
        )

    if direction is PocRespTrendDirection.RESOLVING:
        return (
            "Respiratory pathogen positivity is declining across the sampling "
            "series. The respiratory wave appears to be waning.",
            "Maintain current infection-control measures and continue surveillance "
            "until two consecutive fully-negative results are achieved. Review "
            "antiviral stockpile consumption and replenish if needed.",
        )

    if iav_newly_detected:
        return (
            "Influenza A (pan-IAV) detected in the most-recent respiratory panel "
            "result with no IAV signal in prior samples. This marks the local "
            "onset of an influenza wave. Early-season detection is the optimal "
            "window for antiviral treatment and post-exposure prophylaxis.",
            "Activate influenza-season protocol. Ensure oseltamivir stock is "
            "adequate for anticipated case volume. Implement droplet precautions "
            "for all acute respiratory illness (ARI) presentations. Prioritise "
            "antiviral prescribing for high-risk patients within 48 h of symptom "
            "onset. Notify local public health if this represents the first "
            "seasonal IAV detection in your surveillance network.",
        )

    if direction is PocRespTrendDirection.STABLE_ENDEMIC:
        dual_note = (
            " Both IAV and SARS-CoV-2 are co-circulating (dual-wave / 'twindemic' "
            "pattern) in recent samples, increasing demand for two antiviral classes "
            "simultaneously."
            if dual_wave_active
            else ""
        )
        return (
            "Consistent respiratory viral positivity without clear escalation -- "
            "typical of an ongoing seasonal outbreak." + dual_note,
            "Maintain heightened infection-control precautions (droplet/contact "
            "precautions, cohorting of ARI patients, staff PPE). Review antiviral "
            "prescribing for eligible patients. Increase sampling frequency if "
            "positivity rate rises.",
        )

    # EMERGING
    dual_note = (
        " Dual-wave (IAV + SARS-CoV-2) co-epidemic is contributing to "
        "the escalation -- stock both antiviral classes."
        if dual_wave_active
        else ""
    )
    return (
        "Respiratory pathogen burden is increasing across the sampling series. "
        "A new or escalating outbreak is indicated." + dual_note,
        "Activate clinic respiratory outbreak protocol. Implement cohort "
        "isolation for all ARI patients. Ensure antiviral availability for "
        "high-risk patients. Increase testing frequency and consider contacting "
        "the local public health authority if the surge exceeds clinic capacity.",
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_poc_resp_trend(records: list[PocRespRecord]) -> PocRespTrend:
    """Compute a respiratory surveillance trend from an ordered series of records.

    Uses a split-half burden comparison to classify the trajectory.  IAV
    season-onset overrides the generic EMERGING label with more specific
    season-start guidance.

    Args:
        records: Chronologically ordered list of :class:`PocRespRecord` objects
            (oldest first).  Must contain at least one record.

    Returns:
        :class:`PocRespTrend` with direction, special-event flags, worst alert
        level, and human-readable guidance.

    Raises:
        ValueError: If ``records`` is empty.

    Example::

        from lamp_forge.poc_respiratory_risk import PocRespFlags
        from lamp_forge.poc_respiratory_trend import PocRespRecord, analyse_poc_resp_trend

        recs = [
            PocRespRecord("Jan", PocRespFlags(hrsv=False, iav=True, sars2=False)),
            PocRespRecord("Feb", PocRespFlags(hrsv=False, iav=True, sars2=True)),
            PocRespRecord("Mar", PocRespFlags(hrsv=False, iav=False, sars2=False)),
        ]
        trend = analyse_poc_resp_trend(recs)
        print(trend.direction)        # PocRespTrendDirection.RESOLVING
        print(trend.dual_wave_active) # True
    """
    if not records:
        raise ValueError("analyse_poc_resp_trend requires at least one record.")

    n = len(records)
    worst = _worst_level(records)

    if n < 2:
        return PocRespTrend(
            n_samples=n,
            direction=PocRespTrendDirection.INSUFFICIENT_DATA,
            records=tuple(records),
            iav_newly_detected=False,
            dual_wave_active=False,
            worst_alert_level=worst,
            interpretation=(
                "Insufficient data: at least two samples are required to determine "
                "a respiratory surveillance trend."
            ),
            recommended_action=(
                "Continue sampling and rerun trend analysis once a second result is available."
            ),
        )

    # IAV season-onset: absent in all prior samples, present in most recent.
    previous_any_iav = any(r.flags.iav for r in records[:-1])
    iav_newly_detected = records[-1].flags.iav and not previous_any_iav

    # Dual-wave: both IAV and SARS-CoV-2 in majority of recent half.
    mid = max(n // 2, 1)
    recent_half = records[mid:]
    dual_count = sum(1 for r in recent_half if r.flags.iav and r.flags.sars2)
    dual_wave_active = dual_count > len(recent_half) / 2

    # All-negative shortcut.
    if all(not _any_positive(r.flags) for r in records):
        interp, action = _build_trend_text(PocRespTrendDirection.STABLE_CLEAR, False, False)
        return PocRespTrend(
            n_samples=n,
            direction=PocRespTrendDirection.STABLE_CLEAR,
            records=tuple(records),
            iav_newly_detected=False,
            dual_wave_active=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
        )

    # Split-half burden comparison.
    early_half = records[:mid]
    early_burden = sum(_burden(r.flags) for r in early_half) / len(early_half)
    recent_burden = sum(_burden(r.flags) for r in recent_half) / len(recent_half)
    delta = recent_burden - early_burden

    if delta >= 0.5 or iav_newly_detected:
        direction = PocRespTrendDirection.EMERGING
    elif (not _any_positive(records[-1].flags) and early_burden > 0) or delta <= -0.5:
        direction = PocRespTrendDirection.RESOLVING
    else:
        direction = PocRespTrendDirection.STABLE_ENDEMIC

    interp, action = _build_trend_text(direction, iav_newly_detected, dual_wave_active)
    return PocRespTrend(
        n_samples=n,
        direction=direction,
        records=tuple(records),
        iav_newly_detected=iav_newly_detected,
        dual_wave_active=dual_wave_active,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
    )


# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

_REQUIRED_CSV_COLUMNS = frozenset({"sample_id", "hrsv", "iav", "sars2"})


def records_from_csv(path: Path) -> list[PocRespRecord]:
    """Parse a monitoring CSV into an ordered list of :class:`PocRespRecord`.

    Required header columns: ``sample_id``, ``hrsv``, ``iav``, ``sars2``.
    An optional ``date`` column (ISO-8601) is preserved.  Boolean columns
    accept ``0``/``1`` or ``true``/``false`` (case-insensitive).

    Args:
        path: Path to the CSV file.

    Returns:
        List of :class:`PocRespRecord` in file order.

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
        records: list[PocRespRecord] = []
        for row in reader:
            flags = PocRespFlags(
                hrsv=_parse_bool(row["hrsv"]),
                iav=_parse_bool(row["iav"]),
                sars2=_parse_bool(row["sars2"]),
            )
            date = row.get("date") or None
            records.append(PocRespRecord(sample_id=row["sample_id"], flags=flags, date=date))
    return records


def records_from_resp_json(paths: list[Path]) -> list[PocRespRecord]:
    """Parse ``lamp-forge poc-respiratory-risk --out-json`` files into records.

    Each file must contain a ``panel_flags`` dict with keys ``hrsv``, ``iav``,
    ``sars2``.  The file stem is used as the sample ID; an optional top-level
    ``date`` string is preserved.

    Args:
        paths: Ordered list of JSON paths (chronological order, oldest first).

    Returns:
        List of :class:`PocRespRecord` in paths order.

    Raises:
        ValueError: If a file is missing the required ``panel_flags`` key.
    """
    records: list[PocRespRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, Any] = json.load(fh)
        if "panel_flags" not in data:
            raise ValueError(
                f"{p}: missing 'panel_flags' key. "
                "Generate it with: lamp-forge poc-respiratory-risk --out-json <path>"
            )
        pf = data["panel_flags"]
        flags = PocRespFlags(
            hrsv=bool(pf.get("hrsv", False)),
            iav=bool(pf.get("iav", False)),
            sars2=bool(pf.get("sars2", False)),
        )
        date = data.get("date") or None
        records.append(PocRespRecord(sample_id=p.stem, flags=flags, date=date))
    return records


def write_trend_json(trend: PocRespTrend, path: Path) -> None:
    """Write a :class:`PocRespTrend` to a JSON file.

    Args:
        trend: Trend object to serialise.
        path: Output path; parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: PocRespTrend, path: Path) -> None:
    """Write the per-sample timeline of a :class:`PocRespTrend` to a CSV file.

    Columns: sample_id, date, hrsv, iav, sars2, alert_level.

    Args:
        trend: Trend object whose records will be written.
        path: Output path; parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "hrsv", "iav", "sars2", "alert_level"])
        for r in trend.records:
            lvl = assess_poc_resp_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.hrsv),
                    int(r.flags.iav),
                    int(r.flags.sars2),
                    lvl,
                ]
            )

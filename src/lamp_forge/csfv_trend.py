r"""CSFV NS5B RT-LAMP surveillance trend analysis for on-farm and border-post monitoring.

In BioVind-style farm biosecurity, a single CSFV NS5B RT-LAMP result tells you
whether Classical Swine Fever Virus RNA is present in one sample right now.
Trend analysis across consecutive monitoring intervals reveals whether a detected
CSFV event is newly emerging, represents an ongoing outbreak, or is being
controlled -- the critical decision variable for herd quarantine management and
national-authority reporting timelines.

CSFV is a WOAH-listed Tier 1 notifiable disease; any positive result is already
CRITICAL.  Trend analysis adds the **temporal dimension** that field veterinarians
and national authorities need to distinguish:

* A first detection (initiate emergency response).
* An ongoing outbreak (is containment failing?).
* Apparent containment (can movement restrictions be reviewed?).
* A persistently clean site (routine surveillance continuing).

This module converts a chronological sequence of NS5B RT-LAMP results into a
:class:`CsfvTrend` that captures:

* Overall surveillance **direction** (NEWLY_DETECTED / ONGOING_OUTBREAK /
  CONTROLLED / STABLE_CLEAR / INSUFFICIENT_DATA).
* **First-detection** and **controlled** event flags -- the WOAH-notification
  and quarantine-review inflection points.
* **Sustained outbreak flag** -- NS5B positive in the majority of samples in the
  series, suggesting the virus is actively circulating rather than a transient
  or false positive detection.

Typical workflow::

    lamp-forge csfv-risk --ns5b --context on_farm \
        --out-json results/farm-A/2026-01.json
    lamp-forge csfv-risk --no-ns5b --context on_farm \
        --out-json results/farm-A/2026-02.json

    lamp-forge csfv-trend \
        --csfv-result results/farm-A/2026-01.json \
        --csfv-result results/farm-A/2026-02.json \
        --out-json results/farm-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,ns5b_positive,context
    FarmA-2026-01,2026-01-15,1,on_farm
    FarmA-2026-02,2026-02-15,0,on_farm

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional.  The ``context`` column defaults to
``on_farm`` if absent.

References:
    Wang A, Cai X, Chen H et al. (2018) Rapid and sensitive detection of
        classical swine fever virus by reverse transcription loop-mediated
        isothermal amplification.
        Vet Microbiol 216:66-71. doi:10.1016/j.vetmic.2018.02.003
    Rodriguez-Prieto V et al. (2017) Real-time and conventional RT-LAMP for
        detection of classical swine fever virus.
        J Virol Methods 250:7-12. doi:10.1016/j.jviromet.2017.09.014
    WOAH (2022). Classical swine fever. Terrestrial Manual, Chapter 3.8.3.
    Moennig V et al. (2003) Clinical signs and epidemiology of classical swine
        fever: a review of new knowledge. Vet J 165(1):11-20.
        doi:10.1016/S1090-0233(02)00112-0
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.csfv_risk import (
    CsfvAlertLevel,
    CsfvContext,
    CsfvFlags,
    assess_csfv_risk,
)

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

_LEVEL_ORDER: list[CsfvAlertLevel] = [
    CsfvAlertLevel.NEGATIVE,
    CsfvAlertLevel.CRITICAL,
]

# Fraction of positive samples required to flag sustained_outbreak.
_SUSTAINED_THRESHOLD = 0.75
# Minimum samples required to apply the sustained_outbreak flag.
_SUSTAINED_MIN_N = 3


@dataclass(frozen=True, slots=True)
class CsfvRecord:
    """One chronological CSFV NS5B RT-LAMP result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"FarmA-2026-01"``).
        flags: NS5B result and deployment context for this interval.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the sample
            collection date.  Used only for display; ordering is based on
            sequence position in the input list.
    """

    sample_id: str
    flags: CsfvFlags
    date: str | None = None


class CsfvTrendDirection(StrEnum):
    """CSFV NS5B RT-LAMP surveillance trend direction derived from a time series.

    Attributes:
        NEWLY_DETECTED: NS5B positive in the most recent sample after being
            absent in all prior samples.  Triggers immediate WOAH notification,
            emergency herd quarantine, and trace-back/trace-forward investigation.
        ONGOING_OUTBREAK: NS5B positive in the most recent sample and in at
            least one prior sample.  Active CSFV circulation; containment may
            be failing.  Mandatory continued quarantine and intensified sampling.
        CONTROLLED: NS5B was positive in at least one prior sample but is
            negative in the most recent sample.  Apparent containment; maintain
            quarantine and confirm with repeat sampling before any review.
        STABLE_CLEAR: All samples in the series are NS5B negative.  CSFV not
            detected throughout the monitoring period.
        INSUFFICIENT_DATA: Fewer than two samples; a trend cannot be computed.
    """

    NEWLY_DETECTED = "NEWLY_DETECTED"
    ONGOING_OUTBREAK = "ONGOING_OUTBREAK"
    CONTROLLED = "CONTROLLED"
    STABLE_CLEAR = "STABLE_CLEAR"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class CsfvTrend:
    """CSFV NS5B RT-LAMP surveillance trajectory assessment for a time series.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        newly_detected: True if NS5B is positive in the most recent sample
            after being absent in all prior samples -- the WOAH-notification
            and emergency-quarantine inflection point.
        controlled: True if NS5B was positive in at least one prior sample but
            is negative in the most recent sample, suggesting containment may
            be working.
        sustained_outbreak: True when three or more records are present and
            >= 75% of samples are NS5B positive, indicating active CSFV
            circulation rather than a transient or isolated detection event.
        worst_alert_level: Highest :class:`~lamp_forge.csfv_risk.CsfvAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: CsfvTrendDirection
    newly_detected: bool
    controlled: bool
    sustained_outbreak: bool
    worst_alert_level: CsfvAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[CsfvRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample NS5B result and alert level in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "newly_detected": self.newly_detected,
            "controlled": self.controlled,
            "sustained_outbreak": self.sustained_outbreak,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "ns5b_positive": r.flags.ns5b_positive,
                    "context": r.flags.context.value,
                    "alert_level": assess_csfv_risk(r.flags).alert_level.value,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[CsfvRecord, ...],
) -> tuple[CsfvTrendDirection, bool, bool, bool]:
    """Compute trend direction and key event flags.

    Args:
        records: At least two :class:`CsfvRecord` in chronological order.

    Returns:
        Tuple of (direction, newly_detected, controlled, sustained_outbreak).
    """
    n = len(records)
    last_positive = records[-1].flags.ns5b_positive
    positive_in_prior = any(r.flags.ns5b_positive for r in records[:-1])

    newly_detected = last_positive and not positive_in_prior
    controlled = positive_in_prior and not last_positive

    n_positive = sum(1 for r in records if r.flags.ns5b_positive)
    sustained_outbreak = n >= _SUSTAINED_MIN_N and n_positive / n >= _SUSTAINED_THRESHOLD

    if all(not r.flags.ns5b_positive for r in records):
        return (
            CsfvTrendDirection.STABLE_CLEAR,
            newly_detected,
            controlled,
            sustained_outbreak,
        )

    if newly_detected:
        return (
            CsfvTrendDirection.NEWLY_DETECTED,
            newly_detected,
            controlled,
            sustained_outbreak,
        )

    if controlled:
        return (
            CsfvTrendDirection.CONTROLLED,
            newly_detected,
            controlled,
            sustained_outbreak,
        )

    # Last sample positive and at least one prior positive -> ONGOING_OUTBREAK.
    return (
        CsfvTrendDirection.ONGOING_OUTBREAK,
        newly_detected,
        controlled,
        sustained_outbreak,
    )


def _build_trend_text(
    direction: CsfvTrendDirection,
    n_samples: int,
    newly_detected: bool,
    controlled: bool,
    sustained_outbreak: bool,
    worst_level: CsfvAlertLevel,
    records: tuple[CsfvRecord, ...],
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a CSFV trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        newly_detected: True if NS5B first appeared in the last sample.
        controlled: True if NS5B was present before but is now absent.
        sustained_outbreak: True if >= 75% of >= 3 samples positive.
        worst_level: Highest alert level observed across all records.
        records: Chronological records (used for context inference).

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is CsfvTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the "
            "CSFV NS5B RT-LAMP surveillance trajectory.",
            "Collect the next scheduled monitoring sample and re-run "
            "csfv-trend after the second sample is available.",
        )

    n_positive = sum(1 for r in records if r.flags.ns5b_positive)
    last_ctx = records[-1].flags.context.value if records else "unknown"

    if direction is CsfvTrendDirection.NEWLY_DETECTED:
        interp = (
            f"NS5B RT-LAMP newly positive in the most recent sample "
            f"(context: {last_ctx}) after being negative in all "
            f"{n_samples - 1} prior interval(s). "
            "This is the first-detection event that requires immediate mandatory "
            "reporting under WOAH notification obligations for Classical Swine Fever. "
            "All in-contact pigs must be considered presumptively exposed."
        )
        action = (
            "CRITICAL -- CSFV first detection: initiate emergency herd quarantine "
            "immediately -- prohibit all live-pig, pork-product, and vehicle movement "
            "on and off the premises. "
            "Notify the national veterinary authority (e.g. USDA APHIS VS, national "
            "competent authority) and the state or provincial veterinarian TODAY -- "
            "CSFV is a WOAH-listed disease and same-day reporting is mandatory. "
            "Collect fresh EDTA blood and tonsil swabs from all sick and in-contact "
            "animals and submit to an approved reference laboratory for confirmatory "
            "RT-PCR, virus isolation, and serology. "
            "Do not euthanase or depopulate animals before receiving instructions from "
            "the national authority; confirmatory testing is required to trigger the "
            "official response."
        )
        return interp, action

    if direction is CsfvTrendDirection.ONGOING_OUTBREAK:
        n_pos_str = f"{n_positive}/{n_samples}"
        if sustained_outbreak:
            interp = (
                f"NS5B RT-LAMP positive in {n_pos_str} monitoring intervals "
                f"(context: {last_ctx}). "
                "The high positivity rate indicates a sustained CSFV outbreak with "
                "active virus circulation at this site. "
                "CSFV spreads rapidly via nose-to-nose contact, contaminated feed "
                "and water, and fomites; a sustained high positivity rate suggests "
                "the outbreak source has not been eliminated."
            )
            action = (
                "Maintain strict herd quarantine: all pig, product, and vehicle "
                "movement must remain suspended pending national authority instructions. "
                "Escalate to emergency depopulation assessment if not already initiated "
                "by the national authority. "
                "Review the trace-back/trace-forward investigation for undetected "
                "exposure links. "
                "Intensify sampling frequency (daily if feasible) and submit all "
                "NS5B-positive samples to the OIE Reference Laboratory for CSF for "
                "genotyping and strain characterisation. "
                "Disinfect all vehicles, equipment, and personnel pathways between "
                "sampling events."
            )
        else:
            interp = (
                f"NS5B RT-LAMP positive in {n_pos_str} monitoring intervals "
                f"(context: {last_ctx}), indicating ongoing CSFV detection at this "
                "site. The outbreak is still active; containment measures must remain "
                "in place."
            )
            action = (
                "Maintain herd quarantine and prohibit all movement. "
                "Continue sampling at the scheduled frequency and report all results "
                "to the national veterinary authority. "
                "Ensure full compliance with the national authority's response protocol "
                "-- depopulation, stamping-out, and zonal restriction decisions rest "
                "with the official authority, not the field team. "
                "Submit NS5B-positive samples to the reference laboratory for "
                "confirmatory RT-PCR and sequence genotyping."
            )
        return interp, action

    if direction is CsfvTrendDirection.CONTROLLED:
        interp = (
            f"NS5B RT-LAMP positive in at least one prior interval across "
            f"{n_samples} samples but negative in the most recent sample "
            f"(context: {last_ctx}). "
            "Apparent containment: CSFV RNA is no longer detectable in the most "
            "recent sample, suggesting that quarantine and response measures may be "
            "reducing virus load. "
            "However, a single negative RT-LAMP result is insufficient to confirm "
            "eradication; the high analytical sensitivity of RT-LAMP means a negative "
            "does not exclude low-level residual virus."
        )
        action = (
            "Maintain quarantine and LAMP surveillance at the scheduled or intensified "
            "frequency. "
            "Obtain at least three consecutive negative NS5B LAMP results before "
            "initiating confirmatory serology or requesting a controlled movement "
            "review from the national authority. "
            "Do not lift animal movement restrictions on the basis of a single "
            "negative LAMP result. "
            "Notify the national veterinary authority of the improving status and "
            "confirm their requirements for official clearance testing (typically "
            "RT-PCR on pooled blood samples and neutralisation serology)."
        )
        return interp, action

    # STABLE_CLEAR
    interp = (
        f"All {n_samples} NS5B RT-LAMP samples negative. "
        "CSFV not detected throughout the monitoring period. "
        "This result is consistent with a CSFV-free status for the monitored "
        "premises during this surveillance window."
    )
    action = (
        "Maintain the routine surveillance schedule. "
        "Before introducing new animals from outside the monitored herd, confirm "
        "their CSFV status with the supplier or require quarantine sampling on "
        "arrival. "
        "Repeat NS5B LAMP surveillance promptly if any swine exhibiting acute "
        "febrile illness, hemorrhagic signs, or unexplained mortality are observed. "
        "Ensure an extraction-positive control (e.g. internal RNA spike) is "
        "included in each LAMP run to confirm reagent and extraction integrity."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_csfv_trend(records: Sequence[CsfvRecord]) -> CsfvTrend:
    """Compute a CSFV NS5B RT-LAMP surveillance trajectory from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`CsfvRecord` (oldest first).
            Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`CsfvTrend` with trajectory direction, key event flags, and
        a field-actionable recommendation.

    Example::

        from lamp_forge.csfv_risk import CsfvContext, CsfvFlags
        from lamp_forge.csfv_trend import CsfvRecord, analyse_csfv_trend

        recs = [
            CsfvRecord("Farm-Jan", CsfvFlags(
                ns5b_positive=False, context=CsfvContext.ON_FARM)),
            CsfvRecord("Farm-Feb", CsfvFlags(
                ns5b_positive=True, context=CsfvContext.ON_FARM)),
        ]
        trend = analyse_csfv_trend(recs)
        print(trend.direction)       # CsfvTrendDirection.NEWLY_DETECTED
        print(trend.newly_detected)  # True
    """
    if not records:
        raise ValueError("records must contain at least one CsfvRecord")

    recs = tuple(records)
    n = len(recs)

    worst: CsfvAlertLevel = CsfvAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_csfv_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    if n < 2:
        interp, action = _build_trend_text(
            CsfvTrendDirection.INSUFFICIENT_DATA, n, False, False, False, worst, recs
        )
        return CsfvTrend(
            n_samples=n,
            direction=CsfvTrendDirection.INSUFFICIENT_DATA,
            newly_detected=False,
            controlled=False,
            sustained_outbreak=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, newly_detected, controlled, sustained_outbreak = _compute_direction(recs)
    interp, action = _build_trend_text(
        direction, n, newly_detected, controlled, sustained_outbreak, worst, recs
    )

    return CsfvTrend(
        n_samples=n,
        direction=direction,
        newly_detected=newly_detected,
        controlled=controlled,
        sustained_outbreak=sustained_outbreak,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[CsfvRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`CsfvRecord`.

    Expected columns (order-independent, header required):

        ``sample_id``    -- short identifier (required).
        ``ns5b_positive``-- bool: 0/1 or true/false (required).
        ``date``         -- ISO date string (optional; leave blank or omit column).
        ``context``      -- one of the :class:`~lamp_forge.csfv_risk.CsfvContext`
                           values (optional; defaults to ``on_farm``).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`CsfvRecord` in file order.
    """
    required = {"sample_id", "ns5b_positive"}
    result: list[CsfvRecord] = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"Empty CSV file: {path}")
        field_set = {f.strip().lower() for f in reader.fieldnames}
        missing = required - field_set
        if missing:
            raise ValueError(f"CSV missing required column(s): {', '.join(sorted(missing))}")
        for i, row in enumerate(reader, start=2):
            try:
                sample_id = row["sample_id"].strip()
                if not sample_id:
                    raise ValueError("sample_id must not be empty")
                date_raw = row.get("date", "").strip() or None
                ns5b = _is_truthy(row["ns5b_positive"])
                raw_ctx = row.get("context", "").strip() or CsfvContext.ON_FARM.value
                try:
                    ctx = CsfvContext(raw_ctx.lower())
                except ValueError:
                    valid = ", ".join(c.value for c in CsfvContext)
                    raise ValueError(
                        f"Row {i}: invalid context {raw_ctx!r}. Valid values: {valid}"
                    ) from None
                flags = CsfvFlags(ns5b_positive=ns5b, context=ctx)
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(CsfvRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_csfv_json(paths: Sequence[Path]) -> list[CsfvRecord]:
    """Load a sequence of csfv-risk JSON outputs as :class:`CsfvRecord`.

    Each JSON file (produced by ``lamp-forge csfv-risk --out-json``) supplies
    the ``ns5b_result`` and ``context`` fields for one sample.  The
    ``sample_id`` defaults to the file stem; an optional top-level ``"date"``
    key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``ns5b_result`` key.

    Returns:
        List of :class:`CsfvRecord` in the supplied order.
    """
    result: list[CsfvRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        if "ns5b_result" not in data:
            raise ValueError(
                f"JSON file {p} is missing the 'ns5b_result' key. "
                "Generate it with: lamp-forge csfv-risk --out-json <path>"
            )
        ns5b = str(data["ns5b_result"]).lower() == "positive"
        raw_ctx = str(data.get("context", CsfvContext.ON_FARM.value))
        try:
            ctx = CsfvContext(raw_ctx.lower())
        except ValueError:
            ctx = CsfvContext.ON_FARM
        flags = CsfvFlags(ns5b_positive=ns5b, context=ctx)
        date = str(data["date"]) if "date" in data else None
        result.append(CsfvRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: CsfvTrend, path: Path) -> None:
    """Write the CSFV trend assessment to a JSON file.

    Args:
        trend: :class:`CsfvTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: CsfvTrend, path: Path) -> None:
    """Write the chronological NS5B timeline to a CSV file.

    Each row represents one monitoring interval with the NS5B result, context,
    and per-sample alert level from
    :func:`~lamp_forge.csfv_risk.assess_csfv_risk`.

    Args:
        trend: :class:`CsfvTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "ns5b_positive", "context", "alert_level"])
        for r in trend.records:
            lvl = assess_csfv_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.ns5b_positive),
                    r.flags.context.value,
                    lvl,
                ]
            )

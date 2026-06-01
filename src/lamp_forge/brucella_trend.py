r"""Brucella IS711 LAMP surveillance trend analysis for on-farm and POC monitoring.

In BioVind-style farm biosecurity, a single Brucella IS711 LAMP result tells you
whether Brucella genus DNA is present in one sample right now.  Trend analysis
across consecutive monitoring intervals reveals whether a detected brucellosis
event is being controlled, is persisting, or has cleared -- the critical decision
variable for herd quarantine management and confirmatory testing timelines.

This module converts a chronological sequence of IS711 LAMP results into a
:class:`BrucellaTrend` that captures:

* Overall surveillance **direction** (NEWLY_DETECTED / PERSISTENT /
  RESOLVING / STABLE_CLEAR / INSUFFICIENT_DATA).
* **First-detection** and **clearance** event flags -- the WOAH-notification and
  quarantine-release inflection points.
* **Persistent infection flag** -- IS711 positive in the majority of samples in the
  series, suggesting an ongoing source animal rather than a transient
  environmental contamination event.

Typical workflow::

    lamp-forge brucella-risk --is711 --context cattle \
        --out-json results/herd-A/2026-01.json
    lamp-forge brucella-risk --no-is711 --context cattle \
        --out-json results/herd-A/2026-02.json

    lamp-forge brucella-trend \
        --brucella-result results/herd-A/2026-01.json \
        --brucella-result results/herd-A/2026-02.json \
        --out-json results/herd-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,is711_positive,context
    HerdA-2026-01,2026-01-15,1,cattle
    HerdA-2026-02,2026-02-15,0,cattle

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional.  The ``context`` column defaults to
``cattle`` if absent.

References:
    Abd El Wahed A et al. (2015) Loop-mediated isothermal amplification on a
        chip for detection of Brucella spp. in milk samples.
        PLoS Negl Trop Dis 9(10):e0004124. doi:10.1371/journal.pntd.0004124
    Al-Dahouk S et al. (2010) Evaluation of Brucella IS711-based loop-mediated
        isothermal amplification for human brucellosis.
        Ann Clin Microbiol Antimicrob 9:15. doi:10.1186/1476-0711-9-15
    OIE/WOAH (2021). Bovine brucellosis. Terrestrial Manual, Chapter 3.1.2.
    Godfroid J et al. (2011). Brucellosis at the animal/ecosystem/human interface
        at the beginning of the 21st century.
        Prev Vet Med 102(2):118-131. doi:10.1016/j.prevetmed.2011.04.007
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.brucella_risk import (
    BrucellaAlertLevel,
    BrucellaContext,
    BrucellaFlags,
    assess_brucella_risk,
)

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

_LEVEL_ORDER: list[BrucellaAlertLevel] = [
    BrucellaAlertLevel.NEGATIVE,
    BrucellaAlertLevel.LOW,
    BrucellaAlertLevel.HIGH,
    BrucellaAlertLevel.CRITICAL,
]

# Fraction of positive samples required to flag persistent_infection_likely.
_PERSISTENCE_THRESHOLD = 0.75
# Minimum samples required to apply the persistence flag.
_PERSISTENCE_MIN_N = 3


@dataclass(frozen=True, slots=True)
class BrucellaRecord:
    """One chronological Brucella IS711 LAMP result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"HerdA-2026-01"``).
        flags: IS711 result and sample context for this interval.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the sample
            collection date.  Used only for display; ordering is based on
            sequence position in the input list.
    """

    sample_id: str
    flags: BrucellaFlags
    date: str | None = None


class BrucellaTrendDirection(StrEnum):
    """Brucella IS711 surveillance trend direction derived from a time series.

    Attributes:
        NEWLY_DETECTED: IS711 positive in the most recent sample after being
            absent in all prior samples.  Triggers WOAH notification and
            immediate herd quarantine (livestock contexts) or public-health
            reporting (human context).
        PERSISTENT: IS711 positive in the most recent sample and in at least
            one prior sample.  Ongoing infection; a persistently infected
            source animal is the most likely explanation when multiple
            consecutive positives are observed.
        RESOLVING: IS711 was positive in at least one prior sample but is
            negative in the most recent sample.  The herd may be clearing;
            confirmatory serology required before lifting quarantine.
        STABLE_CLEAR: All samples in the series are IS711 negative.  No
            Brucella detected throughout the monitoring period.
        INSUFFICIENT_DATA: Fewer than two samples; a trend cannot be computed.
    """

    NEWLY_DETECTED = "NEWLY_DETECTED"
    PERSISTENT = "PERSISTENT"
    RESOLVING = "RESOLVING"
    STABLE_CLEAR = "STABLE_CLEAR"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class BrucellaTrend:
    """Brucella IS711 surveillance trajectory assessment for a time series.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        newly_detected: True if IS711 is positive in the most recent sample
            after being absent in all prior samples -- the WOAH-notification
            and quarantine inflection point.
        cleared: True if IS711 was positive in at least one prior sample but
            is negative in the most recent sample.
        persistent_infection_likely: True when three or more records are
            present and >= 75 % of samples are IS711 positive, suggesting a
            persistently shedding source animal rather than a transient
            environmental contamination event.
        worst_alert_level: Highest :class:`~lamp_forge.brucella_risk.BrucellaAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: BrucellaTrendDirection
    newly_detected: bool
    cleared: bool
    persistent_infection_likely: bool
    worst_alert_level: BrucellaAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[BrucellaRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample IS711 result and alert level in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "newly_detected": self.newly_detected,
            "cleared": self.cleared,
            "persistent_infection_likely": self.persistent_infection_likely,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "is711_positive": r.flags.is711_positive,
                    "context": r.flags.context.value,
                    "alert_level": assess_brucella_risk(r.flags).alert_level.value,
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
    records: tuple[BrucellaRecord, ...],
) -> tuple[BrucellaTrendDirection, bool, bool, bool]:
    """Compute trend direction and key event flags.

    Args:
        records: At least two :class:`BrucellaRecord` in chronological order.

    Returns:
        Tuple of (direction, newly_detected, cleared, persistent_infection_likely).
    """
    n = len(records)
    last_positive = records[-1].flags.is711_positive
    positive_in_prior = any(r.flags.is711_positive for r in records[:-1])

    newly_detected = last_positive and not positive_in_prior
    cleared = positive_in_prior and not last_positive

    n_positive = sum(1 for r in records if r.flags.is711_positive)
    persistent_infection_likely = (
        n >= _PERSISTENCE_MIN_N and n_positive / n >= _PERSISTENCE_THRESHOLD
    )

    if all(not r.flags.is711_positive for r in records):
        return (
            BrucellaTrendDirection.STABLE_CLEAR,
            newly_detected,
            cleared,
            persistent_infection_likely,
        )

    if newly_detected:
        return (
            BrucellaTrendDirection.NEWLY_DETECTED,
            newly_detected,
            cleared,
            persistent_infection_likely,
        )

    if cleared:
        return (
            BrucellaTrendDirection.RESOLVING,
            newly_detected,
            cleared,
            persistent_infection_likely,
        )

    # Last sample positive and at least one prior positive -> PERSISTENT.
    return (
        BrucellaTrendDirection.PERSISTENT,
        newly_detected,
        cleared,
        persistent_infection_likely,
    )


def _build_trend_text(
    direction: BrucellaTrendDirection,
    n_samples: int,
    newly_detected: bool,
    cleared: bool,
    persistent_infection_likely: bool,
    worst_level: BrucellaAlertLevel,
    records: tuple[BrucellaRecord, ...],
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a Brucella trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        newly_detected: True if IS711 first appeared in the last sample.
        cleared: True if IS711 was present before but is now absent.
        persistent_infection_likely: True if >= 75 % of >= 3 samples positive.
        worst_level: Highest alert level observed across all records.
        records: Chronological records (used for context inference).

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is BrucellaTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the "
            "Brucella IS711 surveillance trajectory.",
            "Collect the next scheduled monitoring sample and re-run "
            "brucella-trend after the second sample is available.",
        )

    n_positive = sum(1 for r in records if r.flags.is711_positive)
    last_ctx = records[-1].flags.context.value if records else "unknown"

    if direction is BrucellaTrendDirection.NEWLY_DETECTED:
        interp = (
            f"IS711 LAMP newly positive in the most recent sample (context: {last_ctx}) "
            f"after being negative in all {n_samples - 1} prior interval(s). "
            "This is the first-detection event that requires immediate regulatory "
            "action for livestock contexts (WOAH-listed disease) or public-health "
            "notification for human brucellosis."
        )
        if worst_level is BrucellaAlertLevel.CRITICAL:
            action = (
                "CRITICAL -- human Brucella first detection: initiate "
                "doxycycline + rifampin combination therapy immediately after "
                "collecting confirmatory samples. "
                "Notify the national public health authority today. "
                "Interview the patient for animal contact history to identify the "
                "zoonotic source."
            )
        else:
            action = (
                "Immediately quarantine the herd and cease all animal movement. "
                "Notify the national veterinary authority today (e.g. USDA APHIS VS). "
                "Submit blood or milk samples from the positive animal(s) for "
                "confirmatory serology (RBPT + CFT or i-ELISA) and culture at an "
                "accredited reference laboratory. "
                "Advise farm workers and veterinarians of zoonotic exposure risk."
            )
        return interp, action

    if direction is BrucellaTrendDirection.PERSISTENT:
        n_pos_str = f"{n_positive}/{n_samples}"
        if persistent_infection_likely:
            interp = (
                f"IS711 LAMP positive in {n_pos_str} monitoring intervals. "
                "The high positivity rate suggests a persistently infected source "
                "animal (PI animal) that is continuously shedding Brucella into the "
                "herd environment. "
                "Environmental contamination alone does not sustain repeated positive "
                "results at this level."
            )
            action = (
                "Conduct whole-herd serology (RBPT + CFT/i-ELISA) to identify "
                "the PI animal(s) and remove them from the herd. "
                "Maintain quarantine and LAMP surveillance until three consecutive "
                "negative results are achieved after removal of seropositive animals. "
                "Disinfect all shared equipment, pens, and parturition areas. "
                "Do not lift quarantine without confirmatory clearance from an "
                "accredited reference laboratory."
            )
        else:
            interp = (
                f"IS711 LAMP positive in {n_pos_str} monitoring intervals, "
                "indicating ongoing brucellosis in the monitored herd or premises. "
                "A persistently infected source animal is possible but whole-herd "
                "serology is needed to confirm."
            )
            action = (
                "Maintain herd quarantine. "
                "Escalate to whole-herd serology (RBPT + CFT/i-ELISA) to identify "
                "infected individuals. "
                "Continue LAMP surveillance at the scheduled frequency. "
                "Consult the national veterinary authority regarding depopulation or "
                "test-and-remove protocols applicable in your jurisdiction."
            )
        return interp, action

    if direction is BrucellaTrendDirection.RESOLVING:
        if cleared:
            interp = (
                f"IS711 LAMP positive in at least one prior interval across "
                f"{n_samples} samples, but negative in the most recent sample. "
                "The herd may be clearing, but a single negative LAMP result is "
                "insufficient to confirm eradication -- IS711 sensitivity means a "
                "negative result reduces but does not eliminate the probability of "
                "remaining infection."
            )
            action = (
                "Maintain quarantine and LAMP monitoring. "
                "Obtain at least three consecutive negative LAMP results before "
                "initiating confirmatory serology (RBPT + CFT/i-ELISA) for official "
                "clearance. "
                "Do not lift animal movement restrictions on the basis of a single "
                "negative LAMP result. "
                "Notify the national veterinary authority of the improving status and "
                "confirm their clearance testing requirements."
            )
        else:
            interp = (
                f"Brucella IS711 positivity is declining across the {n_samples} "
                "monitoring intervals."
            )
            action = (
                "Maintain quarantine and LAMP surveillance. "
                "Confirm clearance with at least three consecutive negative results "
                "before proceeding to official serology-based clearance testing."
            )
        return interp, action

    # STABLE_CLEAR
    interp = (
        f"All {n_samples} IS711 LAMP samples negative. "
        "No Brucella detected throughout the monitoring period."
    )
    action = (
        "Maintain the routine surveillance schedule. "
        "Before introducing new animals, confirm their Brucella status with "
        "serology or LAMP. "
        "Repeat LAMP surveillance promptly if any reproductive failures "
        "(abortions, weak neonates) or other high-risk events occur."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_brucella_trend(records: Sequence[BrucellaRecord]) -> BrucellaTrend:
    """Compute a Brucella IS711 surveillance trajectory from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`BrucellaRecord` (oldest first).
            Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`BrucellaTrend` with trajectory direction, key event flags, and
        a field-actionable recommendation.

    Example::

        from lamp_forge.brucella_risk import BrucellaContext, BrucellaFlags
        from lamp_forge.brucella_trend import BrucellaRecord, analyse_brucella_trend

        recs = [
            BrucellaRecord("Herd-Jan", BrucellaFlags(
                is711_positive=False, context=BrucellaContext.CATTLE)),
            BrucellaRecord("Herd-Feb", BrucellaFlags(
                is711_positive=True, context=BrucellaContext.CATTLE)),
        ]
        trend = analyse_brucella_trend(recs)
        print(trend.direction)      # BrucellaTrendDirection.NEWLY_DETECTED
        print(trend.newly_detected) # True
    """
    if not records:
        raise ValueError("records must contain at least one BrucellaRecord")

    recs = tuple(records)
    n = len(recs)

    worst: BrucellaAlertLevel = BrucellaAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_brucella_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    if n < 2:
        interp, action = _build_trend_text(
            BrucellaTrendDirection.INSUFFICIENT_DATA, n, False, False, False, worst, recs
        )
        return BrucellaTrend(
            n_samples=n,
            direction=BrucellaTrendDirection.INSUFFICIENT_DATA,
            newly_detected=False,
            cleared=False,
            persistent_infection_likely=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, newly_detected, cleared, persistent_infection_likely = _compute_direction(recs)
    interp, action = _build_trend_text(
        direction, n, newly_detected, cleared, persistent_infection_likely, worst, recs
    )

    return BrucellaTrend(
        n_samples=n,
        direction=direction,
        newly_detected=newly_detected,
        cleared=cleared,
        persistent_infection_likely=persistent_infection_likely,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[BrucellaRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`BrucellaRecord`.

    Expected columns (order-independent, header required):

        ``sample_id``      -- short identifier (required).
        ``is711_positive`` -- bool: 0/1 or true/false (required).
        ``date``           -- ISO date string (optional; leave blank or omit column).
        ``context``        -- one of the :class:`~lamp_forge.brucella_risk.BrucellaContext`
                              values (optional; defaults to ``cattle``).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`BrucellaRecord` in file order.
    """
    required = {"sample_id", "is711_positive"}
    result: list[BrucellaRecord] = []
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
                is711 = _is_truthy(row["is711_positive"])
                raw_ctx = row.get("context", "").strip() or BrucellaContext.CATTLE.value
                try:
                    ctx = BrucellaContext(raw_ctx.lower())
                except ValueError:
                    valid = ", ".join(c.value for c in BrucellaContext)
                    raise ValueError(
                        f"Row {i}: invalid context {raw_ctx!r}. Valid values: {valid}"
                    ) from None
                flags = BrucellaFlags(is711_positive=is711, context=ctx)
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(BrucellaRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_brucella_json(paths: Sequence[Path]) -> list[BrucellaRecord]:
    """Load a sequence of brucella-risk JSON outputs as :class:`BrucellaRecord`.

    Each JSON file (produced by ``lamp-forge brucella-risk --out-json``) supplies
    the ``is711_result`` and ``context`` fields for one sample.  The ``sample_id``
    defaults to the file stem; an optional top-level ``"date"`` key is used if
    present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``is711_result`` key.

    Returns:
        List of :class:`BrucellaRecord` in the supplied order.
    """
    result: list[BrucellaRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        if "is711_result" not in data:
            raise ValueError(
                f"JSON file {p} is missing the 'is711_result' key. "
                "Generate it with: lamp-forge brucella-risk --out-json <path>"
            )
        is711 = str(data["is711_result"]).lower() == "positive"
        raw_ctx = str(data.get("context", BrucellaContext.CATTLE.value))
        try:
            ctx = BrucellaContext(raw_ctx.lower())
        except ValueError:
            ctx = BrucellaContext.CATTLE
        flags = BrucellaFlags(is711_positive=is711, context=ctx)
        date = str(data["date"]) if "date" in data else None
        result.append(BrucellaRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: BrucellaTrend, path: Path) -> None:
    """Write the Brucella trend assessment to a JSON file.

    Args:
        trend: :class:`BrucellaTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: BrucellaTrend, path: Path) -> None:
    """Write the chronological IS711 timeline to a CSV file.

    Each row represents one monitoring interval with the IS711 result, context,
    and per-sample alert level from
    :func:`~lamp_forge.brucella_risk.assess_brucella_risk`.

    Args:
        trend: :class:`BrucellaTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "is711_positive", "context", "alert_level"])
        for r in trend.records:
            lvl = assess_brucella_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.is711_positive),
                    r.flags.context.value,
                    lvl,
                ]
            )

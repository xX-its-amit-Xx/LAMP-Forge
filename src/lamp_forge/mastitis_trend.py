r"""Bovine mastitis LAMP panel surveillance trend analysis for dairy-herd monitoring.

In BioVind-style farm biosecurity, a single mastitis LAMP panel result tells
you which pathogens are present in one milk sample right now.  Trend analysis
across consecutive monitoring intervals reveals whether contagious mastitis
pathogens (*Staphylococcus aureus*, *Streptococcus agalactiae*) are being
eradicated, are persisting, or have newly appeared -- the critical decision
variable for herd quarantine management and milking-hygiene programmes.

This module converts a chronological sequence of bovine mastitis LAMP panel
results into a :class:`MastitisTrend` that captures:

* Overall herd biosecurity **direction** (NEWLY_DETECTED / PERSISTENT /
  RESOLVING / STABLE_CLEAR / STABLE_ENDEMIC / INSUFFICIENT_DATA).
* **Contagious-pathogen first-detection and clearance event flags** --
  the milking-hygiene programme inflection points.
* **Persistent contagious flag** -- contagious pathogen positive in the
  majority of samples in the series, suggesting an unidentified carrier
  cow rather than a transient introduction event.

Typical workflow::

    lamp-forge mastitis-risk --saur \
        --out-json results/herd-A/2026-01.json
    lamp-forge mastitis-risk --saur --sube \
        --out-json results/herd-A/2026-02.json
    lamp-forge mastitis-risk \
        --out-json results/herd-A/2026-03.json

    lamp-forge mastitis-trend \
        --mastitis-result results/herd-A/2026-01.json \
        --mastitis-result results/herd-A/2026-02.json \
        --mastitis-result results/herd-A/2026-03.json \
        --out-json results/herd-A/trend_Q1.json

CSV format (header required)::

    sample_id,date,saur,saga,sube,ecoli
    HerdA-2026-01,2026-01-15,1,0,0,0
    HerdA-2026-02,2026-02-15,1,0,1,0
    HerdA-2026-03,2026-03-15,0,0,0,0

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit or leave blank.

References:
    Hogeveen H et al. (2011) Economic aspects of mastitis: new developments.
        N Z Vet J 59(1):16-23. doi:10.1080/00480169.2011.552541
    Bradley AJ & Green MJ (2004) The importance and control of mastitis.
        Vet Clin N Am Food Anim Pract 20(3):469-491.
        doi:10.1016/j.cvfa.2004.06.010
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

from lamp_forge.mastitis_risk import (
    MastitisAlertLevel,
    MastitisPanelFlags,
    assess_mastitis_risk,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_LEVEL_ORDER: list[MastitisAlertLevel] = [
    MastitisAlertLevel.NEGATIVE,
    MastitisAlertLevel.LOW,
    MastitisAlertLevel.MODERATE,
    MastitisAlertLevel.HIGH,
    MastitisAlertLevel.CRITICAL,
]

# Fraction of samples with contagious pathogen required to flag persistent.
_PERSISTENCE_THRESHOLD = 0.75
# Minimum samples required to apply the persistence flag.
_PERSISTENCE_MIN_N = 3


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MastitisRecord:
    """One chronological bovine mastitis LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"HerdA-2026-01"``).
        flags: Four-pathogen positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the sample
            collection date.  Used only for display; ordering is based on
            sequence position in the input list.
    """

    sample_id: str
    flags: MastitisPanelFlags
    date: str | None = None


class MastitisTrendDirection(StrEnum):
    """Bovine mastitis surveillance trend direction derived from a time series.

    Attributes:
        NEWLY_DETECTED: A contagious pathogen (S. aureus or S. agalactiae)
            appears for the first time in the most recent sample after being
            absent in all prior samples.  Triggers immediate cow segregation
            and whole-herd LAMP screening.
        PERSISTENT: Contagious pathogen present in the most recent sample
            and in at least one prior sample.  Ongoing contagious mastitis
            problem; a carrier cow is likely the uncontrolled source.
        RESOLVING: Contagious pathogen was present in at least one prior
            sample but is absent in the most recent sample.  Eradication
            programme may be working; confirm with three consecutive negatives.
        STABLE_CLEAR: All samples entirely negative for all four targets.
            No mastitis pathogens detected throughout the monitoring period.
        STABLE_ENDEMIC: No contagious pathogen detected in any sample, but
            environmental pathogens (S. uberis or E. coli) recurring across
            at least two samples.  Indicates a husbandry deficiency rather
            than a milking-hygiene failure.
        INSUFFICIENT_DATA: Fewer than two samples; a trend cannot be computed.
    """

    NEWLY_DETECTED = "NEWLY_DETECTED"
    PERSISTENT = "PERSISTENT"
    RESOLVING = "RESOLVING"
    STABLE_CLEAR = "STABLE_CLEAR"
    STABLE_ENDEMIC = "STABLE_ENDEMIC"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class MastitisTrend:
    """Bovine mastitis surveillance trajectory assessment for a time series.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        contagious_newly_detected: True if a contagious pathogen (S. aureus
            or S. agalactiae) appears in the most recent sample after being
            absent in all prior samples.
        contagious_cleared: True if a contagious pathogen was present in at
            least one prior sample but is absent in the most recent sample.
        persistent_contagious_likely: True when three or more records are
            present and >= 75 % of samples have a contagious pathogen,
            suggesting an unidentified carrier cow rather than a transient
            introduction.
        worst_alert_level: Highest
            :class:`~lamp_forge.mastitis_risk.MastitisAlertLevel` observed
            across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: MastitisTrendDirection
    contagious_newly_detected: bool
    contagious_cleared: bool
    persistent_contagious_likely: bool
    worst_alert_level: MastitisAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[MastitisRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample positivity flags and alert level in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "contagious_newly_detected": self.contagious_newly_detected,
            "contagious_cleared": self.contagious_cleared,
            "persistent_contagious_likely": self.persistent_contagious_likely,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "saur": r.flags.saur,
                    "saga": r.flags.saga,
                    "sube": r.flags.sube,
                    "ecoli": r.flags.ecoli,
                    "alert_level": assess_mastitis_risk(r.flags).alert_level.value,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _has_contagious(flags: MastitisPanelFlags) -> bool:
    """Return True if any contagious mastitis pathogen (SAUR or SAGA) is positive."""
    return flags.saur or flags.saga


def _has_environmental(flags: MastitisPanelFlags) -> bool:
    """Return True if any environmental mastitis pathogen (SUBE or ECOLI) is positive."""
    return flags.sube or flags.ecoli


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[MastitisRecord, ...],
) -> tuple[MastitisTrendDirection, bool, bool, bool]:
    """Compute trend direction and key event flags.

    Args:
        records: At least two :class:`MastitisRecord` in chronological order.

    Returns:
        Tuple of (direction, contagious_newly_detected, contagious_cleared,
        persistent_contagious_likely).
    """
    n = len(records)
    last_contagious = _has_contagious(records[-1].flags)
    prior_contagious = any(_has_contagious(r.flags) for r in records[:-1])

    contagious_newly_detected = last_contagious and not prior_contagious
    contagious_cleared = prior_contagious and not last_contagious

    n_contagious = sum(1 for r in records if _has_contagious(r.flags))
    persistent_contagious_likely = (
        n >= _PERSISTENCE_MIN_N and n_contagious / n >= _PERSISTENCE_THRESHOLD
    )

    if not any(_has_contagious(r.flags) for r in records):
        n_env = sum(1 for r in records if _has_environmental(r.flags))
        if n_env >= 2:
            return (
                MastitisTrendDirection.STABLE_ENDEMIC,
                contagious_newly_detected,
                contagious_cleared,
                persistent_contagious_likely,
            )
        return (
            MastitisTrendDirection.STABLE_CLEAR,
            contagious_newly_detected,
            contagious_cleared,
            persistent_contagious_likely,
        )

    if contagious_newly_detected:
        return (
            MastitisTrendDirection.NEWLY_DETECTED,
            contagious_newly_detected,
            contagious_cleared,
            persistent_contagious_likely,
        )

    if contagious_cleared:
        return (
            MastitisTrendDirection.RESOLVING,
            contagious_newly_detected,
            contagious_cleared,
            persistent_contagious_likely,
        )

    return (
        MastitisTrendDirection.PERSISTENT,
        contagious_newly_detected,
        contagious_cleared,
        persistent_contagious_likely,
    )


def _build_trend_text(
    direction: MastitisTrendDirection,
    n_samples: int,
    contagious_newly_detected: bool,
    contagious_cleared: bool,
    persistent_contagious_likely: bool,
    records: tuple[MastitisRecord, ...],
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for a mastitis trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        contagious_newly_detected: True if contagious pathogen first appeared in
            the last sample.
        contagious_cleared: True if contagious pathogen was present before but is
            now absent.
        persistent_contagious_likely: True if >= 75 % of >= 3 samples have a
            contagious pathogen.
        records: Chronological records (used for pathogen-label inference).

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is MastitisTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the "
            "bovine mastitis surveillance trajectory.",
            "Collect the next scheduled monitoring sample and re-run "
            "mastitis-trend after the second sample is available.",
        )

    n_contagious = sum(1 for r in records if _has_contagious(r.flags))
    last_flags = records[-1].flags

    if direction is MastitisTrendDirection.NEWLY_DETECTED:
        labels: list[str] = []
        if last_flags.saur:
            labels.append("S. aureus (nuc)")
        if last_flags.saga:
            labels.append("S. agalactiae (cfb)")
        label_str = " and ".join(labels) if labels else "contagious pathogen"
        interp = (
            f"{label_str} newly detected in the most recent sample after being "
            f"absent in all {n_samples - 1} prior interval(s). "
            "This is a first-detection event for a contagious mastitis pathogen "
            "that spreads cow-to-cow during milking; immediate action is required "
            "to prevent herd-level dissemination."
        )
        action = (
            "Segregate the LAMP-positive cow immediately to the end of the milking "
            "order or to a separate milking unit. Conduct a whole-herd LAMP screen "
            "within 48 hours to identify all contagious-pathogen carriers. "
            "Enforce pre- and post-milking teat disinfection and inspect teat-cup "
            "liner integrity. Notify the herd veterinarian to plan an eradication "
            "programme (targeted dry-cow therapy or culling for S. aureus; "
            "intramammary penicillin for S. agalactiae)."
        )
        return interp, action

    if direction is MastitisTrendDirection.PERSISTENT:
        n_pos_str = f"{n_contagious}/{n_samples}"
        if persistent_contagious_likely:
            interp = (
                f"Contagious mastitis pathogen detected in {n_pos_str} monitoring "
                "intervals. The high frequency suggests an unidentified carrier cow "
                "as the persistent source. Environmental contamination alone does not "
                "sustain repeated contagious-pathogen detection at this rate."
            )
            action = (
                "Escalate to a whole-herd individual-cow LAMP screen to identify "
                "and remove the carrier cow(s). "
                "Maintain milking segregation of known positive cows. "
                "Review and enforce post-milking teat disinfection and teat-cup "
                "back-flush between cows. "
                "Consult the herd veterinarian for a formal eradication programme "
                "using test-and-cull or test-and-treat for S. aureus; intramammary "
                "penicillin eradication is effective for S. agalactiae."
            )
        else:
            interp = (
                f"Contagious mastitis pathogen detected in {n_pos_str} monitoring "
                "intervals, indicating ongoing contagious mastitis in the herd. "
                "A persistent carrier cow is possible but whole-herd screening is "
                "required to confirm."
            )
            action = (
                "Maintain milking-order segregation of LAMP-positive cows. "
                "Conduct a whole-herd LAMP or culture screen to identify all carriers. "
                "Continue mastitis monitoring at the scheduled frequency. "
                "Consult the herd veterinarian for a treatment or culling protocol "
                "appropriate to the detected pathogen."
            )
        return interp, action

    if direction is MastitisTrendDirection.RESOLVING:
        interp = (
            f"Contagious mastitis pathogen was detected in at least one prior interval "
            f"across {n_samples} samples, but is absent in the most recent sample. "
            "The eradication programme may be working, but a single negative result "
            "is insufficient to confirm full herd clearance."
        )
        action = (
            "Maintain milking-hygiene protocols and LAMP monitoring. "
            "Continue surveillance until at least three consecutive negative results "
            "are recorded before declaring the herd free of contagious mastitis. "
            "Verify that previously positive cows have been treated or culled. "
            "Maintain post-milking teat disinfection and teat-cup liner replacement "
            "schedule."
        )
        return interp, action

    if direction is MastitisTrendDirection.STABLE_ENDEMIC:
        interp = (
            f"No contagious mastitis pathogen detected across {n_samples} samples. "
            "Environmental mastitis pathogens (S. uberis or E. coli) detected in "
            "multiple intervals, indicating a recurring husbandry deficiency rather "
            "than a milking-hygiene failure."
        )
        action = (
            "Review cubicle and bedding management: replace contaminated organic "
            "bedding; ensure correct lying-surface depth; disinfect before refilling. "
            "Assess teat-end condition and teat-spray coverage. "
            "Consider blanket dry-cow therapy with teat sealant at dry-off to reduce "
            "the dry-period new infection rate. "
            "No immediate milking-order segregation is required, but treat clinical "
            "cases promptly and record somatic cell count trends."
        )
        return interp, action

    # STABLE_CLEAR
    interp = (
        f"All {n_samples} mastitis LAMP samples negative. "
        "No mastitis pathogens detected throughout the monitoring period."
    )
    action = (
        "Maintain the routine mastitis surveillance schedule. "
        "Confirm sample adequacy with a 16S rRNA or internal extraction control. "
        "Continue post-milking teat disinfection and quarterly mastitis monitoring."
    )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_mastitis_trend(records: Sequence[MastitisRecord]) -> MastitisTrend:
    """Compute a bovine mastitis surveillance trajectory from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`MastitisRecord` (oldest
            first).  Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`MastitisTrend` with trajectory direction, key event flags, and
        a field-actionable recommendation.

    Example::

        from lamp_forge.mastitis_risk import MastitisPanelFlags
        from lamp_forge.mastitis_trend import MastitisRecord, analyse_mastitis_trend

        recs = [
            MastitisRecord("Jan", MastitisPanelFlags(saur=False, saga=False,
                           sube=False, ecoli=False)),
            MastitisRecord("Feb", MastitisPanelFlags(saur=True, saga=False,
                           sube=False, ecoli=False)),
        ]
        trend = analyse_mastitis_trend(recs)
        print(trend.direction)                 # MastitisTrendDirection.NEWLY_DETECTED
        print(trend.contagious_newly_detected)  # True
    """
    if not records:
        raise ValueError("records must contain at least one MastitisRecord")

    recs = tuple(records)
    n = len(recs)

    worst: MastitisAlertLevel = MastitisAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_mastitis_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    if n < 2:
        interp, action = _build_trend_text(
            MastitisTrendDirection.INSUFFICIENT_DATA, n, False, False, False, recs
        )
        return MastitisTrend(
            n_samples=n,
            direction=MastitisTrendDirection.INSUFFICIENT_DATA,
            contagious_newly_detected=False,
            contagious_cleared=False,
            persistent_contagious_likely=False,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, contagious_newly_detected, contagious_cleared, persistent_contagious_likely = (
        _compute_direction(recs)
    )
    interp, action = _build_trend_text(
        direction,
        n,
        contagious_newly_detected,
        contagious_cleared,
        persistent_contagious_likely,
        recs,
    )

    return MastitisTrend(
        n_samples=n,
        direction=direction,
        contagious_newly_detected=contagious_newly_detected,
        contagious_cleared=contagious_cleared,
        persistent_contagious_likely=persistent_contagious_likely,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[MastitisRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`MastitisRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier (required).
        ``saur``      -- bool: 0/1 or true/false (required).
        ``saga``      -- bool: 0/1 or true/false (required).
        ``sube``      -- bool: 0/1 or true/false (optional; defaults to 0).
        ``ecoli``     -- bool: 0/1 or true/false (optional; defaults to 0).
        ``date``      -- ISO date string (optional; leave blank or omit column).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`MastitisRecord` in file order.
    """
    required = {"sample_id", "saur", "saga"}
    result: list[MastitisRecord] = []
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
                saur = _is_truthy(row["saur"])
                saga = _is_truthy(row["saga"])
                sube = _is_truthy(row.get("sube", "0"))
                ecoli = _is_truthy(row.get("ecoli", "0"))
                flags = MastitisPanelFlags(saur=saur, saga=saga, sube=sube, ecoli=ecoli)
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(MastitisRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_mastitis_json(paths: Sequence[Path]) -> list[MastitisRecord]:
    """Load a sequence of mastitis-risk JSON outputs as :class:`MastitisRecord`.

    Each JSON file (produced by ``lamp-forge mastitis-risk --out-json``) supplies
    the ``panel_flags`` dict for one sample.  The ``sample_id`` defaults to the
    file stem; an optional top-level ``"date"`` key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key or
            the value is not a dict.

    Returns:
        List of :class:`MastitisRecord` in the supplied order.
    """
    result: list[MastitisRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        if "panel_flags" not in data:
            raise ValueError(
                f"JSON file {p} is missing the 'panel_flags' key. "
                "Generate it with: lamp-forge mastitis-risk --out-json <path>"
            )
        pf = data["panel_flags"]
        if not isinstance(pf, dict):
            raise ValueError(f"JSON file {p}: 'panel_flags' must be a dict")
        flags = MastitisPanelFlags(
            saur=bool(pf.get("saur", False)),
            saga=bool(pf.get("saga", False)),
            sube=bool(pf.get("sube", False)),
            ecoli=bool(pf.get("ecoli", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(MastitisRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: MastitisTrend, path: Path) -> None:
    """Write the mastitis trend assessment to a JSON file.

    Args:
        trend: :class:`MastitisTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: MastitisTrend, path: Path) -> None:
    """Write the chronological mastitis timeline to a CSV file.

    Each row represents one monitoring interval with the four target flags
    and per-sample alert level from
    :func:`~lamp_forge.mastitis_risk.assess_mastitis_risk`.

    Args:
        trend: :class:`MastitisTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "saur", "saga", "sube", "ecoli", "alert_level"])
        for r in trend.records:
            lvl = assess_mastitis_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.saur),
                    int(r.flags.saga),
                    int(r.flags.sube),
                    int(r.flags.ecoli),
                    lvl,
                ]
            )

r"""Souring-trajectory analysis for time-series oilfield MIC monitoring.

In oilfield MIC management, a single LAMP panel result indicates presence or
absence of each corrosion guild at one point in time.  Trend analysis across
multiple consecutive sample intervals reveals whether the microbial community
-- and therefore the corrosion risk -- is **improving**, **stable**, or
**deteriorating**, which is the critical decision variable for treatment
escalation and biocide programme optimisation.

This module converts a chronological sequence of five-guild LAMP panel results
into a :class:`SouringTrend` that captures:

* Overall souring **trajectory** (IMPROVING / STABLE_LOW / STABLE_HIGH /
  DETERIORATING / INSUFFICIENT_DATA).
* SRB **new-detection** and **clearance** events -- the inflection points that
  trigger immediate operational responses.
* **Treatment effectiveness** tracking: the count of sample intervals in which
  both SRB and NRB are positive indicates nitrate treatment is present but SRB
  is not suppressed (treatment under-dosed or community has adapted).
* Worst risk level observed across the time series.

Typical workflow::

    # Run the LAMP panel monthly and save results
    lamp-forge mic-risk --srb --irb \\
      --out-json results/W1/2026-01.json
    lamp-forge mic-risk --srb --irb --apb \\
      --out-json results/W1/2026-02.json
    lamp-forge mic-risk --srb --irb --apb \\
      --out-json results/W1/2026-03.json

    # Analyse the three-month trend
    lamp-forge souring-trend \\
      --mic-result results/W1/2026-01.json \\
      --mic-result results/W1/2026-02.json \\
      --mic-result results/W1/2026-03.json \\
      --out-json results/W1/trend_Q1.json

    # Or ingest a monitoring-spreadsheet CSV
    lamp-forge souring-trend \\
      --csv monitoring/W1_2026.csv \\
      --out-json results/W1/trend_2026.json

CSV format (header required)::

    sample_id,date,srb,mcr,irb,apb,nrb
    W1-2026-01,2026-01-15,1,0,0,0,0
    W1-2026-02,2026-02-15,1,0,1,0,0
    W1-2026-03,2026-03-15,1,0,1,1,0

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    Hubert C & Voordouw G (2007) Oil field souring control by nitrate-reducing
        bacteria.  Appl Environ Microbiol 73:2644-2652.
    Hamilton WA (2003) Microbially influenced corrosion as a model system for
        the study of metal microbe interactions.  Int Biodeterior Biodegrad
        51(3):151-155.
    Vigneron A et al. (2017) Comammox Nitrospira in oilfield MIC microbiomes.
        ISME J 11:2180-2194.
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.mic_risk import GuildFlags, MICRiskLevel, assess_mic_risk

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class GuildRecord:
    """One chronological LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"W1-2026-01"``).
        flags: Five-guild positivity flags for this sample.
        date: Optional ISO-format date string for the sample collection
            date (``"YYYY-MM-DD"``).  Used only for display and export;
            ordering is based on the sequence position in the input list.
    """

    sample_id: str
    flags: GuildFlags
    date: str | None = None


class TrendDirection(StrEnum):
    """Souring trajectory direction derived from a time series of guild flags.

    Attributes:
        IMPROVING: SRB detection rate is decreasing or SRB has cleared;
            current treatment programme appears effective.
        STABLE_LOW: SRB rarely or never detected; low sustained corrosion
            activity.
        STABLE_HIGH: SRB frequently detected without a clear downward trend;
            chronic souring requiring programme review.
        DETERIORATING: SRB detection rate is increasing or SRB has appeared
            for the first time; immediate treatment response required.
        INSUFFICIENT_DATA: Fewer than two samples; cannot compute a trend.
    """

    IMPROVING = "IMPROVING"
    STABLE_LOW = "STABLE_LOW"
    STABLE_HIGH = "STABLE_HIGH"
    DETERIORATING = "DETERIORATING"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class SouringTrend:
    """Souring-trajectory assessment for a time-series of oilfield panel results.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        srb_detection_count: Number of sample intervals with SRB positive.
        nrb_detection_count: Number of sample intervals with NRB positive
            (indicates nitrate-injection programme activity).
        srb_newly_detected: True when SRB is positive only in the most recent
            sample and was absent in all previous samples -- an inflection
            point requiring immediate intervention.
        srb_clearing: True when SRB was positive in the first sample and is
            negative in the most recent sample -- indicating successful
            suppression.
        treatment_underdosed_count: Number of samples in which both SRB and
            NRB are positive, implying nitrate treatment is present but
            insufficient to suppress sulfate-reduction.
        worst_risk_level: Highest :class:`~lamp_forge.mic_risk.MICRiskLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: TrendDirection
    srb_detection_count: int
    nrb_detection_count: int
    srb_newly_detected: bool
    srb_clearing: bool
    treatment_underdosed_count: int
    worst_risk_level: MICRiskLevel
    interpretation: str
    recommended_action: str
    records: tuple[GuildRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample guild flags in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "srb_detection_count": self.srb_detection_count,
            "nrb_detection_count": self.nrb_detection_count,
            "srb_newly_detected": self.srb_newly_detected,
            "srb_clearing": self.srb_clearing,
            "treatment_underdosed_count": self.treatment_underdosed_count,
            "worst_risk_level": self.worst_risk_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "srb": r.flags.srb,
                    "mcr": r.flags.mcr,
                    "irb": r.flags.irb,
                    "apb": r.flags.apb,
                    "nrb": r.flags.nrb,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Trend direction logic
# ---------------------------------------------------------------------------

# Risk level ordering for worst-level computation.
_LEVEL_ORDER: list[MICRiskLevel] = [
    MICRiskLevel.MINIMAL,
    MICRiskLevel.LOW,
    MICRiskLevel.MODERATE,
    MICRiskLevel.HIGH,
    MICRiskLevel.CRITICAL,
]


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[GuildRecord, ...],
) -> tuple[TrendDirection, bool, bool]:
    """Compute trend direction and key SRB inflection flags.

    The SRB guild is the primary indicator because it is the direct driver of
    H2S production and sulfide corrosion.  MCR, IRB, APB, and NRB contribute
    to risk scoring (via :func:`~lamp_forge.mic_risk.assess_mic_risk`) but the
    trend direction is keyed on SRB dynamics.

    The early half of the record sequence is compared to the recent half.
    A rate change of >=0.4 (40 percentage points) is required to classify a
    directional trend; smaller movements are classified as stable.

    Args:
        records: At least two :class:`GuildRecord` in chronological order.

    Returns:
        Tuple of (direction, srb_newly_detected, srb_clearing).
    """
    n = len(records)
    srb = [r.flags.srb for r in records]

    # SRB appeared for the first time in the most recent sample.
    srb_newly_detected = srb[-1] and all(not s for s in srb[:-1])
    # SRB was present at the start but is absent in the most recent sample.
    srb_clearing = srb[0] and not srb[-1]

    if srb_newly_detected:
        return TrendDirection.DETERIORATING, True, False
    if srb_clearing:
        return TrendDirection.IMPROVING, False, True

    # Compare SRB detection rate between the first half and the second half.
    mid = max(1, n // 2)
    n_early = mid
    n_recent = n - mid
    early_rate = sum(1 for s in srb[:mid] if s) / n_early
    recent_rate = sum(1 for s in srb[mid:] if s) / n_recent

    if recent_rate - early_rate >= 0.4:
        return TrendDirection.DETERIORATING, False, False
    if early_rate - recent_rate >= 0.4 and not srb[-1]:
        return TrendDirection.IMPROVING, False, False

    overall_rate = sum(1 for s in srb if s) / n
    if overall_rate >= 0.5:
        return TrendDirection.STABLE_HIGH, False, False
    return TrendDirection.STABLE_LOW, False, False


def _build_trend_text(
    direction: TrendDirection,
    n_samples: int,
    srb_count: int,
    nrb_count: int,
    underdosed: int,
    newly_detected: bool,
    clearing: bool,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for the trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples.
        srb_count: Number of samples with SRB positive.
        nrb_count: Number of samples with NRB positive.
        underdosed: Number of samples with both SRB and NRB positive.
        newly_detected: True if SRB appeared only in the last sample.
        clearing: True if SRB was positive in the first sample but is now
            negative.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is TrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine souring trajectory.",
            "Collect the next scheduled monitoring sample. Re-run souring-trend "
            "after the second sample is available.",
        )

    # -------------------------------------------------------------------
    # DETERIORATING
    # -------------------------------------------------------------------
    if direction is TrendDirection.DETERIORATING:
        if newly_detected:
            interp = (
                f"SRB newly detected in the most recent sample after "
                f"{n_samples - 1} consecutive negative interval(s). "
                "This is a souring-onset event requiring immediate response."
            )
            action = (
                "Apply biocide treatment immediately. Increase monitoring frequency "
                "to weekly until SRB returns to negative for two consecutive intervals. "
                "Review injection water quality and temperature excursions that may have "
                "promoted SRB establishment."
            )
        else:
            interp = (
                f"SRB detection rate is increasing across the {n_samples}-sample "
                f"time series ({srb_count} of {n_samples} intervals positive). "
                "Corrosion risk is rising."
            )
            action = (
                "Escalate biocide treatment (higher dose or increased frequency). "
                "Inspect pipeline segments for under-deposit pitting. "
                "If nitrate injection is in use, review residence time and dose."
            )
        if underdosed > 0:
            interp += (
                f" NRB co-detected in {underdosed} of these sample(s): nitrate "
                "treatment is active but insufficient to suppress SRB."
            )
            action += (
                " Increase nitrate injection dose -- current level is not achieving "
                "thermodynamic suppression of sulfate-reduction."
            )
        return interp, action

    # -------------------------------------------------------------------
    # IMPROVING
    # -------------------------------------------------------------------
    if direction is TrendDirection.IMPROVING:
        if clearing:
            interp = (
                f"SRB has cleared in the most recent sample after being positive "
                f"in earlier intervals. Current treatment programme appears effective "
                f"({srb_count} of {n_samples} intervals SRB positive, now negative)."
            )
            action = (
                "Maintain current biocide or nitrate-injection programme. "
                "Continue monthly monitoring to confirm sustained suppression. "
                "Do not reduce treatment until SRB has been negative for at least "
                "three consecutive intervals."
            )
        else:
            interp = (
                f"SRB detection rate is decreasing across the {n_samples}-sample "
                f"time series ({srb_count} of {n_samples} intervals positive). "
                "Treatment programme is effective."
            )
            action = (
                "Continue current treatment at the established dose and frequency. "
                "Maintain monitoring schedule; premature treatment reduction risks "
                "SRB rebound."
            )
        return interp, action

    # -------------------------------------------------------------------
    # STABLE_HIGH
    # -------------------------------------------------------------------
    if direction is TrendDirection.STABLE_HIGH:
        interp = (
            f"SRB consistently detected across the monitoring period "
            f"({srb_count} of {n_samples} intervals positive). "
            "Chronic souring; current treatment programme is not achieving suppression."
        )
        action = (
            "Treatment programme review required. Options: increase biocide dose, "
            "switch biocide chemistry (resistance adaptation), or introduce nitrate "
            "injection as a supplementary strategy to shift community towards "
            "NRB dominance. Consult corrosion engineer for treatment optimisation."
        )
        if underdosed > 0:
            interp += (
                f" NRB positive in {underdosed} of {n_samples} interval(s): nitrate "
                "injection is present but the SRB-NRB competition is not resolving "
                "in favour of NRB."
            )
            action += (
                " Nitrate dose is insufficient -- increase injection rate or "
                "consider combined nitrate + biocide protocol."
            )
        return interp, action

    # -------------------------------------------------------------------
    # STABLE_LOW
    # -------------------------------------------------------------------
    interp = (
        f"SRB consistently absent or rarely detected across the monitoring "
        f"period ({srb_count} of {n_samples} intervals positive). "
        "Low sustained corrosion activity."
    )
    action = (
        "Maintain current monitoring schedule. No treatment escalation required. "
        "Re-evaluate if process conditions change (temperature, water injection "
        "rate, sulfate concentration)."
    )
    if nrb_count > 0 and srb_count == 0:
        interp += (
            f" NRB detected in {nrb_count} interval(s) without SRB co-detection: "
            "nitrate injection is suppressing sulfate-reduction effectively."
        )
        action = (
            "Continue nitrate-injection programme at the current dose. "
            "Maintain monthly SRB monitoring to confirm sustained suppression."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_souring_trend(records: Sequence[GuildRecord]) -> SouringTrend:
    """Compute a souring trajectory assessment from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`GuildRecord` (oldest first).
            Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`SouringTrend` with trajectory direction, key event flags,
        and field-actionable recommendation.

    Example::

        from lamp_forge.mic_risk import GuildFlags
        from lamp_forge.souring_trend import GuildRecord, analyse_souring_trend

        recs = [
            GuildRecord("W1-Jan", GuildFlags(srb=False, mcr=False, irb=False,
                        apb=False, nrb=False), date="2026-01-15"),
            GuildRecord("W1-Feb", GuildFlags(srb=False, mcr=False, irb=False,
                        apb=False, nrb=False), date="2026-02-15"),
            GuildRecord("W1-Mar", GuildFlags(srb=True,  mcr=False, irb=False,
                        apb=False, nrb=False), date="2026-03-15"),
        ]
        trend = analyse_souring_trend(recs)
        print(trend.direction)      # TrendDirection.DETERIORATING
        print(trend.srb_newly_detected)  # True
    """
    if not records:
        raise ValueError("records must contain at least one GuildRecord")

    recs = tuple(records)
    n = len(recs)

    srb_count = sum(1 for r in recs if r.flags.srb)
    nrb_count = sum(1 for r in recs if r.flags.nrb)
    underdosed = sum(1 for r in recs if r.flags.srb and r.flags.nrb)

    # Worst risk level across all records.
    worst: MICRiskLevel = MICRiskLevel.MINIMAL
    for r in recs:
        lvl = assess_mic_risk(r.flags).risk_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    if n < 2:
        interp, action = _build_trend_text(
            TrendDirection.INSUFFICIENT_DATA, n, srb_count, nrb_count, underdosed, False, False
        )
        return SouringTrend(
            n_samples=n,
            direction=TrendDirection.INSUFFICIENT_DATA,
            srb_detection_count=srb_count,
            nrb_detection_count=nrb_count,
            srb_newly_detected=False,
            srb_clearing=False,
            treatment_underdosed_count=underdosed,
            worst_risk_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, newly_detected, clearing = _compute_direction(recs)
    interp, action = _build_trend_text(
        direction, n, srb_count, nrb_count, underdosed, newly_detected, clearing
    )

    return SouringTrend(
        n_samples=n,
        direction=direction,
        srb_detection_count=srb_count,
        nrb_detection_count=nrb_count,
        srb_newly_detected=newly_detected,
        srb_clearing=clearing,
        treatment_underdosed_count=underdosed,
        worst_risk_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[GuildRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`GuildRecord`.

    Expected columns (order-independent, header required):

        ``sample_id`` -- short identifier for the interval (required).
        ``date``      -- ISO date string (optional; leave blank or omit column).
        ``srb``, ``mcr``, ``irb``, ``apb``, ``nrb`` -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`GuildRecord` in file order.
    """
    required = {"sample_id", "srb", "mcr", "irb", "apb", "nrb"}
    records: list[GuildRecord] = []
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
                flags = GuildFlags(
                    srb=_is_truthy(row["srb"]),
                    mcr=_is_truthy(row["mcr"]),
                    irb=_is_truthy(row["irb"]),
                    apb=_is_truthy(row["apb"]),
                    nrb=_is_truthy(row["nrb"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            records.append(GuildRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return records


def records_from_mic_json(paths: Sequence[Path]) -> list[GuildRecord]:
    """Load a sequence of mic-risk JSON outputs as :class:`GuildRecord` list.

    Each JSON file (produced by ``lamp-forge mic-risk --out-json``) supplies
    the ``guild_flags`` dict for one sample.  The sample_id defaults to the
    file stem; an optional top-level ``"date"`` key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``guild_flags`` key.

    Returns:
        List of :class:`GuildRecord` in the supplied order.
    """
    records: list[GuildRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        gf = data.get("guild_flags")
        if not isinstance(gf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'guild_flags' dict. "
                "Generate it with: lamp-forge mic-risk --out-json <path>"
            )
        flags = GuildFlags(
            srb=bool(gf.get("srb", False)),
            mcr=bool(gf.get("mcr", False)),
            irb=bool(gf.get("irb", False)),
            apb=bool(gf.get("apb", False)),
            nrb=bool(gf.get("nrb", False)),
        )
        date = str(data["date"]) if "date" in data else None
        records.append(GuildRecord(sample_id=p.stem, flags=flags, date=date))
    return records


def write_trend_json(trend: SouringTrend, path: Path) -> None:
    """Write the trend assessment to a JSON file.

    Args:
        trend: :class:`SouringTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: SouringTrend, path: Path) -> None:
    """Write the chronological timeline of guild flags to a CSV file.

    Each row represents one monitoring interval with per-guild boolean columns
    and the per-sample risk level from
    :func:`~lamp_forge.mic_risk.assess_mic_risk`.

    Args:
        trend: :class:`SouringTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample_id", "date", "srb", "mcr", "irb", "apb", "nrb", "risk_level"])
        for r in trend.records:
            lvl = assess_mic_risk(r.flags).risk_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.srb),
                    int(r.flags.mcr),
                    int(r.flags.irb),
                    int(r.flags.apb),
                    int(r.flags.nrb),
                    lvl,
                ]
            )

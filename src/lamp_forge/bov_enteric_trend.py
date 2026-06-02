r"""Bovine neonatal calf diarrhea LAMP panel trend analysis.

Converts a chronological sequence of bovine neonatal enteric LAMP panel
results into a :class:`BovEntericTrend` that captures the herd-level neonatal
scours trajectory (EMERGING / STABLE_CLEAR / STABLE_ENDEMIC / RESOLVING /
INSUFFICIENT_DATA), ETEC K99 new-detection and clearance events, Salmonella
persistence (continuous zoonotic shedding), and the worst alert level -- the
key signals for immediate IV antibiotic escalation, isolation, and zoonotic
management decisions.

Pathogens tracked
-----------------
    ETEC_K99    faeG     ETEC K99/F5 -- highest mortality, IV antibiotics
    SALMONELLA  invA     Salmonella sp. -- bacteremia risk, zoonotic
    CRYPTO      18S/hsp  Cryptosporidium parvum -- zoonotic, halofuginone
    BCOV        N gene   Bovine coronavirus (enteric) -- supportive care
    ROTA_A      VP6      Bovine rotavirus A -- oral rehydration

Special events tracked
----------------------
* **ETEC newly detected** -- ETEC K99 (faeG) positive for the first time in
  the monitoring series.  A critical inflection event: IV antibiotic therapy
  must begin within 2 hours to prevent near-certain mortality in calves
  <4 days old.  Immediate isolation and pen sanitisation are required.
* **ETEC cleared** -- ETEC K99 was positive in prior intervals but is negative
  in the most recent sample.  Confirms treatment is working; surveillance
  must continue because environmental ETEC contamination persists in pen
  bedding and water troughs.
* **Salmonella persistent** -- Salmonella positive in every monitored interval.
  Suggests herd-level contamination or a persistent shedder; confirmatory
  culture and a herd-level investigation are warranted.
* **Zoonotic risk count** -- number of intervals in which Salmonella or
  Cryptosporidium is positive; both are notifiable food-safety pathogens.
* **Antibiotic required count** -- number of intervals in which ETEC K99 or
  Salmonella is detected; these are the only enteric panel targets for which
  antimicrobial therapy is standard of care.

Typical workflow::

    # Run the enteric panel on the calf group each week
    lamp-forge bov-enteric-risk --rota-a \
        --out-json results/calf-grp-1/week-01.json
    lamp-forge bov-enteric-risk --rota-a --bcov \
        --out-json results/calf-grp-1/week-02.json
    lamp-forge bov-enteric-risk --etec-k99 --rota-a \
        --out-json results/calf-grp-1/week-03.json

    # Analyse the three-week trend
    lamp-forge bov-enteric-trend \
        --enteric-result results/calf-grp-1/week-01.json \
        --enteric-result results/calf-grp-1/week-02.json \
        --enteric-result results/calf-grp-1/week-03.json \
        --out-json results/calf-grp-1/trend_wk1-3.json

    # Or ingest a monitoring-spreadsheet CSV
    lamp-forge bov-enteric-trend \
        --csv monitoring/calf_grp1_2026.csv \
        --out-json results/calf-grp-1/trend_2026.json

CSV format (header required)::

    sample_id,date,etec_k99,salmonella,crypto,bcov,rota_a
    CG1-week01,2026-01-07,0,0,0,0,1
    CG1-week02,2026-01-14,0,0,0,1,1
    CG1-week03,2026-01-21,1,0,0,1,1

Boolean columns accept 0/1 or true/false (case-insensitive).
The ``date`` column is optional; omit the column or leave values blank.

References:
    Cho Y & Yoon K (2014) An overview of calf diarrhea -- infectious etiology,
        diagnosis, and intervention. J Vet Sci 15(1):1-17.
        doi:10.4142/jvs.2014.15.1.1
    Gulliksen SM et al. (2009) Enteropathogens and risk factors for diarrhea
        in Norwegian dairy calves. J Dairy Sci 92(10):5057-5066.
        doi:10.3168/jds.2009-2080
    Smith DR (2012) Field diseases of beef calves. Vet Clin Food Anim
        28(1):57-82. doi:10.1016/j.cvfa.2012.01.001
"""

from __future__ import annotations

import csv
import json
from collections.abc import Sequence
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path

from lamp_forge.bov_enteric_risk import (
    BovEntericAlertLevel,
    BovEntericFlags,
    assess_bov_enteric_risk,
)

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

_LEVEL_ORDER: list[BovEntericAlertLevel] = [
    BovEntericAlertLevel.NEGATIVE,
    BovEntericAlertLevel.LOW,
    BovEntericAlertLevel.MODERATE,
    BovEntericAlertLevel.HIGH,
    BovEntericAlertLevel.CRITICAL,
]


@dataclass(frozen=True, slots=True)
class BovEntericRecord:
    """One chronological neonatal enteric LAMP panel result for trend analysis.

    Attributes:
        sample_id: Short identifier for this monitoring interval
            (e.g. ``"CG1-week01"``).
        flags: Five-target positivity flags for this sample.
        date: Optional ISO-format date string ``"YYYY-MM-DD"`` for the
            sample collection date.  Used only for display and export;
            ordering is based on sequence position in the input list.
    """

    sample_id: str
    flags: BovEntericFlags
    date: str | None = None


class BovEntericTrendDirection(StrEnum):
    """Neonatal calf diarrhea trend direction from a time series of panel results.

    Attributes:
        EMERGING: Pathogen burden increasing or ETEC K99 newly detected --
            outbreak risk rising; immediate escalation required.
        STABLE_CLEAR: All panels consistently negative.
        STABLE_ENDEMIC: Pathogens consistently detected without a clear trend;
            chronic scours requiring herd-management review.
        RESOLVING: Pathogen burden decreasing or ETEC K99 cleared in the most
            recent interval; control measures appear effective.
        INSUFFICIENT_DATA: Fewer than two samples; trend cannot be computed.
    """

    EMERGING = "EMERGING"
    STABLE_CLEAR = "STABLE_CLEAR"
    STABLE_ENDEMIC = "STABLE_ENDEMIC"
    RESOLVING = "RESOLVING"
    INSUFFICIENT_DATA = "INSUFFICIENT_DATA"


@dataclass(frozen=True, slots=True)
class BovEntericTrend:
    """Neonatal calf diarrhea trajectory for a time series of enteric panel results.

    Attributes:
        n_samples: Number of chronological records included.
        direction: Overall trend direction.
        etec_newly_detected: True if ETEC K99 was absent in all prior samples
            but positive in the most recent one -- an emergency event
            requiring IV antibiotic therapy within 2 hours.
        etec_cleared: True if ETEC K99 was positive in at least one earlier
            sample but is negative in the most recent sample.
        salmonella_persistent: True if Salmonella was positive in every
            monitored interval -- suggests herd-level contamination or a
            persistent shedder.
        zoonotic_risk_count: Number of intervals in which Salmonella or
            Cryptosporidium is positive.  Both are notifiable zoonotic
            pathogens.
        antibiotic_required_count: Number of intervals in which ETEC K99 or
            Salmonella is detected.  Antimicrobial therapy is indicated for
            both.
        worst_alert_level: Highest :class:`~lamp_forge.bov_enteric_risk.BovEntericAlertLevel`
            observed across all records.
        interpretation: Concise description of the trend pattern.
        recommended_action: Field-actionable recommendation for this trend.
        records: The chronological input records used to compute this trend.
    """

    n_samples: int
    direction: BovEntericTrendDirection
    etec_newly_detected: bool
    etec_cleared: bool
    salmonella_persistent: bool
    zoonotic_risk_count: int
    antibiotic_required_count: int
    worst_alert_level: BovEntericAlertLevel
    interpretation: str
    recommended_action: str
    records: tuple[BovEntericRecord, ...]

    def to_dict(self) -> dict[str, object]:
        """Serialise to a plain dict suitable for JSON output.

        Returns:
            Dict with all trend fields; the ``timeline`` list preserves the
            per-sample pathogen flags in chronological order.
        """
        return {
            "n_samples": self.n_samples,
            "direction": self.direction.value,
            "etec_newly_detected": self.etec_newly_detected,
            "etec_cleared": self.etec_cleared,
            "salmonella_persistent": self.salmonella_persistent,
            "zoonotic_risk_count": self.zoonotic_risk_count,
            "antibiotic_required_count": self.antibiotic_required_count,
            "worst_alert_level": self.worst_alert_level.value,
            "interpretation": self.interpretation,
            "recommended_action": self.recommended_action,
            "timeline": [
                {
                    "sample_id": r.sample_id,
                    "date": r.date,
                    "etec_k99": r.flags.etec_k99,
                    "salmonella": r.flags.salmonella,
                    "crypto": r.flags.crypto,
                    "bcov": r.flags.bcov,
                    "rota_a": r.flags.rota_a,
                }
                for r in self.records
            ],
        }


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _burden(flags: BovEntericFlags) -> int:
    """Count the number of positive targets in a single panel result."""
    return (
        int(flags.etec_k99)
        + int(flags.salmonella)
        + int(flags.crypto)
        + int(flags.bcov)
        + int(flags.rota_a)
    )


def _is_truthy(val: str) -> bool:
    """Parse a boolean cell from a CSV row (0/1, true/false, yes/no)."""
    return val.strip().lower() in ("1", "true", "yes")


def _compute_direction(
    records: tuple[BovEntericRecord, ...],
) -> tuple[BovEntericTrendDirection, bool, bool]:
    """Compute trend direction and ETEC event flags.

    ETEC K99 is the dominant signal because it is the leading cause of
    near-certain neonatal mortality and the primary antibiotic-urgency
    trigger.  A new ETEC detection always yields EMERGING regardless of
    overall burden change; ETEC clearance with a reduced overall burden
    yields RESOLVING.

    Overall burden (positive target count) is compared between the first
    and second half of the time series to classify non-ETEC-event trends.
    A mean change of >=0.5 targets per interval is required to classify a
    directional trend.

    Args:
        records: At least two :class:`BovEntericRecord` in chronological order.

    Returns:
        Tuple of (direction, etec_newly_detected, etec_cleared).
    """
    n = len(records)
    last_flags = records[-1].flags

    etec_in_prior = any(r.flags.etec_k99 for r in records[:-1])
    etec_newly_detected = last_flags.etec_k99 and not etec_in_prior
    etec_cleared = etec_in_prior and not last_flags.etec_k99

    if etec_newly_detected:
        return BovEntericTrendDirection.EMERGING, True, False

    burdens = [_burden(r.flags) for r in records]

    if all(b == 0 for b in burdens):
        return BovEntericTrendDirection.STABLE_CLEAR, False, False

    if burdens[-1] == 0 and burdens[0] > 0:
        return BovEntericTrendDirection.RESOLVING, False, etec_cleared

    mid = max(1, n // 2)
    early_mean = sum(burdens[:mid]) / mid
    recent_mean = sum(burdens[mid:]) / (n - mid)

    if recent_mean - early_mean >= 0.5:
        return BovEntericTrendDirection.EMERGING, False, etec_cleared
    if early_mean - recent_mean >= 0.5:
        return BovEntericTrendDirection.RESOLVING, False, etec_cleared

    return BovEntericTrendDirection.STABLE_ENDEMIC, False, etec_cleared


def _build_trend_text(
    direction: BovEntericTrendDirection,
    n_samples: int,
    etec_newly_detected: bool,
    etec_cleared: bool,
    salmonella_persistent: bool,
    zoonotic_risk_count: int,
    antibiotic_required_count: int,
    worst_level: BovEntericAlertLevel,
) -> tuple[str, str]:
    """Build (interpretation, recommended_action) strings for an enteric trend.

    Args:
        direction: Pre-computed trend direction.
        n_samples: Total number of samples in the time series.
        etec_newly_detected: True if ETEC K99 first appeared in the last sample.
        etec_cleared: True if ETEC K99 was positive before but is now negative.
        salmonella_persistent: True if Salmonella positive in all intervals.
        zoonotic_risk_count: Number of intervals with Salmonella or Crypto.
        antibiotic_required_count: Number of intervals with ETEC or Salmonella.
        worst_level: Highest alert level across all records.

    Returns:
        Tuple of (interpretation, recommended_action) as non-empty strings.
    """
    if direction is BovEntericTrendDirection.INSUFFICIENT_DATA:
        return (
            "Insufficient data for trend analysis (only 1 sample collected). "
            "At least 2 consecutive samples are required to determine the "
            "neonatal scours trajectory.",
            "Collect the next scheduled monitoring sample and re-run "
            "bov-enteric-trend after the second sample is available. "
            "If neonatal diarrhea is observed clinically, act on the current "
            "single-sample result without waiting for trend data.",
        )

    salmonella_note = (
        " Salmonella detected in every monitoring interval: persistent herd-level "
        "contamination or a chronic shedder is suspected."
        if salmonella_persistent
        else ""
    )
    zoonotic_note = (
        f" Zoonotic risk active in {zoonotic_risk_count} interval(s) "
        "(Salmonella and/or Cryptosporidium): enforce strict hygiene with gloves "
        "and handwashing for all personnel handling affected calves."
        if zoonotic_risk_count > 0
        else ""
    )
    abx_note = (
        f" Antibiotics indicated in {antibiotic_required_count} interval(s) "
        "(ETEC K99 and/or Salmonella detected)."
        if antibiotic_required_count > 0
        else ""
    )

    if direction is BovEntericTrendDirection.EMERGING:
        if etec_newly_detected:
            interp = (
                "ETEC K99 (faeG) newly detected in the most recent sample after "
                f"{n_samples - 1} consecutive negative interval(s). "
                "ETEC K99 causes near-certain mortality (>90%) in calves <4 days old "
                "without IV antibiotic therapy and fluid replacement within 24-48 h."
                + salmonella_note
                + zoonotic_note
            )
            action = (
                "CRITICAL: Start IV or oral electrolyte rehydration immediately. "
                "Administer IV antibiotics (ampicillin or amoxicillin-clavulanate) "
                "within 2 h of diagnosis -- do NOT wait for sensitivity results. "
                "Isolate affected calves. Sanitise pens. Review colostrum management "
                "and dam vaccination (K99 fimbrial vaccine programme). Notify the "
                "farm veterinarian. Increase monitoring to daily for the affected group."
            )
        else:
            interp = (
                f"Enteric pathogen burden is increasing across the {n_samples}-sample "
                f"time series (worst alert: {worst_level.value}). "
                "Scours outbreak risk is rising in the calf group."
                + salmonella_note
                + zoonotic_note
                + abx_note
            )
            action = (
                "Increase monitoring frequency. Review colostrum management, pen hygiene, "
                "and calf nutrition. If ETEC K99 or Salmonella is active, initiate "
                "antibiotic therapy and isolate affected calves. Consult the farm "
                "veterinarian for treatment protocols."
            )
        return interp, action

    if direction is BovEntericTrendDirection.RESOLVING:
        if etec_cleared:
            interp = (
                "ETEC K99 (faeG) was detected in earlier samples but is negative "
                f"in the most recent sample across {n_samples} intervals. "
                "Current treatment and hygiene measures appear effective."
                + salmonella_note
                + zoonotic_note
            )
            action = (
                "Maintain current antibiotic, rehydration, and hygiene protocols. "
                "Continue weekly monitoring until at least three consecutive ETEC-negative "
                "results are confirmed before relaxing pen-sanitation procedures. "
                "Environmental ETEC contamination persists in bedding and water troughs "
                "after clinical resolution; thorough disinfection of all calf-contact "
                "surfaces is essential."
            )
        else:
            interp = (
                f"Enteric pathogen burden is decreasing across the {n_samples}-sample "
                "time series. Control measures appear to be working."
                + salmonella_note
                + zoonotic_note
                + abx_note
            )
            action = (
                "Maintain current management, treatment, and hygiene programme. "
                "Continue surveillance until at least three consecutive negative "
                "results are achieved before relaxing controls."
            )
        return interp, action

    if direction is BovEntericTrendDirection.STABLE_CLEAR:
        interp = (
            f"All enteric targets negative across {n_samples} monitoring interval(s). "
            "No active neonatal calf diarrhea pathogens detected."
        )
        action = (
            "Maintain the routine surveillance schedule. Ensure sample-adequacy "
            "controls (e.g. 16S rRNA) are positive to confirm extraction quality. "
            "Review colostrum programme and dam vaccination before peak calving season."
        )
        return interp, action

    # STABLE_ENDEMIC
    if worst_level in (BovEntericAlertLevel.CRITICAL, BovEntericAlertLevel.HIGH):
        interp = (
            f"High-alert enteric pathogens detected across multiple of the {n_samples} "
            f"monitoring intervals without a clear downward trend "
            f"(worst alert: {worst_level.value})." + salmonella_note + zoonotic_note + abx_note
        )
        action = (
            "Reassess the calf-management programme with the farm veterinarian. "
            "Submit faecal samples to a diagnostic laboratory for culture and "
            "sensitivity testing. Review antimicrobial therapy selection. "
            "Implement pen-all-in all-out protocols to break the pathogen cycle. "
            "Audit colostrum programme and maternal vaccination coverage."
        )
    else:
        interp = (
            f"Enteric pathogens detected in one or more of the {n_samples} monitoring "
            f"intervals without a clear trend (worst alert: {worst_level.value}). "
            "Endemic low-level scours present in the calf group."
            + salmonella_note
            + zoonotic_note
            + abx_note
        )
        action = (
            "Maintain current biosecurity, nutrition, and hygiene protocols. "
            "Consult the farm veterinarian if clinical scours is observed. "
            "Review colostrum quality and dam vaccination programme. "
            "Consider diagnostic culture if the same pathogen persists across "
            "multiple intervals."
        )
    return interp, action


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def analyse_bov_enteric_trend(records: Sequence[BovEntericRecord]) -> BovEntericTrend:
    """Compute a neonatal calf diarrhea trajectory from a chronological sequence.

    The input sequence must be ordered chronologically (oldest first).
    Two or more records are required for trend analysis; a single record
    returns direction=INSUFFICIENT_DATA.

    Args:
        records: Chronological sequence of :class:`BovEntericRecord`
            (oldest first).  Must not be empty.

    Raises:
        ValueError: If ``records`` is empty.

    Returns:
        :class:`BovEntericTrend` with trajectory direction, ETEC event flags,
        Salmonella persistence, zoonotic risk count, and a field-actionable
        recommendation.

    Example::

        from lamp_forge.bov_enteric_risk import BovEntericFlags
        from lamp_forge.bov_enteric_trend import (
            BovEntericRecord,
            analyse_bov_enteric_trend,
        )

        recs = [
            BovEntericRecord("CG1-wk1",
                BovEntericFlags(etec_k99=False, salmonella=False,
                                crypto=False, bcov=False, rota_a=True)),
            BovEntericRecord("CG1-wk2",
                BovEntericFlags(etec_k99=True,  salmonella=False,
                                crypto=False, bcov=False, rota_a=True)),
        ]
        trend = analyse_bov_enteric_trend(recs)
        print(trend.direction)           # BovEntericTrendDirection.EMERGING
        print(trend.etec_newly_detected) # True
    """
    if not records:
        raise ValueError("records must contain at least one BovEntericRecord")

    recs = tuple(records)
    n = len(recs)

    worst: BovEntericAlertLevel = BovEntericAlertLevel.NEGATIVE
    for r in recs:
        lvl = assess_bov_enteric_risk(r.flags).alert_level
        if _LEVEL_ORDER.index(lvl) > _LEVEL_ORDER.index(worst):
            worst = lvl

    salmonella_persistent = all(r.flags.salmonella for r in recs)
    zoonotic_risk_count = sum(1 for r in recs if r.flags.salmonella or r.flags.crypto)
    antibiotic_required_count = sum(1 for r in recs if r.flags.etec_k99 or r.flags.salmonella)

    if n < 2:
        interp, action = _build_trend_text(
            BovEntericTrendDirection.INSUFFICIENT_DATA,
            n,
            False,
            False,
            salmonella_persistent,
            zoonotic_risk_count,
            antibiotic_required_count,
            worst,
        )
        return BovEntericTrend(
            n_samples=n,
            direction=BovEntericTrendDirection.INSUFFICIENT_DATA,
            etec_newly_detected=False,
            etec_cleared=False,
            salmonella_persistent=salmonella_persistent,
            zoonotic_risk_count=zoonotic_risk_count,
            antibiotic_required_count=antibiotic_required_count,
            worst_alert_level=worst,
            interpretation=interp,
            recommended_action=action,
            records=recs,
        )

    direction, etec_newly_detected, etec_cleared = _compute_direction(recs)
    interp, action = _build_trend_text(
        direction,
        n,
        etec_newly_detected,
        etec_cleared,
        salmonella_persistent,
        zoonotic_risk_count,
        antibiotic_required_count,
        worst,
    )

    return BovEntericTrend(
        n_samples=n,
        direction=direction,
        etec_newly_detected=etec_newly_detected,
        etec_cleared=etec_cleared,
        salmonella_persistent=salmonella_persistent,
        zoonotic_risk_count=zoonotic_risk_count,
        antibiotic_required_count=antibiotic_required_count,
        worst_alert_level=worst,
        interpretation=interp,
        recommended_action=action,
        records=recs,
    )


def records_from_csv(path: Path) -> list[BovEntericRecord]:
    """Read a monitoring-spreadsheet CSV into a list of :class:`BovEntericRecord`.

    Expected columns (order-independent, header required):

        ``sample_id``   -- short identifier for the interval (required).
        ``date``        -- ISO date string (optional; leave blank or omit column).
        ``etec_k99``, ``salmonella``, ``crypto``, ``bcov``, ``rota_a``
                        -- bool (0/1 or true/false).

    Rows are returned in file order (oldest first by convention).

    Args:
        path: Path to the CSV file.

    Raises:
        ValueError: If a required column is missing or a row cannot be parsed.

    Returns:
        List of :class:`BovEntericRecord` in file order.
    """
    required = {"sample_id", "etec_k99", "salmonella", "crypto", "bcov", "rota_a"}
    result: list[BovEntericRecord] = []
    with path.open(encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        if reader.fieldnames is None:
            raise ValueError(f"Empty CSV file: {path}")
        col_names = {f.strip().lower() for f in reader.fieldnames}
        missing = required - col_names
        if missing:
            raise ValueError(f"CSV missing required column(s): {', '.join(sorted(missing))}")
        for i, row in enumerate(reader, start=2):
            try:
                sample_id = row["sample_id"].strip()
                if not sample_id:
                    raise ValueError("sample_id must not be empty")
                date_raw = row.get("date", "").strip() or None
                flags = BovEntericFlags(
                    etec_k99=_is_truthy(row["etec_k99"]),
                    salmonella=_is_truthy(row["salmonella"]),
                    crypto=_is_truthy(row["crypto"]),
                    bcov=_is_truthy(row["bcov"]),
                    rota_a=_is_truthy(row["rota_a"]),
                )
            except KeyError as exc:
                raise ValueError(f"Row {i}: missing column {exc}") from exc
            result.append(BovEntericRecord(sample_id=sample_id, flags=flags, date=date_raw))
    return result


def records_from_enteric_json(paths: Sequence[Path]) -> list[BovEntericRecord]:
    """Load a sequence of bov-enteric-risk JSON outputs as :class:`BovEntericRecord` list.

    Each JSON file (produced by ``lamp-forge bov-enteric-risk --out-json``) supplies
    the ``panel_flags`` dict for one sample.  The ``sample_id`` defaults to the
    file stem; an optional top-level ``"date"`` key is used if present.

    Files are added in the order given; callers must supply them oldest-first.

    Args:
        paths: Sequence of paths to JSON files (oldest first).

    Raises:
        ValueError: If a JSON file is missing the ``panel_flags`` key.

    Returns:
        List of :class:`BovEntericRecord` in the supplied order.
    """
    result: list[BovEntericRecord] = []
    for p in paths:
        with p.open(encoding="utf-8") as fh:
            data: dict[str, object] = json.load(fh)
        pf = data.get("panel_flags")
        if not isinstance(pf, dict):
            raise ValueError(
                f"JSON file {p} is missing a 'panel_flags' dict. "
                "Generate it with: lamp-forge bov-enteric-risk --out-json <path>"
            )
        flags = BovEntericFlags(
            etec_k99=bool(pf.get("etec_k99", False)),
            salmonella=bool(pf.get("salmonella", False)),
            crypto=bool(pf.get("crypto", False)),
            bcov=bool(pf.get("bcov", False)),
            rota_a=bool(pf.get("rota_a", False)),
        )
        date = str(data["date"]) if "date" in data else None
        result.append(BovEntericRecord(sample_id=p.stem, flags=flags, date=date))
    return result


def write_trend_json(trend: BovEntericTrend, path: Path) -> None:
    """Write the trend assessment to a JSON file.

    Args:
        trend: :class:`BovEntericTrend` to serialise.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        json.dump(trend.to_dict(), fh, indent=2)
        fh.write("\n")


def write_trend_csv(trend: BovEntericTrend, path: Path) -> None:
    """Write the chronological timeline of pathogen flags to a CSV file.

    Each row represents one monitoring interval with per-pathogen boolean
    columns and the per-sample alert level from
    :func:`~lamp_forge.bov_enteric_risk.assess_bov_enteric_risk`.

    Args:
        trend: :class:`BovEntericTrend` to export.
        path: Output path.  Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "sample_id",
                "date",
                "etec_k99",
                "salmonella",
                "crypto",
                "bcov",
                "rota_a",
                "alert_level",
            ]
        )
        for r in trend.records:
            lvl = assess_bov_enteric_risk(r.flags).alert_level.value
            writer.writerow(
                [
                    r.sample_id,
                    r.date or "",
                    int(r.flags.etec_k99),
                    int(r.flags.salmonella),
                    int(r.flags.crypto),
                    int(r.flags.bcov),
                    int(r.flags.rota_a),
                    lvl,
                ]
            )

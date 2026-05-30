"""Time-to-positive (TTP) estimator for LAMP / RT-LAMP reactions.

Portable LAMP platforms (e.g. BioVind BioID, 30-60 min read-out) need to
know whether an assay will turn positive before the device time window
closes.  This module uses the empirically validated linear relationship
between TTP and log10(initial copy count) to predict the expected read-out
time across a range of sample loads.

Model
-----
::

    TTP(N) = ttp_one_copy_min - slope_min_per_decade * log10(max(N, 1.0))
    TTP(N) = max(min_ttp_min, TTP(N))

where:

* ``ttp_one_copy_min`` (min): threshold time at a mean of 1 copy per reaction.
  Based on Poisson single-molecule experiments with Bst 2.0 WarmStart at 65 degC.
* ``slope_min_per_decade`` (min/decade): TTP acceleration per 10-fold increase
  in initial copy count (template abundance effect).
* ``min_ttp_min`` (min): physical floor -- LAMP cannot amplify faster than
  protocol minimum incubation time regardless of copy count.

Default parameters (Bst 2.0 WarmStart, 65 degC, 25 uL reaction, LAMP buffer B;
NEB WarmStart LAMP Kit E1700 / Tanner 2012)::

    ttp_one_copy_min         = 55 min   (conservative 95th-percentile single-copy TTP)
    slope_min_per_decade     = 6.0 min/decade
    min_ttp_min              = 10 min   (protocol-mandated minimum incubation)
    device_window_min        = 60 min   (BioVind BioID 30-60 min run window)

At these defaults:

    1 copy       -> 55 min  (borderline in 60-min window)
    10 copies    -> 49 min  (inside window)
    100 copies   -> 43 min
    1 000 copies -> 37 min
    10 000 copies -> 31 min
    100 000 copies -> 25 min

One-step RT-LAMP (RTx + Bst 2.0) typically adds 3-8 min compared to
pure DNA-LAMP due to the reverse-transcription lag; use the RNA preset
(``TtpPreset.RT_LAMP``) for RNA-target assays (PRRSV, FMDV, avian influenza).

Integration with LOD estimation
--------------------------------
:func:`ttp_from_sample_concentration` combines this model with the extraction
chain from :mod:`lamp_forge.lod` so users can answer: "If my produced-water
sample contains X GE/mL of SRB, will my assay turn positive within 60 min?"

References:
    Tanner NA et al. (2012) BioTechniques 53(2):81-89.
    Notomi T et al. (2000) Nucleic Acids Res 28(12):e63.
    Dao Thi VL et al. (2020) Science 370(6518):914-917.
    NEB WarmStart LAMP Kit manual (E1700).
    NEB WarmStart RTx Reverse Transcriptase protocol.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from lamp_forge.lod import ExtractionParams


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------


class TtpPreset(StrEnum):
    """Named parameter presets for common LAMP chemistry / instrument combinations.

    Attributes:
        DNA_LAMP: Pure DNA-LAMP with Bst 2.0 WarmStart at 65 degC. Default for
            bacterial (SRB, MCR, IRB, APB, NRB) and DNA-virus (ASFV) assays.
        RT_LAMP: One-step RT-LAMP (NEB RTx + Bst 2.0 WarmStart, 63-65 degC).
            Adds ~5 min TTP at low copies compared to DNA-LAMP due to reverse-
            transcription lag. Use for PRRSV, FMDV, avian influenza, NDV.
        FAST_LAMP: Optimised DNA-LAMP with elevated enzyme loading and early
            fluorescence detection at 65 degC. Roughly 5 min faster than
            DNA_LAMP at every copy count. Useful for sensitivity benchmarking.
    """

    DNA_LAMP = "dna-lamp"
    RT_LAMP = "rt-lamp"
    FAST_LAMP = "fast-lamp"


# Keyed by TtpPreset: (ttp_one_copy_min, slope_min_per_decade, min_ttp_min)
_PRESET_PARAMS: dict[TtpPreset, tuple[float, float, float]] = {
    TtpPreset.DNA_LAMP: (55.0, 6.0, 10.0),
    TtpPreset.RT_LAMP: (60.0, 6.0, 12.0),
    TtpPreset.FAST_LAMP: (50.0, 6.0, 8.0),
}


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class TtpParams:
    """Parameters for the linear TTP-vs-log10(copies) model.

    Attributes:
        ttp_one_copy_min: Expected TTP (minutes) at a mean of 1 copy per
            reaction; the conservative 95th-percentile threshold from Poisson
            single-molecule LAMP experiments. Default 55 min (Bst 2.0, 65 degC).
        slope_min_per_decade: TTP reduction (minutes) per 10-fold increase in
            initial copy count. Default 6 min/decade.
        device_window_min: Maximum TTP for a result to be reportable within the
            device run time. Default 60 min (BioVind BioID 30-60 min window).
        min_ttp_min: Physical lower bound on TTP (minutes). The reaction cannot
            amplify faster than the protocol minimum incubation time. Default
            10 min for DNA-LAMP; 12 min for RT-LAMP.
    """

    ttp_one_copy_min: float = 55.0
    slope_min_per_decade: float = 6.0
    device_window_min: float = 60.0
    min_ttp_min: float = 10.0

    def __post_init__(self) -> None:  # noqa: D105
        if self.ttp_one_copy_min <= 0.0:
            raise ValueError(f"ttp_one_copy_min must be positive, got {self.ttp_one_copy_min!r}")
        if self.slope_min_per_decade <= 0.0:
            raise ValueError(
                f"slope_min_per_decade must be positive, got {self.slope_min_per_decade!r}"
            )
        if self.device_window_min <= 0.0:
            raise ValueError(f"device_window_min must be positive, got {self.device_window_min!r}")
        if self.min_ttp_min < 0.0:
            raise ValueError(f"min_ttp_min must be >= 0, got {self.min_ttp_min!r}")

    @classmethod
    def from_preset(cls, preset: TtpPreset, device_window_min: float = 60.0) -> TtpParams:
        """Build TtpParams from a named chemistry preset.

        Args:
            preset: Named chemistry / instrument combination.
            device_window_min: Device run time in minutes (default 60).

        Returns:
            TtpParams populated from the preset.
        """
        ttp_one, slope, min_ttp = _PRESET_PARAMS[preset]
        return cls(
            ttp_one_copy_min=ttp_one,
            slope_min_per_decade=slope,
            device_window_min=device_window_min,
            min_ttp_min=min_ttp,
        )


@dataclass(frozen=True, slots=True)
class TtpPoint:
    """TTP prediction at a single initial copy count.

    Attributes:
        copies_per_reaction: Mean initial target copies per reaction.
        ttp_minutes: Predicted time-to-positive (minutes).
        in_device_window: True if ttp_minutes <= TtpParams.device_window_min.
    """

    copies_per_reaction: float
    ttp_minutes: float
    in_device_window: bool


@dataclass(frozen=True, slots=True)
class TtpTable:
    """TTP predictions across a log-spaced range of copy counts.

    Attributes:
        params: TtpParams used to build the table.
        points: TtpPoint list, ordered ascending by copies_per_reaction.
        min_detectable_copies: Smallest copies_per_reaction where
            in_device_window is True. None if no point in the table passes.
    """

    params: TtpParams
    points: list[TtpPoint]
    min_detectable_copies: float | None


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def ttp_at_copies(copies_per_reaction: float, params: TtpParams) -> TtpPoint:
    """Predict TTP for a given mean copy count per reaction.

    Applies the linear log10 model::

        TTP = ttp_one_copy_min - slope * log10(max(copies, 1.0))
        TTP = max(min_ttp_min, TTP)

    Args:
        copies_per_reaction: Mean target copies added to the reaction (> 0).
        params: TTP model parameters.

    Returns:
        TtpPoint with predicted TTP and window-pass flag.

    Raises:
        ValueError: If copies_per_reaction <= 0.
    """
    if copies_per_reaction <= 0.0:
        raise ValueError(f"copies_per_reaction must be positive, got {copies_per_reaction!r}")
    log10_copies = math.log10(max(copies_per_reaction, 1.0))
    raw_ttp = params.ttp_one_copy_min - params.slope_min_per_decade * log10_copies
    ttp = max(params.min_ttp_min, raw_ttp)
    return TtpPoint(
        copies_per_reaction=copies_per_reaction,
        ttp_minutes=ttp,
        in_device_window=ttp <= params.device_window_min,
    )


def ttp_table(
    params: TtpParams,
    copies_range: tuple[float, float] = (1.0, 1e6),
    n_points: int = 13,
) -> TtpTable:
    """Build a TTP table across a log-spaced range of copy counts.

    Default of 13 points over 1-10^6 gives a half-decade step (10^0, 10^0.5,
    10^1, ...) that is conventionally reported in molecular diagnostics papers.

    Args:
        params: TTP model parameters.
        copies_range: (min_copies, max_copies) inclusive. Both must be > 0 and
            min_copies <= max_copies.
        n_points: Number of log-spaced evaluation points. Must be >= 2.

    Returns:
        TtpTable with predictions and the minimum detectable copy count.

    Raises:
        ValueError: If copies_range is invalid or n_points < 2.
    """
    lo, hi = copies_range
    if lo <= 0.0 or hi <= 0.0 or lo > hi:
        raise ValueError(f"copies_range must satisfy 0 < lo <= hi, got {copies_range!r}")
    if n_points < 2:
        raise ValueError(f"n_points must be >= 2, got {n_points!r}")

    log_lo = math.log10(lo)
    log_hi = math.log10(hi)
    step = (log_hi - log_lo) / (n_points - 1) if n_points > 1 else 0.0

    points: list[TtpPoint] = [
        ttp_at_copies(10.0 ** (log_lo + i * step), params) for i in range(n_points)
    ]

    min_det: float | None = next(
        (pt.copies_per_reaction for pt in points if pt.in_device_window), None
    )
    return TtpTable(params=params, points=points, min_detectable_copies=min_det)


def ttp_from_sample_concentration(
    conc_copies_per_ml: float,
    extraction_params: ExtractionParams,
    ttp_params: TtpParams,
) -> TtpPoint:
    """Predict TTP for a sample at a known genome-equivalent concentration.

    Bridges the extraction chain (sample -> eluate -> reaction) to the TTP
    model.  Answers: "If my produced-water sample has X GE/mL of SRB, will
    the BioVind assay turn positive within the 60-min device window?"

    Args:
        conc_copies_per_ml: Target GE concentration in the original sample
            (copies/mL). Must be > 0.
        extraction_params: Sample preparation chain parameters (volumes +
            efficiency); see :class:`lamp_forge.lod.ExtractionParams`.
        ttp_params: TTP model parameters.

    Returns:
        TtpPoint with predicted TTP and window-pass flag.

    Raises:
        ValueError: If conc_copies_per_ml <= 0.
    """
    if conc_copies_per_ml <= 0.0:
        raise ValueError(f"conc_copies_per_ml must be positive, got {conc_copies_per_ml!r}")
    copies_per_ul = conc_copies_per_ml / 1000.0
    copies_in_rxn = copies_per_ul * extraction_params.copies_per_rxn_per_copy_per_ul
    return ttp_at_copies(copies_in_rxn, ttp_params)


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_ttp_csv(table: TtpTable, path: Path) -> None:
    """Write a TTP table to a CSV file.

    Args:
        table: TtpTable to serialise.
        path: Output file path. Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "copies_per_reaction",
                "ttp_minutes",
                "in_device_window",
                "device_window_min",
                "ttp_one_copy_min",
                "slope_min_per_decade",
                "min_ttp_min",
            ]
        )
        for pt in table.points:
            writer.writerow(
                [
                    f"{pt.copies_per_reaction:.2f}",
                    f"{pt.ttp_minutes:.1f}",
                    str(pt.in_device_window).lower(),
                    f"{table.params.device_window_min:.1f}",
                    f"{table.params.ttp_one_copy_min:.1f}",
                    f"{table.params.slope_min_per_decade:.2f}",
                    f"{table.params.min_ttp_min:.1f}",
                ]
            )

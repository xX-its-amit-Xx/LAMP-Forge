"""Assay limit-of-detection (LOD) and copies-per-reaction estimator.

Portable LAMP diagnostics live at the bottom of what nucleic-acid
amplification can detect. This module quantifies that limit given the
sample-to-reaction chain: sample -> extraction -> eluate -> LAMP reaction.

The underlying model is the **Poisson single-molecule sampling** law:

    P(>= 1 copy | mean copies = lambda) = 1 - exp(-lambda)

Inverting gives the critical mean copy count per reaction at a given
detection probability::

    lambda_LOD = -ln(1 - P_target)

    LOD_95 => lambda ~ 2.996 copies / reaction
    LOD_99 => lambda ~ 4.605 copies / reaction

The *back-calculated* LOD in the original sample follows from the
extraction chain:

    copies_in_rxn = C_sample [copies/uL]
                    * V_sample [uL]
                    * eta_extraction
                    * (V_rxn_input [uL] / V_eluate [uL])

Solving for C_sample at lambda_LOD gives the reportable LOD in copies/uL
(or copies/mL) of the original sample.

References:
    Dingle & Basler (2012) Poisson statistics and end-point digital LAMP.
    Bio-Rad dPCR application notes on Poisson correction.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ExtractionParams:
    """Physical parameters of the sample-preparation-to-reaction chain.

    All volumes are in microlitres (uL). Typical field / clinical values:

    * sample_volume_ul: 200 (200 uL blood, nasal-swab eluent, or produced water)
    * extraction_efficiency: 0.50 (50% is a conservative column-kit estimate)
    * eluate_volume_ul: 50-100 (column kits typically elute in 50-100 uL buffer)
    * reaction_input_ul: 2-5 (fraction of eluate added to the 25 uL LAMP reaction)

    Attributes:
        sample_volume_ul: uL of original sample input to the extraction step.
        extraction_efficiency: Fraction of target copies recovered (0 < eta <= 1).
        eluate_volume_ul: Final extraction eluate volume in uL.
        reaction_input_ul: uL of eluate added to the LAMP reaction.
    """

    sample_volume_ul: float
    extraction_efficiency: float
    eluate_volume_ul: float
    reaction_input_ul: float

    def __post_init__(self) -> None:  # noqa: D105
        if self.sample_volume_ul <= 0.0:
            raise ValueError(f"sample_volume_ul must be positive, got {self.sample_volume_ul!r}")
        if not (0.0 < self.extraction_efficiency <= 1.0):
            raise ValueError(
                f"extraction_efficiency must be in (0, 1], got {self.extraction_efficiency!r}"
            )
        if self.eluate_volume_ul <= 0.0:
            raise ValueError(f"eluate_volume_ul must be positive, got {self.eluate_volume_ul!r}")
        if self.reaction_input_ul <= 0.0:
            raise ValueError(f"reaction_input_ul must be positive, got {self.reaction_input_ul!r}")
        if self.reaction_input_ul > self.eluate_volume_ul:
            raise ValueError(
                f"reaction_input_ul ({self.reaction_input_ul}) must not exceed "
                f"eluate_volume_ul ({self.eluate_volume_ul})"
            )

    @property
    def copies_per_rxn_per_copy_per_ul(self) -> float:
        """Conversion factor: copies in reaction per copy/uL in original sample.

        Equals V_sample * eta * (V_rxn_input / V_eluate), with units of uL.
        This is the *effective sample volume* represented in one reaction after
        all extraction and dilution losses.
        """
        return (
            self.sample_volume_ul
            * self.extraction_efficiency
            * self.reaction_input_ul
            / self.eluate_volume_ul
        )


@dataclass(frozen=True, slots=True)
class LodEstimate:
    """LOD summary for one extraction chain at one detection probability.

    Attributes:
        extraction: The extraction-chain parameters used.
        detection_probability: Requested P(detect | target present at LOD).
        lod_copies_per_reaction: Mean copies per reaction at the LOD (lambda_LOD).
        lod_copies_per_ul: Back-calculated LOD in the original sample (copies/uL).
        lod_copies_per_ml: Same converted to copies/mL (the conventional clinical
            and industrial unit for liquid samples).
        lod_genome_eq_per_ml: Alias of lod_copies_per_ml; "genome equivalents
            per mL" (GE/mL) is preferred nomenclature in molecular diagnostics.
        copies_per_cell: Gene copy number per organism cell (default 1 for
            single-copy genes).  Use 4-10 for 16S rRNA (multi-copy ribosomal
            locus) or the plasmid copy number for resistance-gene targets.
        lod_cells_per_ml: LOD expressed as cells (organisms) per mL of original
            sample, equal to ``lod_copies_per_ml / copies_per_cell``.  For
            single-copy genes this equals ``lod_copies_per_ml``.
    """

    extraction: ExtractionParams
    detection_probability: float
    lod_copies_per_reaction: float
    lod_copies_per_ul: float
    lod_copies_per_ml: float
    lod_genome_eq_per_ml: float
    copies_per_cell: int
    lod_cells_per_ml: float


# ---------------------------------------------------------------------------
# Core functions
# ---------------------------------------------------------------------------


def poisson_lod(detection_probability: float) -> float:
    """Mean copies per reaction at which P(at least one copy) = detection_probability.

    Uses the exact Poisson inverse: lambda = -ln(1 - P).

    Args:
        detection_probability: Required detection probability, e.g. 0.95 for
            LOD_95. Must be in the open interval (0, 1).

    Returns:
        lambda_LOD: mean copies per reaction at the LOD.

    Raises:
        ValueError: If ``detection_probability`` is not in (0, 1).
    """
    if not (0.0 < detection_probability < 1.0):
        raise ValueError(f"detection_probability must be in (0, 1), got {detection_probability!r}")
    return -math.log(1.0 - detection_probability)


def detection_probability_at(mean_copies: float) -> float:
    """P(at least one copy detected) given a mean copy count per reaction.

    Args:
        mean_copies: Mean number of target copies per reaction (lambda >= 0).

    Returns:
        Detection probability in [0, 1).

    Raises:
        ValueError: If ``mean_copies`` is negative.
    """
    if mean_copies < 0.0:
        raise ValueError(f"mean_copies must be non-negative, got {mean_copies!r}")
    return 1.0 - math.exp(-mean_copies)


def estimate_lod(
    extraction: ExtractionParams,
    detection_probability: float = 0.95,
    copies_per_cell: int = 1,
) -> LodEstimate:
    """Estimate the assay LOD for a given extraction chain and target probability.

    Args:
        extraction: Parameters describing the sample preparation chain.
        detection_probability: Required P(detect) at the LOD. Default 0.95.
        copies_per_cell: Number of target-gene copies per organism cell.
            Default 1 for single-copy genes (most functional markers: dsrB,
            mcrA, narG, omcA, fthfs, rpoB).  Use 4–10 for the 16S rRNA
            control channel, or the known plasmid copy number for resistance
            genes (e.g. blaKPC on a multi-copy plasmid).  When > 1, the
            returned ``lod_cells_per_ml`` is lower than ``lod_copies_per_ml``
            by the same factor, reflecting the analytical advantage of
            targeting a multi-copy locus.

    Returns:
        :class:`LodEstimate` with lambda_LOD and back-calculated sample LOD.

    Raises:
        ValueError: If ``detection_probability`` is not in (0, 1), or
            ``copies_per_cell`` is less than 1.
    """
    if copies_per_cell < 1:
        raise ValueError(f"copies_per_cell must be >= 1, got {copies_per_cell!r}")
    lod_rxn = poisson_lod(detection_probability)
    conv = extraction.copies_per_rxn_per_copy_per_ul
    lod_per_ul = lod_rxn / conv
    lod_per_ml = lod_per_ul * 1000.0
    lod_cells = lod_per_ml / copies_per_cell
    return LodEstimate(
        extraction=extraction,
        detection_probability=detection_probability,
        lod_copies_per_reaction=lod_rxn,
        lod_copies_per_ul=lod_per_ul,
        lod_copies_per_ml=lod_per_ml,
        lod_genome_eq_per_ml=lod_per_ml,
        copies_per_cell=copies_per_cell,
        lod_cells_per_ml=lod_cells,
    )


def lod_table(
    extraction: ExtractionParams,
    probabilities: tuple[float, ...] = (0.90, 0.95, 0.99, 0.999),
    copies_per_cell: int = 1,
) -> list[LodEstimate]:
    """Return LOD estimates at multiple detection probabilities.

    A convenience wrapper around :func:`estimate_lod` that produces a
    table useful for comparing LOD_90 / LOD_95 / LOD_99 / LOD_99.9 —
    the set of thresholds most diagnostic guidance documents request.

    Args:
        extraction: Sample preparation chain parameters.
        probabilities: Detection probabilities to evaluate.
        copies_per_cell: Gene copy number per organism cell; forwarded to
            :func:`estimate_lod`.  Default 1 (single-copy genes).

    Returns:
        List of :class:`LodEstimate`, one per probability, in input order.
    """
    return [estimate_lod(extraction, p, copies_per_cell) for p in probabilities]


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_lod_csv(estimates: list[LodEstimate], path: Path) -> None:
    """Write an LOD table to a CSV file.

    Args:
        estimates: List of :class:`LodEstimate` to serialise.
        path: Output file path. Parent directories are created if needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "detection_probability",
                "lod_copies_per_reaction",
                "lod_copies_per_ul_sample",
                "lod_copies_per_ml_sample",
                "lod_genome_eq_per_ml",
                "copies_per_cell",
                "lod_cells_per_ml_sample",
                "sample_volume_ul",
                "extraction_efficiency",
                "eluate_volume_ul",
                "reaction_input_ul",
                "effective_sample_ul_in_rxn",
            ]
        )
        for e in estimates:
            writer.writerow(
                [
                    f"{e.detection_probability:.4f}",
                    f"{e.lod_copies_per_reaction:.4f}",
                    f"{e.lod_copies_per_ul:.4f}",
                    f"{e.lod_copies_per_ml:.2f}",
                    f"{e.lod_genome_eq_per_ml:.2f}",
                    f"{e.copies_per_cell}",
                    f"{e.lod_cells_per_ml:.2f}",
                    f"{e.extraction.sample_volume_ul:.1f}",
                    f"{e.extraction.extraction_efficiency:.3f}",
                    f"{e.extraction.eluate_volume_ul:.1f}",
                    f"{e.extraction.reaction_input_ul:.1f}",
                    f"{e.extraction.copies_per_rxn_per_copy_per_ul:.3f}",
                ]
            )

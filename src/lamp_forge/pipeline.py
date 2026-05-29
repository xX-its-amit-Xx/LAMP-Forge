"""End-to-end pipeline orchestration as a single function call.

The Snakefile wires up the same steps for distributed/cluster execution; this
module exists so the CLI (and notebook walkthrough) can run the whole thing
locally in one process without spawning snakemake.
"""

from __future__ import annotations

import logging
from pathlib import Path

from lamp_forge import align, conserve, fetch, primer_design, report, specificity
from lamp_forge.types import LampPrimerSet, PipelineConfig

logger = logging.getLogger(__name__)


def run_pipeline(config: PipelineConfig) -> list[LampPrimerSet]:
    """Execute fetch → align → conserve → design → screen → report.

    All intermediate artefacts are written under ``config.output_dir`` so that
    a failed run can be inspected. The function returns the ranked primer
    sets; the same data is also persisted as JSON/CSV/HTML.
    """
    out = config.output_dir
    out.mkdir(parents=True, exist_ok=True)
    work = out / "work"
    work.mkdir(exist_ok=True)

    # ---------- 1. fetch ----------
    logger.info("[1/6] Fetching sequences for %s", config.target_name)
    raw_fasta = work / "sequences.fasta"
    records = fetch.cached_or_fetch(raw_fasta, config)
    input_manifest = fetch.manifest(records)
    logger.info("Retrieved %d sequence(s)", len(records))

    # ---------- 2. align ----------
    logger.info("[2/6] Running MAFFT alignment")
    aligned_fasta = work / "alignment.fasta"
    alignment = align.align_records(raw_fasta, aligned_fasta)
    logger.info("Alignment width: %d bp", alignment.get_alignment_length())

    # ---------- 3. conserve ----------
    logger.info("[3/6] Computing conservation track")
    track = conserve.compute_track(alignment, window_size=config.window_size)
    regions = conserve.find_conserved_regions(
        track,
        entropy_threshold=config.entropy_threshold,
        min_region_length=config.min_region_length,
    )
    report.write_conservation_tsv(track, out / "conservation.tsv")
    if not regions:
        logger.warning(
            "No conserved regions met threshold=%.2f / min_length=%d. "
            "Try relaxing conservation.entropy_threshold or conservation.min_region_length.",
            config.entropy_threshold,
            config.min_region_length,
        )

    # ---------- 4. design ----------
    logger.info("[4/6] Designing LAMP primer sets across %d region(s)", len(regions))
    primer_sets = primer_design.design_all(regions, config)
    logger.info("Designed %d candidate primer sets total", len(primer_sets))

    # ---------- 5. specificity ----------
    n_off_targets = 0
    if primer_sets:
        off_target_files = specificity.discover_off_targets(config.off_target_dir)
        n_off_targets = len(off_target_files)
        if off_target_files:
            logger.info("[5/6] Building BLAST DB from %d off-target genome(s)", n_off_targets)
            db_path = work / "off_targets_db"
            specificity.build_blast_db(off_target_files, db_path)
            specificity.screen_all(primer_sets, db_path, config)
            specificity.write_specificity_tsv(primer_sets, out / "specificity.tsv")
        else:
            logger.warning(
                "No off-target FASTA files found in %s — skipping specificity screen "
                "(specificity_score will default to 1.0 for all sets)",
                config.off_target_dir,
            )

    # ---------- 6. report ----------
    logger.info("[6/6] Writing outputs")
    report.write_json(primer_sets, out / "primer_sets.json")
    if config.generate_csv:
        report.write_csv(primer_sets, out / "primer_sets.csv")
    if config.generate_html:
        report.render_html(
            primer_sets,
            regions,
            track,
            config,
            n_input_seqs=len(records),
            n_off_targets=n_off_targets,
            path=out / "lamp_forge_report.html",
        )
    report.write_manifest(
        config,
        input_manifest,
        n_regions=len(regions),
        n_sets=len(primer_sets),
        path=out / "run_manifest.json",
    )

    return primer_sets


def _ensure_results_dir(path: Path) -> Path:
    """Create ``path`` if missing, return it. Convenience for callers."""
    path.mkdir(parents=True, exist_ok=True)
    return path

"""arq task definitions: the async wrapper around the pipeline.

The pipeline itself is sync (subprocess calls to MAFFT/BLAST + numpy).
We run it in an executor thread to avoid blocking the worker's event loop.
Progress events are persisted to the DB between stages so the UI can poll
without needing a websocket.
"""

from __future__ import annotations

import asyncio
import os
from collections.abc import Callable
from typing import Any

from lamp_forge import align, conserve, fetch, primer_design, report, specificity
from lamp_forge.types import PipelineConfig
from lamp_forge.web.config import get_settings
from lamp_forge.web.db import (
    JobStatus,
    append_progress,
    finalize_metrics,
    get_job,
    get_session_factory,
    update_status,
)
from lamp_forge.web.logging_setup import get_logger

logger = get_logger("lamp_forge.web.tasks")


def _to_pipeline_config(job_id: str, raw: dict[str, Any]) -> PipelineConfig:
    """Translate the API JobCreate schema into a :class:`PipelineConfig`.

    The schema layer prevents impossible values; this is a pure mapping —
    no further validation needed.
    """
    settings = get_settings()
    job_results = settings.results_dir / job_id
    job_input = settings.input_dir / job_id / "off_targets"
    target = raw["target"]
    primer = raw["primer"]
    conservation = raw["conservation"]
    off_targets = raw["off_targets"]
    output = raw["output"]
    amp = primer["amplicon_size"]
    return PipelineConfig(
        target_name=target["name"],
        taxon_id=target.get("taxon_id"),
        accessions=list(target.get("accessions") or []),
        gene=target.get("gene"),
        max_sequences=int(target.get("max_sequences", 20)),
        email=target["email"],
        off_target_dir=job_input,
        min_identity_threshold=float(off_targets["min_identity_threshold"]),
        min_coverage_threshold=float(off_targets["min_coverage_threshold"]),
        window_size=int(conservation["window_size"]),
        entropy_threshold=float(conservation["entropy_threshold"]),
        min_region_length=int(conservation["min_region_length"]),
        tm_min=float(primer["tm_min"]),
        tm_max=float(primer["tm_max"]),
        tm_match_tolerance=float(primer["tm_match_tolerance"]),
        gc_min=float(primer["gc_min"]),
        gc_max=float(primer["gc_max"]),
        hairpin_dg_threshold=float(primer["hairpin_dg_threshold"]),
        dimer_dg_threshold=float(primer["dimer_dg_threshold"]),
        f2_b2_min=int(amp["f2_b2_min"]),
        f2_b2_max=int(amp["f2_b2_max"]),
        output_dir=job_results,
        top_n=int(output["top_n"]),
        generate_html=bool(output["generate_html"]),
        generate_csv=bool(output["generate_csv"]),
        cache_dir=settings.cache_dir,
    )


async def _push_progress(
    job_id: str, stage: str, message: str, percent: float | None = None
) -> None:
    """Persist a progress event to the job row."""
    factory = get_session_factory()
    async with factory() as session:
        await append_progress(session, job_id, stage, message, percent)
        await session.commit()


async def _check_cancelled(job_id: str) -> bool:
    """Read the DB to check whether the user has requested cancellation."""
    factory = get_session_factory()
    async with factory() as session:
        job = await get_job(session, job_id)
        return job is not None and job.status == JobStatus.cancelled


class _CancelledError(RuntimeError):
    """Raised mid-pipeline when the user cancels — caught and turned into a clean exit."""


async def _run_blocking(func: Callable[..., Any], *args: Any, **kwargs: Any) -> Any:
    """Run a sync function in the default executor."""
    loop = asyncio.get_running_loop()
    return await loop.run_in_executor(None, lambda: func(*args, **kwargs))


# ---------------------------------------------------------------------------
# The actual task
# ---------------------------------------------------------------------------


async def run_pipeline(ctx: dict[str, Any], job_id: str) -> dict[str, Any]:
    """Arq task entry — run the LAMP-Forge pipeline for ``job_id``.

    Stages mirror :func:`lamp_forge.pipeline.run_pipeline` but interleave
    DB progress writes and cancellation checks so the UI can show a live
    timeline.
    """
    settings = get_settings()
    if settings.ncbi_api_key:
        os.environ.setdefault("NCBI_API_KEY", settings.ncbi_api_key)

    factory = get_session_factory()

    async with factory() as session:
        job = await get_job(session, job_id)
        if job is None:
            logger.error("job_not_found", job_id=job_id)
            return {"ok": False, "reason": "job_not_found"}
        cfg_raw = job.config
        await update_status(session, job_id, JobStatus.running)
        await session.commit()

    try:
        config = _to_pipeline_config(job_id, cfg_raw)
        config.output_dir.mkdir(parents=True, exist_ok=True)

        # --- 1. Fetch ---
        await _push_progress(job_id, "fetch", "Querying NCBI Entrez", 5)
        if await _check_cancelled(job_id):
            raise _CancelledError()
        records = await _run_blocking(fetch.fetch_for_config, config)
        raw_fasta = config.output_dir / "work" / "sequences.fasta"
        raw_fasta.parent.mkdir(parents=True, exist_ok=True)
        await _run_blocking(fetch.write_fasta, records, raw_fasta)
        await _push_progress(job_id, "fetch", f"Retrieved {len(records)} sequence(s)", 15)

        # --- 2. Align ---
        await _push_progress(job_id, "align", "Running MAFFT alignment", 20)
        if await _check_cancelled(job_id):
            raise _CancelledError()
        aligned = config.output_dir / "work" / "alignment.fasta"
        alignment = await _run_blocking(align.align_records, raw_fasta, aligned)
        await _push_progress(
            job_id,
            "align",
            f"Alignment width: {alignment.get_alignment_length()} bp",
            35,
        )

        # --- 3. Conserve ---
        await _push_progress(job_id, "conserve", "Computing conservation track", 40)
        if await _check_cancelled(job_id):
            raise _CancelledError()
        track = await _run_blocking(conserve.compute_track, alignment, config.window_size)
        regions = await _run_blocking(
            conserve.find_conserved_regions,
            track,
            entropy_threshold=config.entropy_threshold,
            min_region_length=config.min_region_length,
        )
        report.write_conservation_tsv(track, config.output_dir / "conservation.tsv")
        await _push_progress(job_id, "conserve", f"Detected {len(regions)} conserved region(s)", 55)

        # --- 4. Design ---
        await _push_progress(job_id, "design", "Generating candidate LAMP primer sets", 60)
        if await _check_cancelled(job_id):
            raise _CancelledError()
        primer_sets = await _run_blocking(primer_design.design_all, regions, config)
        await _push_progress(
            job_id,
            "design",
            f"Designed {len(primer_sets)} primer set(s)",
            75,
        )

        # --- 5. Specificity ---
        n_off_targets = 0
        if primer_sets:
            off_files = specificity.discover_off_targets(config.off_target_dir)
            n_off_targets = len(off_files)
            if off_files:
                await _push_progress(
                    job_id,
                    "specificity",
                    f"BLAST screening against {n_off_targets} off-target(s)",
                    80,
                )
                if await _check_cancelled(job_id):
                    raise _CancelledError()
                db_path = config.output_dir / "work" / "off_targets_db"
                await _run_blocking(specificity.build_blast_db, off_files, db_path)
                await _run_blocking(specificity.screen_all, primer_sets, db_path, config)
                specificity.write_specificity_tsv(
                    primer_sets, config.output_dir / "specificity.tsv"
                )
                await _push_progress(job_id, "specificity", "Screening complete", 90)
            else:
                await _push_progress(
                    job_id,
                    "specificity",
                    "No off-target genomes supplied; skipping screen",
                    85,
                )

        # --- 6. Report ---
        await _push_progress(job_id, "report", "Rendering report", 95)
        if await _check_cancelled(job_id):
            raise _CancelledError()
        report.write_json(primer_sets, config.output_dir / "primer_sets.json")
        if config.generate_csv:
            report.write_csv(primer_sets, config.output_dir / "primer_sets.csv")
        if config.generate_html:
            report.render_html(
                primer_sets,
                regions,
                track,
                config,
                n_input_seqs=len(records),
                n_off_targets=n_off_targets,
                path=config.output_dir / "lamp_forge_report.html",
            )
        report.write_manifest(
            config,
            fetch.manifest(records),
            n_regions=len(regions),
            n_sets=len(primer_sets),
            path=config.output_dir / "run_manifest.json",
        )

        # --- finalise ---
        async with factory() as session:
            await finalize_metrics(session, job_id, len(primer_sets), len(regions))
            await update_status(session, job_id, JobStatus.succeeded)
            await session.commit()
        await _push_progress(job_id, "complete", f"Job complete: {len(primer_sets)} sets", 100)
        logger.info(
            "job_completed",
            job_id=job_id,
            n_sets=len(primer_sets),
            n_regions=len(regions),
        )
        return {"ok": True, "n_sets": len(primer_sets), "n_regions": len(regions)}

    except _CancelledError:
        logger.info("job_cancelled_mid_run", job_id=job_id)
        # Status is already 'cancelled' (that's what _check_cancelled detected).
        return {"ok": False, "reason": "cancelled"}

    except Exception as e:
        logger.exception("job_failed", job_id=job_id, error=str(e))
        async with factory() as session:
            await update_status(
                session,
                job_id,
                JobStatus.failed,
                error=f"{type(e).__name__}: {e}",
            )
            await session.commit()
        await _push_progress(job_id, "error", f"Pipeline failed: {e}", 100)
        # Re-raise so arq records the failure in its own job state too.
        raise

"""Job submission, listing, retrieval, and cancellation."""

from __future__ import annotations

import shutil
import uuid
from pathlib import Path
from typing import Annotated

from fastapi import APIRouter, Depends, File, HTTPException, Request, UploadFile, status
from sqlalchemy.ext.asyncio import AsyncSession

from lamp_forge.web.auth import ApiKeyPrincipal, require_admin, require_api_key
from lamp_forge.web.config import Settings, get_settings
from lamp_forge.web.db import (
    Job,
    JobStatus,
    create_job,
    delete_job,
    get_job,
    list_jobs,
    session_dep,
    update_status,
)
from lamp_forge.web.logging_setup import get_logger
from lamp_forge.web.queue import enqueue_run_pipeline
from lamp_forge.web.schemas import (
    FileUploadResponse,
    JobCreate,
    JobCreateResponse,
    JobDetail,
    JobSummary,
    ProgressEvent,
)

router = APIRouter()
logger = get_logger("lamp_forge.web.routes.jobs")

MAX_UPLOAD_BYTES = 50 * 1024 * 1024  # 50 MB per file
ALLOWED_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".ffn"}


def _job_dirs(settings: Settings, job_id: str) -> tuple[Path, Path]:
    """Return (input_dir, output_dir) for a job, creating them."""
    in_dir = settings.input_dir / job_id / "off_targets"
    out_dir = settings.results_dir / job_id
    in_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)
    return in_dir, out_dir


def _to_summary(job: Job) -> JobSummary:
    return JobSummary(
        id=job.id,
        target_name=job.target_name,
        status=job.status,
        created_at=job.created_at,
        updated_at=job.updated_at,
        started_at=job.started_at,
        finished_at=job.finished_at,
        n_primer_sets=job.n_primer_sets,
        n_conserved_regions=job.n_conserved_regions,
    )


def _result_urls(settings: Settings, job: Job) -> dict[str, str]:
    """Build the public result URL map for a succeeded job."""
    if job.status != JobStatus.succeeded:
        return {}
    base = f"{settings.base_url.rstrip('/')}/api/v1/jobs/{job.id}/results"
    return {
        "html": f"{base}/lamp_forge_report.html",
        "json": f"{base}/primer_sets.json",
        "csv": f"{base}/primer_sets.csv",
        "manifest": f"{base}/run_manifest.json",
        "conservation_tsv": f"{base}/conservation.tsv",
        "specificity_tsv": f"{base}/specificity.tsv",
    }


def _to_detail(settings: Settings, job: Job) -> JobDetail:
    return JobDetail(
        id=job.id,
        target_name=job.target_name,
        status=job.status,
        created_at=job.created_at,
        updated_at=job.updated_at,
        started_at=job.started_at,
        finished_at=job.finished_at,
        n_primer_sets=job.n_primer_sets,
        n_conserved_regions=job.n_conserved_regions,
        config=job.config,
        progress=[ProgressEvent(**e) for e in (job.progress or [])],
        error=job.error,
        result_urls=_result_urls(settings, job),
    )


# ---------------------------------------------------------------------------
# Upload endpoint — clients POST off-target FASTAs here before submitting a job.
# ---------------------------------------------------------------------------


@router.post(
    "/uploads",
    response_model=FileUploadResponse,
    status_code=status.HTTP_201_CREATED,
    summary="Upload an off-target FASTA file for later use in a job",
)
async def upload(
    file: Annotated[UploadFile, File(description="A FASTA file (.fa / .fasta / .fna / .ffn)")],
    principal: ApiKeyPrincipal = Depends(require_api_key),
    settings: Settings = Depends(get_settings),
) -> FileUploadResponse:
    """Persist an uploaded FASTA file under the user's upload area.

    Returns an ``upload_id`` that can be referenced in a subsequent
    ``POST /jobs`` body's ``off_targets.fasta_files``.
    """
    if not file.filename:
        raise HTTPException(status_code=400, detail="Missing filename")
    suffix = Path(file.filename).suffix.lower()
    if suffix not in ALLOWED_FASTA_SUFFIXES:
        raise HTTPException(
            status_code=400,
            detail=f"Unsupported file type {suffix!r}. Allowed: {sorted(ALLOWED_FASTA_SUFFIXES)}",
        )

    upload_id = uuid.uuid4().hex
    user_uploads = settings.input_dir / "uploads" / principal.key_id
    user_uploads.mkdir(parents=True, exist_ok=True)
    target = user_uploads / f"{upload_id}_{Path(file.filename).name}"

    total = 0
    with target.open("wb") as out:
        while chunk := await file.read(64 * 1024):
            total += len(chunk)
            if total > MAX_UPLOAD_BYTES:
                out.close()
                target.unlink(missing_ok=True)
                raise HTTPException(
                    status_code=413,
                    detail=f"File exceeds {MAX_UPLOAD_BYTES // (1024 * 1024)} MB limit",
                )
            out.write(chunk)

    logger.info(
        "upload_received",
        upload_id=upload_id,
        filename=file.filename,
        size=total,
        api_key_id=principal.key_id,
    )
    return FileUploadResponse(upload_id=upload_id, filename=file.filename, size_bytes=total)


# ---------------------------------------------------------------------------
# Job CRUD
# ---------------------------------------------------------------------------


@router.post(
    "/jobs",
    response_model=JobCreateResponse,
    status_code=status.HTTP_202_ACCEPTED,
    summary="Submit a new primer-design job",
)
async def submit_job(
    body: JobCreate,
    request: Request,
    principal: ApiKeyPrincipal = Depends(require_api_key),
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> JobCreateResponse:
    """Validate, persist, and enqueue a primer-design job.

    Returns 202 with a job id you can poll via ``GET /api/v1/jobs/{id}``.
    Off-target FASTAs referenced by upload_id are copied from the upload
    area into the job's input directory at submission time so subsequent
    deletes don't disturb the in-flight job.
    """
    if not body.target.taxon_id and not body.target.accessions:
        raise HTTPException(
            status_code=422,
            detail="target.taxon_id or target.accessions must be provided.",
        )

    job_id = uuid.uuid4().hex
    in_dir, _ = _job_dirs(settings, job_id)

    # Copy any pre-uploaded files referenced by upload_id into the job.
    user_uploads = settings.input_dir / "uploads" / principal.key_id
    copied: list[str] = []
    for ref in body.off_targets.fasta_files:
        # ref may be an upload_id prefix (e.g. "abc123") or a full filename.
        candidates = (
            list(user_uploads.glob(f"{ref}_*")) if user_uploads.exists() else []
        )
        if not candidates:
            # Fallback: literal filename lookup.
            literal = user_uploads / ref
            if literal.exists():
                candidates = [literal]
        if not candidates:
            raise HTTPException(
                status_code=400,
                detail=(
                    f"Off-target reference {ref!r} not found in uploads. "
                    "POST the file to /api/v1/uploads first."
                ),
            )
        dest = in_dir / candidates[0].name
        shutil.copy2(candidates[0], dest)
        copied.append(dest.name)

    job = Job(
        id=job_id,
        api_key_id=principal.key_id,
        target_name=body.target.name,
        status=JobStatus.pending,
        config=body.model_dump(mode="json"),
        progress=[],
    )
    await create_job(session, job)

    arq_id = await enqueue_run_pipeline(request.app.state.redis_pool, job_id)
    job.arq_job_id = arq_id

    logger.info(
        "job_submitted",
        job_id=job_id,
        api_key_id=principal.key_id,
        target=body.target.name,
        off_target_files=copied,
    )

    return JobCreateResponse(
        id=job_id,
        status=JobStatus.pending,
        detail_url=f"{settings.base_url.rstrip('/')}/api/v1/jobs/{job_id}",
    )


@router.get("/jobs", response_model=list[JobSummary], summary="List jobs")
async def list_my_jobs(
    principal: ApiKeyPrincipal = Depends(require_api_key),
    session: AsyncSession = Depends(session_dep),
    limit: int = 50,
    offset: int = 0,
    status_filter: str | None = None,
) -> list[JobSummary]:
    """List jobs submitted by the calling API key, newest first.

    Admin keys see every job.
    """
    api_key_id = None if principal.is_admin else principal.key_id
    jobs = await list_jobs(
        session,
        api_key_id=api_key_id,
        limit=min(limit, 200),
        offset=max(offset, 0),
        status=status_filter,
    )
    return [_to_summary(j) for j in jobs]


@router.get("/jobs/{job_id}", response_model=JobDetail, summary="Job detail")
async def get_job_detail(
    job_id: str,
    principal: ApiKeyPrincipal = Depends(require_api_key),
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> JobDetail:
    """Return full job detail including config and progress timeline."""
    job = await get_job(session, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if not principal.is_admin and job.api_key_id != principal.key_id:
        raise HTTPException(status_code=404, detail="Job not found")
    return _to_detail(settings, job)


@router.post(
    "/jobs/{job_id}/cancel",
    response_model=JobDetail,
    summary="Request cancellation of a queued/running job",
)
async def cancel_job(
    job_id: str,
    request: Request,
    principal: ApiKeyPrincipal = Depends(require_api_key),
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> JobDetail:
    """Mark a job cancelled and abort the arq task if it's queued.

    For a job already running, the worker will stop at the next stage
    boundary (we check the DB status between stages).
    """
    job = await get_job(session, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if not principal.is_admin and job.api_key_id != principal.key_id:
        raise HTTPException(status_code=404, detail="Job not found")
    if job.status in (JobStatus.succeeded, JobStatus.failed, JobStatus.cancelled):
        return _to_detail(settings, job)

    await update_status(session, job_id, JobStatus.cancelled, error="Cancelled by user")
    # Best-effort: tell arq to abort if not yet picked up.
    try:
        from arq.jobs import Job as ArqJob

        if job.arq_job_id:
            arq_job = ArqJob(job.arq_job_id, request.app.state.redis_pool)
            await arq_job.abort(timeout=2)
    except Exception as e:  # noqa: BLE001
        logger.warning("arq_abort_failed", job_id=job_id, error=str(e))

    refreshed = await get_job(session, job_id)
    assert refreshed is not None
    return _to_detail(settings, refreshed)


@router.delete(
    "/jobs/{job_id}",
    status_code=status.HTTP_204_NO_CONTENT,
    summary="Delete a job and its artefacts",
)
async def delete_job_endpoint(
    job_id: str,
    principal: ApiKeyPrincipal = Depends(require_admin),
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> None:
    """Hard-delete a job's DB row and its result + input directories.

    Admin only — regular users can cancel but not delete (so audit history
    remains intact).
    """
    job = await get_job(session, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    await delete_job(session, job_id)
    for d in (settings.results_dir / job_id, settings.input_dir / job_id):
        if d.exists():
            shutil.rmtree(d, ignore_errors=True)
    logger.info("job_deleted", job_id=job_id, by=principal.key_id)

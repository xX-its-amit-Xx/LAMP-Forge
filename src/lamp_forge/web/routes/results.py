"""Download endpoints for completed-job artefacts."""

from __future__ import annotations

from fastapi import APIRouter, Depends, HTTPException
from fastapi.responses import FileResponse
from sqlalchemy.ext.asyncio import AsyncSession

from lamp_forge.web.auth import ApiKeyPrincipal, require_api_key
from lamp_forge.web.config import Settings, get_settings
from lamp_forge.web.db import JobStatus, get_job, session_dep

router = APIRouter()

# Whitelist of files a user may download. Keep this tight — the worker
# writes intermediate files we don't want clients pulling.
ALLOWED_FILES: dict[str, tuple[str, str]] = {
    "lamp_forge_report.html": ("text/html", "report"),
    "primer_sets.json": ("application/json", "json"),
    "primer_sets.csv": ("text/csv", "csv"),
    "conservation.tsv": ("text/tab-separated-values", "conservation"),
    "specificity.tsv": ("text/tab-separated-values", "specificity"),
    "run_manifest.json": ("application/json", "manifest"),
}


@router.get(
    "/jobs/{job_id}/results/{filename}",
    summary="Download a job result artefact",
    response_class=FileResponse,
)
async def download_artefact(
    job_id: str,
    filename: str,
    principal: ApiKeyPrincipal = Depends(require_api_key),
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> FileResponse:
    """Serve a single artefact file for a successful job.

    Returns 404 if the job doesn't belong to the calling key, or if the
    requested filename isn't on the whitelist.
    """
    job = await get_job(session, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    if not principal.is_admin and job.api_key_id != principal.key_id:
        raise HTTPException(status_code=404, detail="Job not found")
    if filename not in ALLOWED_FILES:
        raise HTTPException(
            status_code=404,
            detail=f"Artefact {filename!r} not available. Allowed: {sorted(ALLOWED_FILES)}",
        )
    if job.status != JobStatus.succeeded:
        raise HTTPException(
            status_code=409,
            detail=f"Job is {job.status}; results unavailable.",
        )

    media_type, _ = ALLOWED_FILES[filename]
    path = settings.results_dir / job_id / filename
    if not path.exists():
        raise HTTPException(status_code=404, detail=f"Artefact file missing on disk: {filename}")
    return FileResponse(path=path, media_type=media_type, filename=filename)

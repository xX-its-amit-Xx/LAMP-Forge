"""Server-rendered HTML routes (the user-facing dashboard).

We use Jinja2 templates + HTMX for progressive interactivity, no SPA. Three
reasons: (1) faster to ship and audit, (2) no frontend build pipeline, (3)
single deployable artefact. Tailwind via CDN keeps the templates readable.
"""

from __future__ import annotations

from pathlib import Path

from fastapi import APIRouter, Depends, HTTPException, Request
from fastapi.responses import HTMLResponse
from fastapi.templating import Jinja2Templates
from sqlalchemy.ext.asyncio import AsyncSession

from lamp_forge.web.config import Settings, get_settings
from lamp_forge.web.db import get_job, list_jobs, session_dep

router = APIRouter()

TEMPLATES_DIR = Path(__file__).resolve().parent.parent / "templates"
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))


@router.get("/", response_class=HTMLResponse, include_in_schema=False)
async def index(
    request: Request, settings: Settings = Depends(get_settings)
) -> HTMLResponse:
    """Landing page with links to docs, new-job form, and recent jobs."""
    return templates.TemplateResponse(
        request,
        "index.html",
        {"settings": settings, "title": "LAMP-Forge"},
    )


@router.get("/jobs", response_class=HTMLResponse, include_in_schema=False)
async def jobs_page(
    request: Request,
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> HTMLResponse:
    """Recent-jobs dashboard. Auth is delegated to the browser session: in
    dev mode (no API keys), shows all jobs; in prod, the user is expected
    to pass their API key in via the dashboard form.
    """
    api_key_id = request.query_params.get("key_id")
    jobs = await list_jobs(session, api_key_id=api_key_id, limit=50)
    return templates.TemplateResponse(
        request,
        "jobs_list.html",
        {
            "settings": settings,
            "jobs": jobs,
            "title": "Jobs — LAMP-Forge",
            "filter_key_id": api_key_id,
        },
    )


@router.get("/jobs/new", response_class=HTMLResponse, include_in_schema=False)
async def new_job_form(
    request: Request, settings: Settings = Depends(get_settings)
) -> HTMLResponse:
    """Show the new-job submission form."""
    return templates.TemplateResponse(
        request,
        "new_job.html",
        {"settings": settings, "title": "New job — LAMP-Forge"},
    )


@router.get("/jobs/{job_id}", response_class=HTMLResponse, include_in_schema=False)
async def job_detail_page(
    job_id: str,
    request: Request,
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> HTMLResponse:
    """Per-job detail page with live progress timeline."""
    job = await get_job(session, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return templates.TemplateResponse(
        request,
        "job_detail.html",
        {
            "settings": settings,
            "job": job,
            "title": f"Job {job_id[:8]} — LAMP-Forge",
        },
    )


@router.get(
    "/jobs/{job_id}/progress",
    response_class=HTMLResponse,
    include_in_schema=False,
)
async def job_progress_fragment(
    job_id: str,
    request: Request,
    session: AsyncSession = Depends(session_dep),
) -> HTMLResponse:
    """HTMX-polled fragment that renders the live progress timeline.

    The job detail page polls this endpoint every 2 seconds via HTMX until
    the job hits a terminal state, then stops polling.
    """
    job = await get_job(session, job_id)
    if job is None:
        raise HTTPException(status_code=404, detail="Job not found")
    return templates.TemplateResponse(
        request, "_progress_fragment.html", {"job": job}
    )

"""Liveness + readiness endpoints.

- ``/api/v1/health`` is liveness — always 200 when the process is running.
- ``/api/v1/ready`` checks every external dependency (Redis, DB, mafft,
  blast) and returns 503 if any are down. Wire this to your load balancer's
  health check.
"""

from __future__ import annotations

import shutil

from fastapi import APIRouter, Depends, Request
from sqlalchemy import text
from sqlalchemy.ext.asyncio import AsyncSession

from lamp_forge.web.config import Settings, get_settings
from lamp_forge.web.db import get_engine, session_dep
from lamp_forge.web.schemas import HealthResponse, ReadinessResponse

router = APIRouter()


@router.get("/health", response_model=HealthResponse)
async def health(settings: Settings = Depends(get_settings)) -> HealthResponse:
    """Liveness probe — succeeds as long as the process can serve HTTP."""
    return HealthResponse(version=settings.version, env=settings.env)


@router.get("/ready", response_model=ReadinessResponse, responses={503: {}})
async def ready(
    request: Request,
    settings: Settings = Depends(get_settings),
    session: AsyncSession = Depends(session_dep),
) -> ReadinessResponse:
    """Readiness probe — checks Redis, DB, and CLI tool availability."""
    statuses: dict[str, str] = {}

    # --- Redis ---
    try:
        pool = request.app.state.redis_pool
        await pool.ping()
        statuses["redis"] = "ok"
    except Exception as e:  # noqa: BLE001 — readiness must catch anything
        statuses["redis"] = f"down: {e.__class__.__name__}"

    # --- DB ---
    try:
        await session.execute(text("SELECT 1"))
        statuses["database"] = "ok"
    except Exception as e:  # noqa: BLE001
        statuses["database"] = f"down: {e.__class__.__name__}"

    # --- CLI binaries ---
    statuses["mafft"] = "ok" if shutil.which("mafft") else "missing"
    statuses["blast"] = (
        "ok"
        if (shutil.which("blastn") and shutil.which("makeblastdb"))
        else "missing"
    )

    overall = "ok" if all(v == "ok" for v in statuses.values()) else "degraded"
    response = ReadinessResponse(status=overall, **statuses)

    if overall != "ok":
        # Hint the load balancer to drain this replica without crashing it.
        # 503 makes ECS/Kubernetes pull the pod out of rotation.
        from fastapi import HTTPException

        raise HTTPException(status_code=503, detail=response.model_dump())

    return response

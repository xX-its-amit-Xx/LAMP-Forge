"""SQLite-backed job-state persistence (async via SQLAlchemy 2.x + aiosqlite).

We deliberately use SQLite, not Postgres, as the default. Justification:
- LAMP-Forge jobs are coarse-grained (one row written per state change, not
  per primer). Even at 100 jobs/day this is well under SQLite's write
  throughput.
- Zero-ops deployment matters for the typical user (a single research lab).
- WAL mode + a single writer (the FastAPI process plus the arq worker, both
  hitting the same file via aiosqlite) handles concurrency fine for this load.

For a multi-tenant SaaS deployment, swap to Postgres by setting
``LAMP_FORGE_DB_URL=postgresql+asyncpg://...`` — the model layer is
backend-agnostic.
"""

from __future__ import annotations

from collections.abc import AsyncGenerator
from datetime import UTC, datetime
from enum import StrEnum
from typing import Any

from sqlalchemy import JSON, DateTime, Integer, String, Text, select
from sqlalchemy.ext.asyncio import (
    AsyncEngine,
    AsyncSession,
    async_sessionmaker,
    create_async_engine,
)
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column

from lamp_forge.web.config import get_settings


class JobStatus(StrEnum):
    """Lifecycle of a primer-design job.

    Transitions: pending → running → (succeeded | failed | cancelled).
    """

    pending = "pending"
    running = "running"
    succeeded = "succeeded"
    failed = "failed"
    cancelled = "cancelled"


class Base(DeclarativeBase):
    """SQLAlchemy declarative base."""


def _utcnow() -> datetime:
    return datetime.now(UTC)


class Job(Base):
    """A LAMP-Forge pipeline job submitted by a user.

    One row per submission. The ``config`` column stores the validated
    YAML config as JSON so the same job can be replayed from history.
    Progress events are appended to ``progress`` so the UI can render a
    timeline without polling for log lines.
    """

    __tablename__ = "jobs"

    id: Mapped[str] = mapped_column(String(36), primary_key=True)
    api_key_id: Mapped[str] = mapped_column(String(64), index=True)
    target_name: Mapped[str] = mapped_column(String(200))
    status: Mapped[str] = mapped_column(String(20), index=True, default=JobStatus.pending)
    config: Mapped[dict[str, Any]] = mapped_column(JSON)
    progress: Mapped[list[dict[str, Any]]] = mapped_column(JSON, default=list)
    error: Mapped[str | None] = mapped_column(Text, nullable=True)
    n_primer_sets: Mapped[int | None] = mapped_column(Integer, nullable=True)
    n_conserved_regions: Mapped[int | None] = mapped_column(Integer, nullable=True)
    arq_job_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime, default=_utcnow, index=True)
    updated_at: Mapped[datetime] = mapped_column(DateTime, default=_utcnow, onupdate=_utcnow)
    started_at: Mapped[datetime | None] = mapped_column(DateTime, nullable=True)
    finished_at: Mapped[datetime | None] = mapped_column(DateTime, nullable=True)


# ---------------------------------------------------------------------------
# Engine / session management
# ---------------------------------------------------------------------------


_engine: AsyncEngine | None = None
_session_factory: async_sessionmaker[AsyncSession] | None = None


def _build_engine() -> AsyncEngine:
    settings = get_settings()
    settings.ensure_dirs()
    url = settings.resolved_db_url
    # echo=False because we have structlog covering query observability;
    # connect_args for sqlite are needed so we can share the connection
    # across async tasks safely.
    connect_args: dict[str, Any] = {}
    if url.startswith("sqlite"):
        connect_args["check_same_thread"] = False
    return create_async_engine(url, echo=False, future=True, connect_args=connect_args)


def get_engine() -> AsyncEngine:
    """Return the (lazily initialised) SQLAlchemy async engine."""
    global _engine
    if _engine is None:
        _engine = _build_engine()
    return _engine


def get_session_factory() -> async_sessionmaker[AsyncSession]:
    """Return the async session factory bound to the global engine."""
    global _session_factory
    if _session_factory is None:
        _session_factory = async_sessionmaker(
            get_engine(), expire_on_commit=False, class_=AsyncSession
        )
    return _session_factory


async def init_db() -> None:
    """Create tables if they don't exist. Idempotent — safe at every boot."""
    engine = get_engine()
    async with engine.begin() as conn:
        await conn.run_sync(Base.metadata.create_all)
        # Enable WAL for better concurrency between the FastAPI process and
        # the arq worker hitting the same sqlite file.
        if engine.url.get_backend_name().startswith("sqlite"):
            await conn.exec_driver_sql("PRAGMA journal_mode=WAL;")
            await conn.exec_driver_sql("PRAGMA synchronous=NORMAL;")


async def session_dep() -> AsyncGenerator[AsyncSession, None]:
    """FastAPI dependency yielding a session.

    Used as ``session: AsyncSession = Depends(session_dep)`` in routes.
    """
    factory = get_session_factory()
    async with factory() as session:
        try:
            yield session
            await session.commit()
        except Exception:
            await session.rollback()
            raise


# ---------------------------------------------------------------------------
# CRUD helpers
# ---------------------------------------------------------------------------


async def create_job(session: AsyncSession, job: Job) -> Job:
    """Insert a new job row and flush so the id is populated."""
    session.add(job)
    await session.flush()
    return job


async def get_job(session: AsyncSession, job_id: str) -> Job | None:
    """Look up a single job by id. Returns None if it doesn't exist."""
    return await session.get(Job, job_id)


async def list_jobs(
    session: AsyncSession,
    api_key_id: str | None,
    *,
    limit: int = 50,
    offset: int = 0,
    status: str | None = None,
) -> list[Job]:
    """List jobs for an API key (or all jobs if api_key_id is None — admin).

    Results are ordered newest-first.
    """
    stmt = select(Job).order_by(Job.created_at.desc())
    if api_key_id is not None:
        stmt = stmt.where(Job.api_key_id == api_key_id)
    if status:
        stmt = stmt.where(Job.status == status)
    stmt = stmt.limit(limit).offset(offset)
    res = await session.execute(stmt)
    return list(res.scalars().all())


async def update_status(
    session: AsyncSession,
    job_id: str,
    status: JobStatus,
    *,
    error: str | None = None,
) -> Job | None:
    """Atomically update a job's status (and optional error message)."""
    job = await session.get(Job, job_id)
    if job is None:
        return None
    job.status = status
    if status == JobStatus.running and job.started_at is None:
        job.started_at = _utcnow()
    if status in (JobStatus.succeeded, JobStatus.failed, JobStatus.cancelled):
        job.finished_at = _utcnow()
    if error is not None:
        job.error = error
    job.updated_at = _utcnow()
    return job


async def append_progress(
    session: AsyncSession,
    job_id: str,
    stage: str,
    message: str,
    percent: float | None = None,
) -> Job | None:
    """Append a progress event to the job's timeline."""
    job = await session.get(Job, job_id)
    if job is None:
        return None
    events = list(job.progress or [])
    events.append(
        {
            "ts": _utcnow().isoformat(),
            "stage": stage,
            "message": message,
            "percent": percent,
        }
    )
    job.progress = events
    job.updated_at = _utcnow()
    return job


async def finalize_metrics(
    session: AsyncSession,
    job_id: str,
    n_primer_sets: int,
    n_conserved_regions: int,
) -> Job | None:
    """Persist final-counts metrics for a successful job."""
    job = await session.get(Job, job_id)
    if job is None:
        return None
    job.n_primer_sets = n_primer_sets
    job.n_conserved_regions = n_conserved_regions
    job.updated_at = _utcnow()
    return job


async def delete_job(session: AsyncSession, job_id: str) -> bool:
    """Hard-delete a job row. Returns True if a row was deleted."""
    job = await session.get(Job, job_id)
    if job is None:
        return False
    await session.delete(job)
    return True

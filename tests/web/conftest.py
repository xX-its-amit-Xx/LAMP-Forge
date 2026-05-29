"""Shared fixtures for web tests.

We override the queue + DB + auth so unit tests run offline without a real
Redis or arq worker. The pipeline-task code path is exercised separately in
``test_tasks.py`` via direct invocation.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any
from unittest.mock import AsyncMock, MagicMock

import pytest
import pytest_asyncio
from httpx import ASGITransport, AsyncClient


@pytest.fixture(autouse=True)
def _set_dev_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Force a clean, isolated config for every web test."""
    monkeypatch.setenv("LAMP_FORGE_ENV", "development")
    monkeypatch.setenv("LAMP_FORGE_DATA_DIR", str(tmp_path / "data"))
    monkeypatch.setenv("LAMP_FORGE_DB_URL", f"sqlite+aiosqlite:///{tmp_path}/jobs.db")
    monkeypatch.setenv("LAMP_FORGE_REDIS_URL", "redis://localhost:6379/15")
    monkeypatch.setenv("LAMP_FORGE_API_KEYS", "")
    monkeypatch.setenv("LAMP_FORGE_LOG_JSON", "false")
    monkeypatch.setenv("LAMP_FORGE_ENABLE_METRICS", "true")
    # Clear cached settings between tests.
    from lamp_forge.web.config import get_settings

    get_settings.cache_clear()


@pytest_asyncio.fixture
async def app(monkeypatch: pytest.MonkeyPatch):
    """Build a FastAPI app with the arq pool mocked out."""
    from lamp_forge.web import app as app_module
    from lamp_forge.web.db import init_db

    # Avoid touching real Redis: patch the queue accessor + ping.
    fake_pool = MagicMock()
    fake_pool.ping = AsyncMock(return_value=True)
    fake_pool.aclose = AsyncMock()
    fake_pool.enqueue_job = AsyncMock(
        return_value=MagicMock(job_id="arq-test-job"), spec=None
    )

    async def _fake_get_pool() -> Any:
        return fake_pool

    monkeypatch.setattr(app_module, "get_arq_pool", _fake_get_pool)
    # Also patch the queue helper used by route handlers.
    monkeypatch.setattr(
        "lamp_forge.web.routes.jobs.enqueue_run_pipeline",
        AsyncMock(return_value="arq-test-job"),
    )

    fastapi_app = app_module.create_app()
    # Run startup hooks manually for tests that don't use the lifespan.
    await init_db()
    fastapi_app.state.redis_pool = fake_pool
    return fastapi_app


@pytest_asyncio.fixture
async def client(app):
    """An async httpx test client wired to the in-process app."""
    transport = ASGITransport(app=app)
    async with AsyncClient(transport=transport, base_url="http://test") as ac:
        yield ac

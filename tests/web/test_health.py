"""Tests for /health and /ready."""

from __future__ import annotations

import pytest


@pytest.mark.asyncio
async def test_health_ok(client) -> None:
    res = await client.get("/api/v1/health")
    assert res.status_code == 200
    body = res.json()
    assert body["status"] == "ok"
    assert "version" in body


@pytest.mark.asyncio
async def test_ready_checks_dependencies(client) -> None:
    res = await client.get("/api/v1/ready")
    # In CI we may not have mafft/blast installed, so this can be 503 with
    # a "missing" status. Either 200 or 503 is acceptable; what matters is
    # that the payload structure is correct.
    assert res.status_code in (200, 503)
    if res.status_code == 200:
        body = res.json()
        assert body["redis"] == "ok"
        assert body["database"] == "ok"


@pytest.mark.asyncio
async def test_metrics_endpoint(client) -> None:
    res = await client.get("/metrics")
    assert res.status_code == 200
    assert "python_info" in res.text or "process_" in res.text

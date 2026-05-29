"""Smoke tests for the HTML UI."""

from __future__ import annotations

import pytest


@pytest.mark.asyncio
async def test_index_renders(client) -> None:
    res = await client.get("/")
    assert res.status_code == 200
    assert "LAMP-Forge" in res.text


@pytest.mark.asyncio
async def test_new_job_form_renders(client) -> None:
    res = await client.get("/jobs/new")
    assert res.status_code == 200
    assert "New job" in res.text
    assert "Target name" in res.text


@pytest.mark.asyncio
async def test_jobs_list_renders_empty(client) -> None:
    res = await client.get("/jobs")
    assert res.status_code == 200
    assert "Jobs" in res.text
    assert "No jobs yet" in res.text


@pytest.mark.asyncio
async def test_unknown_job_detail_404(client) -> None:
    res = await client.get("/jobs/does-not-exist")
    assert res.status_code == 404


@pytest.mark.asyncio
async def test_openapi_schema_available(client) -> None:
    res = await client.get("/api/openapi.json")
    assert res.status_code == 200
    assert res.json()["info"]["title"] == "LAMP-Forge"

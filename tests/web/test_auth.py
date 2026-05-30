"""Auth tests — API key enforcement in production mode."""

from __future__ import annotations

import pytest


@pytest.fixture
def _prod_with_keys(monkeypatch: pytest.MonkeyPatch) -> str:
    """Switch to production-style auth with two keys (1 user, 1 admin)."""
    monkeypatch.setenv("LAMP_FORGE_ENV", "production")
    monkeypatch.setenv("LAMP_FORGE_API_KEYS", "test-user-key,test-admin-key")
    monkeypatch.setenv("LAMP_FORGE_ADMIN_API_KEYS", "test-admin-key")
    from lamp_forge.web.config import get_settings

    get_settings.cache_clear()
    return "test-user-key"


@pytest.mark.asyncio
async def test_unauthenticated_rejected_in_prod(_prod_with_keys, client) -> None:
    res = await client.get("/api/v1/jobs")
    assert res.status_code == 401


@pytest.mark.asyncio
async def test_invalid_key_rejected(_prod_with_keys, client) -> None:
    res = await client.get("/api/v1/jobs", headers={"X-API-Key": "wrong"})
    assert res.status_code == 401


@pytest.mark.asyncio
async def test_valid_key_accepted(_prod_with_keys, client) -> None:
    res = await client.get("/api/v1/jobs", headers={"Authorization": "Bearer test-user-key"})
    assert res.status_code == 200


@pytest.mark.asyncio
async def test_admin_delete_user_cannot(_prod_with_keys, client) -> None:
    # User key creates a job (using x-api-key header form for variety).
    body = {
        "target": {
            "name": "t",
            "taxon_id": 1773,
            "max_sequences": 1,
            "email": "a@b.com",
        }
    }
    create = await client.post("/api/v1/jobs", json=body, headers={"X-API-Key": "test-user-key"})
    assert create.status_code == 202
    job_id = create.json()["id"]

    # Non-admin DELETE is forbidden.
    res = await client.delete(f"/api/v1/jobs/{job_id}", headers={"X-API-Key": "test-user-key"})
    assert res.status_code == 403

    # Admin DELETE works.
    res = await client.delete(f"/api/v1/jobs/{job_id}", headers={"X-API-Key": "test-admin-key"})
    assert res.status_code == 204

"""Tests for the jobs API.

We don't actually run the pipeline here — that's covered separately. What
matters is that the HTTP layer validates input, persists rows, and
enforces auth correctly.
"""

from __future__ import annotations

import io

import pytest


def _valid_job_body() -> dict:
    return {
        "target": {
            "name": "test_target",
            "taxon_id": 1773,
            "gene": "rpoB",
            "max_sequences": 5,
            "email": "user@test.edu",
        },
        "off_targets": {
            "min_identity_threshold": 0.80,
            "min_coverage_threshold": 0.80,
            "fasta_files": [],
        },
        "conservation": {
            "window_size": 30,
            "entropy_threshold": 0.25,
            "min_region_length": 250,
        },
        "primer": {
            "tm_min": 60.0,
            "tm_max": 65.0,
            "tm_match_tolerance": 2.0,
            "gc_min": 40.0,
            "gc_max": 65.0,
            "hairpin_dg_threshold": -2.0,
            "dimer_dg_threshold": -5.0,
            "amplicon_size": {"f2_b2_min": 120, "f2_b2_max": 160},
        },
        "output": {"top_n": 10, "generate_html": True, "generate_csv": True},
    }


@pytest.mark.asyncio
async def test_submit_job_returns_202(client) -> None:
    res = await client.post("/api/v1/jobs", json=_valid_job_body())
    assert res.status_code == 202, res.text
    data = res.json()
    assert "id" in data
    assert data["status"] == "pending"
    assert data["detail_url"].endswith(f"/api/v1/jobs/{data['id']}")


@pytest.mark.asyncio
async def test_submit_then_get_detail(client) -> None:
    create = await client.post("/api/v1/jobs", json=_valid_job_body())
    job_id = create.json()["id"]
    detail = await client.get(f"/api/v1/jobs/{job_id}")
    assert detail.status_code == 200
    body = detail.json()
    assert body["target_name"] == "test_target"
    assert body["status"] in ("pending", "running")
    assert body["config"]["target"]["taxon_id"] == 1773
    assert body["progress"] == []


@pytest.mark.asyncio
async def test_submit_rejects_missing_target_identifier(client) -> None:
    body = _valid_job_body()
    body["target"]["taxon_id"] = None
    body["target"]["accessions"] = []
    res = await client.post("/api/v1/jobs", json=body)
    assert res.status_code == 422
    assert "taxon_id" in res.text


@pytest.mark.asyncio
async def test_submit_rejects_inverted_tm_band(client) -> None:
    body = _valid_job_body()
    body["primer"]["tm_min"] = 70.0
    res = await client.post("/api/v1/jobs", json=body)
    assert res.status_code == 422


@pytest.mark.asyncio
async def test_list_jobs(client) -> None:
    await client.post("/api/v1/jobs", json=_valid_job_body())
    res = await client.get("/api/v1/jobs")
    assert res.status_code == 200
    items = res.json()
    assert len(items) >= 1


@pytest.mark.asyncio
async def test_get_unknown_job_returns_404(client) -> None:
    res = await client.get("/api/v1/jobs/does-not-exist")
    assert res.status_code == 404


@pytest.mark.asyncio
async def test_cancel_job(client) -> None:
    create = await client.post("/api/v1/jobs", json=_valid_job_body())
    job_id = create.json()["id"]
    cancel = await client.post(f"/api/v1/jobs/{job_id}/cancel")
    assert cancel.status_code == 200
    assert cancel.json()["status"] == "cancelled"


@pytest.mark.asyncio
async def test_upload_fasta_file(client) -> None:
    fasta = b">test\nACGTACGTACGTACGT\n"
    files = {"file": ("test.fasta", io.BytesIO(fasta), "text/plain")}
    res = await client.post("/api/v1/uploads", files=files)
    assert res.status_code == 201
    body = res.json()
    assert body["filename"] == "test.fasta"
    assert body["size_bytes"] == len(fasta)
    assert "upload_id" in body


@pytest.mark.asyncio
async def test_upload_rejects_bad_extension(client) -> None:
    files = {"file": ("test.txt", io.BytesIO(b"not a fasta"), "text/plain")}
    res = await client.post("/api/v1/uploads", files=files)
    assert res.status_code == 400


@pytest.mark.asyncio
async def test_download_artefact_requires_completed_job(client) -> None:
    create = await client.post("/api/v1/jobs", json=_valid_job_body())
    job_id = create.json()["id"]
    res = await client.get(f"/api/v1/jobs/{job_id}/results/primer_sets.csv")
    assert res.status_code == 409  # not succeeded yet


@pytest.mark.asyncio
async def test_download_unknown_artefact_rejected(client) -> None:
    create = await client.post("/api/v1/jobs", json=_valid_job_body())
    job_id = create.json()["id"]
    res = await client.get(f"/api/v1/jobs/{job_id}/results/etc_passwd")
    assert res.status_code == 404

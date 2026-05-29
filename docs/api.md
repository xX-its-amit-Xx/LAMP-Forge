# LAMP-Forge HTTP API

OpenAPI / interactive: `https://<host>/api/docs` (Swagger) or `/api/redoc`.

All endpoints under `/api/v1/*` are versioned. The HTML UI under `/` is
unversioned and may change without notice.

## Auth

Pass your API key either way:

```
Authorization: Bearer lf_xxxxxxxxx...
X-API-Key: lf_xxxxxxxxx...
```

In development (no `LAMP_FORGE_API_KEYS` configured), auth is bypassed and
requests are attributed to an "anonymous" principal.

Admin keys (configured separately via `LAMP_FORGE_ADMIN_API_KEYS`) can
DELETE jobs and list jobs across all users.

## Rate limits

- Default: `30 requests/minute` per API key (configurable).
- Job submission: `20/hour` per API key.
- 429 responses include the standard `Retry-After` header.

## Endpoints

### `POST /api/v1/uploads`

Upload an off-target FASTA file. Returns an `upload_id` to reference in a
later job submission.

```bash
curl -X POST https://lamp-forge.example.org/api/v1/uploads \
     -H "Authorization: Bearer $LF_KEY" \
     -F "file=@m_smegmatis.fasta"
```

```json
{
  "upload_id": "8a3f...",
  "filename": "m_smegmatis.fasta",
  "size_bytes": 2451323
}
```

Allowed extensions: `.fa`, `.fasta`, `.fna`, `.ffn`. Max 50 MB per file.

### `POST /api/v1/jobs`

Submit a new primer-design job.

```bash
curl -X POST https://lamp-forge.example.org/api/v1/jobs \
     -H "Authorization: Bearer $LF_KEY" \
     -H "Content-Type: application/json" \
     -d @body.json
```

Body schema (full reference in OpenAPI):

```jsonc
{
  "target": {
    "name": "mtb_rpoB",
    "taxon_id": 1773,           // OR "accessions": ["NC_000962.3", ...]
    "gene": "rpoB",
    "max_sequences": 20,
    "email": "you@inst.edu"
  },
  "off_targets": {
    "min_identity_threshold": 0.80,
    "min_coverage_threshold": 0.80,
    "fasta_files": ["8a3f...", "9b4e..."]   // upload_ids from /uploads
  },
  "conservation": {
    "window_size": 30,
    "entropy_threshold": 0.25,
    "min_region_length": 250
  },
  "primer": {
    "tm_min": 60.0, "tm_max": 65.0, "tm_match_tolerance": 2.0,
    "gc_min": 40.0, "gc_max": 65.0,
    "hairpin_dg_threshold": -2.0, "dimer_dg_threshold": -5.0,
    "amplicon_size": { "f2_b2_min": 120, "f2_b2_max": 160 }
  },
  "output": {
    "top_n": 10,
    "generate_html": true,
    "generate_csv": true
  }
}
```

Response (202 Accepted):

```json
{
  "id": "ab1c2d3e4f...",
  "status": "pending",
  "detail_url": "https://lamp-forge.example.org/api/v1/jobs/ab1c2d3e4f..."
}
```

### `GET /api/v1/jobs/{id}`

Get full detail including config, progress timeline, and result URLs.

```jsonc
{
  "id": "ab1c2d3e...",
  "target_name": "mtb_rpoB",
  "status": "succeeded",         // pending | running | succeeded | failed | cancelled
  "created_at": "2026-05-28T14:33:01Z",
  "started_at": "2026-05-28T14:33:04Z",
  "finished_at": "2026-05-28T14:35:48Z",
  "n_primer_sets": 47,
  "n_conserved_regions": 8,
  "config": { /* full config snapshot */ },
  "progress": [
    { "ts": "...", "stage": "fetch", "message": "...", "percent": 5 },
    { "ts": "...", "stage": "align", "message": "...", "percent": 20 },
    /* ... */
  ],
  "error": null,
  "result_urls": {
    "html": "https://.../jobs/ab1.../results/lamp_forge_report.html",
    "json": "https://.../jobs/ab1.../results/primer_sets.json",
    "csv":  "https://.../jobs/ab1.../results/primer_sets.csv",
    "manifest": "https://.../jobs/ab1.../results/run_manifest.json",
    "conservation_tsv": "https://.../jobs/ab1.../results/conservation.tsv",
    "specificity_tsv":  "https://.../jobs/ab1.../results/specificity.tsv"
  }
}
```

### `GET /api/v1/jobs`

List jobs for the calling API key (or all jobs, for admin keys). Newest first.

Query params:
- `limit` (default 50, max 200)
- `offset`
- `status_filter` (`pending` | `running` | `succeeded` | `failed` | `cancelled`)

### `POST /api/v1/jobs/{id}/cancel`

Request cancellation. If the job is still queued, it's removed from the
queue immediately; if it's running, the worker stops at the next stage
boundary (typically within a few seconds).

### `DELETE /api/v1/jobs/{id}`

Admin-only. Hard-delete the job row and its on-disk artefacts.

### `GET /api/v1/jobs/{id}/results/{filename}`

Download a single artefact. Allowed filenames:

- `lamp_forge_report.html`
- `primer_sets.json`
- `primer_sets.csv`
- `conservation.tsv`
- `specificity.tsv`
- `run_manifest.json`

Returns 409 if the job hasn't succeeded yet. Returns 404 if `filename` is
not in the whitelist (deliberate — protects against path-traversal abuse).

### `GET /api/v1/health`

Liveness. Always 200 if the process is up.

### `GET /api/v1/ready`

Readiness. Returns 200 only if Redis, the DB, MAFFT, and BLAST+ are all
reachable. 503 otherwise — wire to your load balancer's health check.

### `GET /metrics`

Prometheus scrape endpoint. Disabled if `LAMP_FORGE_ENABLE_METRICS=false`.

## Errors

All errors use:

```json
{
  "detail": "Human-readable message",
  "code": "optional_machine_code"
}
```

Common codes:

| Status | When |
|---|---|
| 400 | Malformed upload or unknown reference id |
| 401 | Missing / invalid API key |
| 403 | Authenticated but lacking permission (admin endpoints) |
| 404 | Resource doesn't exist or is not yours |
| 409 | Job not in a state that supports the requested action |
| 413 | Upload exceeds size limit |
| 422 | Pydantic validation failure (with field paths in the detail) |
| 429 | Rate limit hit |
| 503 | Readiness probe failed; one or more dependencies down |

## Python client example

```python
import httpx, time

BASE = "https://lamp-forge.example.org"
KEY = "lf_xxxxx..."
H = {"Authorization": f"Bearer {KEY}"}

with httpx.Client(base_url=BASE, headers=H, timeout=30) as c:
    # 1. Upload off-target panels.
    upload_ids = []
    for p in ["m_smegmatis.fasta", "m_avium.fasta"]:
        r = c.post("/api/v1/uploads", files={"file": open(p, "rb")})
        r.raise_for_status()
        upload_ids.append(r.json()["upload_id"])

    # 2. Submit job.
    r = c.post("/api/v1/jobs", json={
        "target": {
            "name": "mtb_rpoB", "taxon_id": 1773,
            "gene": "rpoB", "max_sequences": 20,
            "email": "you@inst.edu",
        },
        "off_targets": {"fasta_files": upload_ids},
    })
    r.raise_for_status()
    job_id = r.json()["id"]

    # 3. Poll.
    while True:
        d = c.get(f"/api/v1/jobs/{job_id}").json()
        print(d["status"], d["progress"][-1] if d["progress"] else "")
        if d["status"] in {"succeeded", "failed", "cancelled"}:
            break
        time.sleep(3)

    # 4. Download CSV.
    if d["status"] == "succeeded":
        r = c.get(d["result_urls"]["csv"])
        open(f"{job_id}.csv", "wb").write(r.content)
```

# LAMP-Forge deployment guide

The web service is built so you can stand it up on a single VM for a lab,
or scale it horizontally across a kubernetes cluster for an institute-wide
deployment. This guide covers both.

## Architecture refresher

```
                   ┌──────────────┐
   HTTPS users ───►│  nginx       │── TLS termination, rate limit, static
                   │  (sidecar)   │
                   └──────┬───────┘
                          │ HTTP
                   ┌──────▼───────┐
                   │  FastAPI web │── REST + HTMX UI, lifespan-managed
                   │  (1-N pods)  │   Redis pool, SQLAlchemy session
                   └──────┬───────┘
                          │
                   ┌──────▼───────┐         ┌─────────────────┐
                   │   Redis      │◄────────┤ arq workers (N) │── run pipeline
                   │   (queue)    │         │                 │   stages
                   └──────────────┘         └─────────────────┘
                                                      │
                                                      ▼
                                         ┌────────────────────────┐
                                         │ Persistent volume      │── SQLite,
                                         │ /work/data             │   cache,
                                         │   /results/<job_id>/   │   results
                                         │   /cache/              │
                                         │   /input/<job_id>/     │
                                         │   /jobs.db (SQLite)    │
                                         └────────────────────────┘
```

Three runtime images, all built from the same `Dockerfile`:

| Image | Built with | Default ENTRYPOINT |
|---|---|---|
| `lamp-forge:cli` | `--target cli` | `lamp-forge` (single-shot) |
| `lamp-forge:web` | `--target web` | `lamp-forge-web` on :8000 |
| `lamp-forge:worker` | `--target worker` | `lamp-forge-worker` |

## Single-VM Docker Compose (the common case)

```bash
cp .env.example .env
# Edit .env: set LAMP_FORGE_API_KEYS to one or more strong random strings.
# Example:
python -c "import secrets; print(secrets.token_urlsafe(32))"

docker compose up -d redis web worker
```

Add nginx + TLS for production:

```bash
# 1. Drop a TLS cert into deploy/nginx/certs/ (fullchain.pem + privkey.pem)
# 2. Uncomment the HTTPS server block in deploy/nginx/nginx.conf
# 3. Up the prod profile.
docker compose --profile prod up -d
```

Scale workers:

```bash
LAMP_FORGE_WORKER_REPLICAS=4 docker compose up -d worker
```

## Kubernetes

A minimal manifest set lives in [deploy/k8s/](../deploy/k8s/) (skeletons —
adapt to your cluster's ingress / secret / volume conventions):

- `web-deployment.yaml` — 2-replica HPA-friendly Deployment
- `worker-deployment.yaml` — N-replica Deployment, no HPA (queue-driven)
- `redis-statefulset.yaml` — single-replica StatefulSet with PV
- `pvc.yaml` — PersistentVolumeClaim for `/work/data`
- `ingress.yaml` — Ingress with TLS via cert-manager
- `configmap.yaml`, `secret.yaml` — config + API keys

Key things to remember in k8s:

1. **Mount the same PVC** to the web and worker deployments (`/work/data`).
   Results land here; both must see them.
2. **Web is HPA-friendly**, workers aren't — workers are queue-driven, so
   scale them with a fixed replica count or KEDA based on Redis queue depth.
3. **Use a managed Postgres for the DB in multi-tenant prod** — set
   `LAMP_FORGE_DB_URL=postgresql+asyncpg://...`. SQLite is fine for a
   single-instance lab deployment but doesn't scale across multiple web
   pods cleanly.
4. **API keys go in a Secret** mounted as `LAMP_FORGE_API_KEYS`.
5. **The readiness probe** is `GET /api/v1/ready` — it checks Redis, the DB,
   MAFFT, and BLAST+. Use it as both readinessProbe and startup gate.

## Observability

### Logs

All services emit structured JSON via structlog. Every log line carries
`request_id` (HTTP requests) and `job_id` (worker tasks) so you can trace a
single submission end-to-end across the web + worker logs:

```bash
# example with jq + docker compose logs
docker compose logs --tail=200 web worker \
  | jq 'select(.job_id == "abc123def456")'
```

### Metrics

`GET /metrics` exposes a Prometheus scrape endpoint with FastAPI's default
metrics (request counts, latencies, in-flight) plus python_info. Add
custom metrics in `lamp_forge.web.app` if you want per-stage timings.

A starter Grafana dashboard JSON lives in
[deploy/grafana/dashboard.json](../deploy/grafana/dashboard.json).

### Health probes

- `GET /api/v1/health` — liveness (200 if process alive)
- `GET /api/v1/ready` — readiness (200 if all deps OK, 503 otherwise)

## Security checklist before going to production

- [ ] `LAMP_FORGE_ENV=production` is set.
- [ ] `LAMP_FORGE_API_KEYS` has at least one strong random value;
      `LAMP_FORGE_ADMIN_API_KEYS` is set separately and limited.
- [ ] TLS termination is configured (nginx with valid cert, or load balancer).
- [ ] `LAMP_FORGE_CORS_ORIGINS` is restricted to your actual frontend
      origins, not `*`.
- [ ] Rate limits (`LAMP_FORGE_RATE_LIMIT_PER_MINUTE`) match your usage.
- [ ] The persistent volume is backed up.
- [ ] Old jobs are cleaned up — wire up the retention cron (see
      [deploy/cron/cleanup.sh](../deploy/cron/cleanup.sh)).
- [ ] NCBI_API_KEY is set so the worker doesn't hit the 3 req/sec limit.
- [ ] The container runs as the non-root `forge` user (verify with
      `docker exec lamp-forge-web id`).

## Backup and disaster recovery

The only stateful data is in the `/work/data` volume:

- `jobs.db` — SQLite job metadata
- `cache/` — NCBI fetch cache (recoverable, just slow to rebuild)
- `results/` — per-job artefacts
- `input/` — per-job uploaded off-target panels

Back up the whole volume nightly. To restore, mount the snapshot at
`/work/data` and start the stack — the schema is idempotent and the workers
will pick up any queued-but-unfinished jobs.

## Upgrading

LAMP-Forge follows semver. Within a major version, upgrades are
schema-compatible — just pull the new image and `docker compose up -d`.

For major-version upgrades, run `lamp-forge-web --check-migrations` first
(it dry-runs the schema diff). If a migration is required, follow the
notes in `CHANGELOG.md`.

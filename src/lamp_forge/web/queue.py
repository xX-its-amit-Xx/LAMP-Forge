"""arq queue client + helpers.

The same Redis pool is shared by the API process (enqueues jobs, reads
results) and the worker process (consumes the queue). Keeping the
construction in one place makes it easy to mock for tests.
"""

from __future__ import annotations

from arq.connections import ArqRedis, RedisSettings, create_pool

from lamp_forge.web.config import Settings, get_settings


def redis_settings(settings: Settings | None = None) -> RedisSettings:
    """Build an arq RedisSettings from the web service config."""
    s = settings or get_settings()
    return RedisSettings.from_dsn(s.redis_url)


async def get_arq_pool() -> ArqRedis:
    """Create an arq Redis pool. Cache the result at app scope."""
    return await create_pool(redis_settings())


async def enqueue_run_pipeline(pool: ArqRedis, job_id: str) -> str:
    """Submit a pipeline-run task to the worker queue.

    Returns the arq job id (distinct from our LAMP-Forge job_id) so the API
    can store it for later cancel/inspect calls.
    """
    arq_job = await pool.enqueue_job("run_pipeline", job_id, _job_id=f"lf-{job_id}")
    if arq_job is None:
        raise RuntimeError(
            "arq returned None when enqueueing job; this usually means the "
            "job id is already in use. Check the Redis queue state."
        )
    return arq_job.job_id

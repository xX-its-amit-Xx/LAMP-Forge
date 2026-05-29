"""arq worker entry point.

Run with:
    arq lamp_forge.web.worker.WorkerSettings
or via the bundled console script ``lamp-forge-worker``.

We initialise the DB at startup so the worker can serve from a fresh
container without depending on the API process to have run first.
"""

from __future__ import annotations

import asyncio

from arq.connections import RedisSettings

from lamp_forge.web.config import get_settings
from lamp_forge.web.db import init_db
from lamp_forge.web.logging_setup import configure_logging, get_logger
from lamp_forge.web.queue import redis_settings
from lamp_forge.web.tasks import run_pipeline

logger = get_logger("lamp_forge.web.worker")


async def on_startup(ctx: dict) -> None:  # type: ignore[type-arg]
    """arq startup hook — bootstrap logging + DB once per worker process."""
    settings = get_settings()
    configure_logging(level=settings.log_level, json_logs=settings.log_json)
    settings.ensure_dirs()
    await init_db()
    logger.info("worker_startup", env=settings.env)


async def on_shutdown(ctx: dict) -> None:  # type: ignore[type-arg]
    """arq shutdown hook."""
    logger.info("worker_shutdown")


class WorkerSettings:
    """arq worker configuration.

    Tunables worth knowing about:
    - ``max_jobs``: how many concurrent pipeline runs a single worker handles.
      The pipeline is CPU-light but memory-medium (MAFFT, BLAST, primer3);
      4 concurrent is comfortable on a 4-vCPU node. Scale horizontally by
      adding worker replicas instead of cranking this up.
    - ``job_timeout``: hard ceiling per job; matches the configured
      ``LAMP_FORGE_JOB_TIMEOUT_SECONDS``.
    """

    functions = [run_pipeline]
    on_startup = on_startup
    on_shutdown = on_shutdown
    redis_settings: RedisSettings = redis_settings()
    max_jobs = 4
    keep_result = 3600  # how long arq retains the result blob
    job_timeout = get_settings().job_timeout_seconds


def main() -> None:
    """Console-script entrypoint that runs the arq worker.

    Equivalent to ``arq lamp_forge.web.worker.WorkerSettings`` but doesn't
    require the ``arq`` CLI to be on PATH.
    """
    from arq.worker import run_worker

    async def _run() -> None:
        run_worker(WorkerSettings)  # type: ignore[arg-type]

    asyncio.run(_run())


if __name__ == "__main__":
    main()

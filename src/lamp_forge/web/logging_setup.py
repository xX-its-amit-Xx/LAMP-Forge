"""Structured logging configuration via structlog.

Production gets JSON logs (ingestible by Datadog/Splunk/Loki); dev gets a
human-readable console renderer. Either way, every log line carries the
request_id, job_id (when applicable), and api_key_id (hashed) so a single
job can be traced end-to-end from HTTP request to worker completion.
"""

from __future__ import annotations

import logging
import sys
from typing import cast

import structlog


def configure_logging(level: str = "INFO", json_logs: bool = True) -> None:
    """Configure structlog + stdlib logging.

    Args:
        level: Root log level (DEBUG / INFO / WARNING / ERROR).
        json_logs: If True emit JSON; if False emit pretty console output.
    """
    timestamper = structlog.processors.TimeStamper(fmt="iso", utc=True)
    shared_processors: list[structlog.types.Processor] = [
        structlog.contextvars.merge_contextvars,
        structlog.processors.add_log_level,
        structlog.stdlib.add_logger_name,
        timestamper,
        structlog.processors.StackInfoRenderer(),
    ]

    if json_logs:
        renderer: structlog.types.Processor = structlog.processors.JSONRenderer()
    else:
        renderer = structlog.dev.ConsoleRenderer(colors=True)

    structlog.configure(
        processors=[
            *shared_processors,
            structlog.processors.format_exc_info,
            renderer,
        ],
        wrapper_class=structlog.stdlib.BoundLogger,
        context_class=dict,
        logger_factory=structlog.stdlib.LoggerFactory(),
        cache_logger_on_first_use=True,
    )

    # Route stdlib logs through the same handler so third-party libs
    # (uvicorn, sqlalchemy, biopython) emit consistently formatted lines.
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(logging.Formatter("%(message)s"))
    root = logging.getLogger()
    root.handlers = [handler]
    root.setLevel(level.upper())

    # Quiet down a couple of noisy loggers at INFO; raise to DEBUG in dev.
    for noisy in ("uvicorn.access", "sqlalchemy.engine.Engine"):
        logging.getLogger(noisy).setLevel(logging.WARNING)


def get_logger(name: str | None = None) -> structlog.stdlib.BoundLogger:
    """Convenience wrapper around ``structlog.get_logger``."""
    return cast(structlog.stdlib.BoundLogger, structlog.get_logger(name))

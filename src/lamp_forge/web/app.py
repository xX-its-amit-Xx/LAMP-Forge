"""FastAPI application factory.

We use a factory rather than a module-level ``app = FastAPI()`` so tests can
build a fresh app per test with overridden settings, and so the same code
path constructs prod and test instances identically.
"""

from __future__ import annotations

import uuid
from collections.abc import AsyncIterator
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from fastapi.staticfiles import StaticFiles
from prometheus_client import CONTENT_TYPE_LATEST, generate_latest
from slowapi import Limiter
from slowapi.errors import RateLimitExceeded
from slowapi.util import get_remote_address
from starlette.middleware.base import RequestResponseEndpoint
from starlette.responses import Response

from lamp_forge.web.config import Settings, get_settings
from lamp_forge.web.db import init_db
from lamp_forge.web.logging_setup import configure_logging, get_logger
from lamp_forge.web.queue import get_arq_pool
from lamp_forge.web.routes import health, jobs, results, ui

logger = get_logger("lamp_forge.web.app")


@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncIterator[None]:
    """Startup: init DB, open Redis pool. Shutdown: close pool cleanly."""
    settings: Settings = get_settings()
    configure_logging(level=settings.log_level, json_logs=settings.log_json)
    settings.ensure_dirs()
    await init_db()
    app.state.redis_pool = await get_arq_pool()
    logger.info("web_startup", env=settings.env, version=settings.version)
    try:
        yield
    finally:
        await app.state.redis_pool.aclose()
        logger.info("web_shutdown")


def create_app(settings: Settings | None = None) -> FastAPI:
    """Build a configured FastAPI app.

    Args:
        settings: Optional explicit settings (used by tests). Production
            should let this default to :func:`get_settings`.
    """
    s = settings or get_settings()
    app = FastAPI(
        title="LAMP-Forge",
        description=(
            "Reproducible LAMP primer design as a service. Submit a target "
            "organism + off-target panel, get back a ranked, ready-to-order "
            "set of LAMP primer sets."
        ),
        version=s.version,
        debug=s.debug,
        lifespan=lifespan,
        docs_url="/api/docs",
        redoc_url="/api/redoc",
        openapi_url="/api/openapi.json",
    )

    # --- Middleware ---
    app.add_middleware(
        CORSMiddleware,
        allow_origins=s.cors_origins,
        allow_credentials=True,
        allow_methods=["*"],
        allow_headers=["*"],
    )

    # Rate limiting: keyed on API key when present, falling back to client IP.
    def _rate_key(request: Request) -> str:
        api_key = request.headers.get("X-API-Key") or request.headers.get("Authorization", "")
        return api_key or get_remote_address(request)

    limiter = Limiter(key_func=_rate_key, default_limits=[f"{s.rate_limit_per_minute}/minute"])
    app.state.limiter = limiter

    @app.exception_handler(RateLimitExceeded)
    async def _rate_limit_handler(request: Request, exc: RateLimitExceeded) -> Response:
        return JSONResponse(
            status_code=429,
            content={
                "detail": f"Rate limit exceeded: {exc.detail}",
                "code": "rate_limit_exceeded",
            },
        )

    # --- Request ID middleware (sets a correlation id on every request) ---
    @app.middleware("http")
    async def _request_id_middleware(
        request: Request, call_next: RequestResponseEndpoint
    ) -> Response:
        request_id = request.headers.get("X-Request-ID") or str(uuid.uuid4())
        request.state.request_id = request_id
        # Bind to structlog contextvars so every log line in this request
        # automatically carries the request_id.
        from structlog import contextvars

        contextvars.clear_contextvars()
        contextvars.bind_contextvars(
            request_id=request_id, path=request.url.path, method=request.method
        )
        try:
            response = await call_next(request)
        finally:
            contextvars.clear_contextvars()
        response.headers["X-Request-ID"] = request_id
        return response

    # --- Routes ---
    app.include_router(health.router, prefix="/api/v1", tags=["health"])
    app.include_router(jobs.router, prefix="/api/v1", tags=["jobs"])
    app.include_router(results.router, prefix="/api/v1", tags=["results"])
    app.include_router(ui.router, tags=["ui"])

    # Static assets for the HTML UI.
    static_dir = __import__("pathlib").Path(__file__).resolve().parent / "static"
    if static_dir.exists():
        app.mount("/static", StaticFiles(directory=static_dir), name="static")

    # --- Metrics ---
    if s.enable_metrics:

        @app.get("/metrics", include_in_schema=False)
        async def metrics() -> Response:
            return Response(content=generate_latest(), media_type=CONTENT_TYPE_LATEST)

    return app

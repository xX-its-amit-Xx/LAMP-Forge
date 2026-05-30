"""Web service configuration via environment variables.

All settings are loaded from the environment (or a ``.env`` file) so the
same container image can be deployed across dev / staging / prod with only
env overrides. Defaults are dev-safe; production deployments MUST set
``API_KEYS`` to at least one strong random value.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Annotated

from pydantic import Field, field_validator
from pydantic_settings import BaseSettings, NoDecode, SettingsConfigDict


class Settings(BaseSettings):
    """Runtime settings for the LAMP-Forge web service.

    All fields are populated from the environment with the ``LAMP_FORGE_``
    prefix (e.g. ``LAMP_FORGE_REDIS_URL``). A ``.env`` file in the working
    directory is also loaded if present.
    """

    model_config = SettingsConfigDict(
        env_prefix="LAMP_FORGE_",
        env_file=".env",
        env_file_encoding="utf-8",
        extra="ignore",
        case_sensitive=False,
    )

    # ---- General ----
    env: str = Field("development", description="Deployment environment label.")
    debug: bool = Field(False, description="Enable FastAPI debug mode and verbose tracebacks.")
    version: str = "0.1.0"
    base_url: str = Field(
        "http://localhost:8000",
        description="Public base URL; used for absolute links in emails/notifications.",
    )

    # ---- Auth ----
    api_keys: Annotated[
        list[str],
        NoDecode,
        Field(
            default_factory=list,
            description=(
                "Whitelist of accepted API keys. Pass as comma-separated env var. "
                "If empty AND env != 'production', auth is disabled (dev convenience)."
            ),
        ),
    ]
    admin_api_keys: Annotated[
        list[str],
        NoDecode,
        Field(
            default_factory=list,
            description="Subset of api_keys with admin privileges (delete jobs, view all).",
        ),
    ]

    # ---- Storage ----
    data_dir: Path = Field(
        Path("/work/data"),
        description="Persistent root for the job DB, fetch cache, and results.",
    )
    db_url: str | None = Field(
        None,
        description="SQLAlchemy URL override. Defaults to sqlite at {data_dir}/jobs.db.",
    )

    # ---- Redis / queue ----
    redis_url: str = Field("redis://localhost:6379/0", description="Redis URL for arq queue.")
    job_timeout_seconds: int = Field(1800, description="Hard timeout per job (30 min by default).")
    job_retention_days: int = Field(
        14, description="How long completed job artefacts are kept before cleanup."
    )

    # ---- Rate limiting ----
    rate_limit_per_minute: int = Field(
        30, description="Per-API-key request ceiling, minute window."
    )
    job_submit_rate_limit_per_hour: int = Field(
        20, description="Per-API-key job submissions, hour window."
    )

    # ---- Observability ----
    log_level: str = "INFO"
    log_json: bool = Field(
        True,
        description="Emit structured JSON logs. Disable for human-readable dev logs.",
    )
    enable_metrics: bool = Field(True, description="Expose Prometheus /metrics endpoint.")

    # ---- CORS ----
    cors_origins: Annotated[
        list[str],
        NoDecode,
        Field(
            default_factory=lambda: ["http://localhost:8000"],
            description="Allowed CORS origins.",
        ),
    ]

    # ---- NCBI passthrough ----
    ncbi_email: str = ""
    ncbi_api_key: str = ""

    @field_validator("api_keys", "admin_api_keys", "cors_origins", mode="before")
    @classmethod
    def _split_comma_separated(cls, v: object) -> object:
        """Parse comma-separated env vars into list[str].

        pydantic-settings 2.x otherwise insists on JSON for list-typed env
        values, but ``LAMP_FORGE_API_KEYS=key1,key2`` is what operators
        actually set.
        """
        if isinstance(v, str):
            return [x.strip() for x in v.split(",") if x.strip()]
        return v

    @property
    def resolved_db_url(self) -> str:
        """Return the configured DB URL or the default sqlite file URL."""
        if self.db_url:
            return self.db_url
        db_path = self.data_dir / "jobs.db"
        return f"sqlite+aiosqlite:///{db_path}"

    @property
    def results_dir(self) -> Path:
        """Where per-job result directories live."""
        return self.data_dir / "results"

    @property
    def cache_dir(self) -> Path:
        """Where the NCBI fetch cache lives."""
        return self.data_dir / "cache"

    @property
    def input_dir(self) -> Path:
        """Where per-job upload bundles (off-target FASTAs) are staged."""
        return self.data_dir / "input"

    @property
    def auth_enabled(self) -> bool:
        """Auth is enforced when keys are configured or env == production."""
        return bool(self.api_keys) or self.env == "production"

    def ensure_dirs(self) -> None:
        """Create all data directories. Safe to call at startup."""
        for d in (self.data_dir, self.results_dir, self.cache_dir, self.input_dir):
            d.mkdir(parents=True, exist_ok=True)


@lru_cache(maxsize=1)
def get_settings() -> Settings:
    """Return a cached Settings instance.

    Cached at process scope; tests use ``get_settings.cache_clear()`` plus a
    monkeypatched env to swap config.
    """
    return Settings()  # type: ignore[call-arg]

"""Pydantic request/response schemas for the HTTP API.

These are the public-facing types. The internal :class:`PipelineConfig`
dataclass (from :mod:`lamp_forge.types`) is a separate concern — the schema
here is what consumers of the JSON API see, and we keep it stable across
internal refactors via this translation layer.
"""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, ConfigDict, EmailStr, Field, field_validator


class TargetSpec(BaseModel):
    """How the user identifies the organism / gene to design against."""

    name: str = Field(..., min_length=1, max_length=200, description="Human-readable label.")
    taxon_id: int | None = Field(
        None, ge=1, description="NCBI taxonomy ID. Required if accessions is empty."
    )
    accessions: list[str] = Field(
        default_factory=list,
        max_length=200,
        description="Explicit NCBI accession.version list. Overrides taxon_id.",
    )
    gene: str | None = Field(None, max_length=50, description="Gene name filter, e.g. rpoB.")
    max_sequences: int = Field(20, ge=1, le=200)
    email: EmailStr = Field(..., description="NCBI courtesy email; required by Entrez.")

    @field_validator("accessions")
    @classmethod
    def _strip_accessions(cls, v: list[str]) -> list[str]:
        return [a.strip() for a in v if a.strip()]


class OffTargetsSpec(BaseModel):
    """Off-target screening parameters."""

    min_identity_threshold: float = Field(0.80, ge=0.0, le=1.0)
    min_coverage_threshold: float = Field(0.80, ge=0.0, le=1.0)
    # off-target FASTAs are uploaded separately and referenced by name on the job
    fasta_files: list[str] = Field(
        default_factory=list,
        description="Names of previously-uploaded off-target FASTA files for this job.",
    )


class ConservationSpec(BaseModel):
    window_size: int = Field(30, ge=10, le=100)
    entropy_threshold: float = Field(0.25, ge=0.0, le=2.0)
    min_region_length: int = Field(250, ge=120, le=2000)


class AmpliconSizeSpec(BaseModel):
    f2_b2_min: int = Field(120, ge=80, le=500)
    f2_b2_max: int = Field(160, ge=80, le=500)


class PrimerSpec(BaseModel):
    tm_min: float = Field(60.0, ge=40.0, le=80.0)
    tm_max: float = Field(65.0, ge=40.0, le=80.0)
    tm_match_tolerance: float = Field(2.0, ge=0.0, le=10.0)
    gc_min: float = Field(40.0, ge=0.0, le=100.0)
    gc_max: float = Field(65.0, ge=0.0, le=100.0)
    hairpin_dg_threshold: float = Field(-2.0)
    dimer_dg_threshold: float = Field(-5.0)
    amplicon_size: AmpliconSizeSpec = Field(default_factory=AmpliconSizeSpec)


class OutputSpec(BaseModel):
    top_n: int = Field(10, ge=1, le=100)
    generate_html: bool = True
    generate_csv: bool = True


class JobCreate(BaseModel):
    """Body of POST /api/v1/jobs."""

    model_config = ConfigDict(extra="forbid")

    target: TargetSpec
    off_targets: OffTargetsSpec = Field(default_factory=OffTargetsSpec)
    conservation: ConservationSpec = Field(default_factory=ConservationSpec)
    primer: PrimerSpec = Field(default_factory=PrimerSpec)
    output: OutputSpec = Field(default_factory=OutputSpec)

    @field_validator("primer")
    @classmethod
    def _check_tm_band(cls, v: PrimerSpec) -> PrimerSpec:
        if v.tm_min >= v.tm_max:
            raise ValueError(f"tm_min ({v.tm_min}) must be < tm_max ({v.tm_max})")
        if v.gc_min >= v.gc_max:
            raise ValueError(f"gc_min ({v.gc_min}) must be < gc_max ({v.gc_max})")
        if v.amplicon_size.f2_b2_min >= v.amplicon_size.f2_b2_max:
            raise ValueError("amplicon_size.f2_b2_min must be < f2_b2_max")
        return v


class ProgressEvent(BaseModel):
    """A single event in a job's progress timeline."""

    ts: datetime
    stage: str
    message: str
    percent: float | None = None


class JobSummary(BaseModel):
    """Compact job representation for list endpoints."""

    id: str
    target_name: str
    status: str
    created_at: datetime
    updated_at: datetime
    started_at: datetime | None
    finished_at: datetime | None
    n_primer_sets: int | None
    n_conserved_regions: int | None


class JobDetail(JobSummary):
    """Full job representation including config and progress timeline."""

    config: dict[str, Any]
    progress: list[ProgressEvent]
    error: str | None
    result_urls: dict[str, str] = Field(
        default_factory=dict,
        description=(
            "Map of artefact name → absolute URL. Populated only when status "
            "is `succeeded`."
        ),
    )


class JobCreateResponse(BaseModel):
    """Response from POST /api/v1/jobs."""

    id: str
    status: str
    detail_url: str


class FileUploadResponse(BaseModel):
    """Response from POST /api/v1/uploads."""

    upload_id: str
    filename: str
    size_bytes: int


class HealthResponse(BaseModel):
    """Liveness response."""

    status: str = "ok"
    version: str
    env: str


class ReadinessResponse(BaseModel):
    """Readiness response — reports per-dependency status."""

    status: str
    redis: str
    database: str
    mafft: str
    blast: str


class ErrorResponse(BaseModel):
    """Standard error envelope."""

    detail: str
    code: str | None = None

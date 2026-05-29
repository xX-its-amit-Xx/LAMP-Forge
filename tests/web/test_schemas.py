"""Pydantic schema validation tests."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from lamp_forge.web.schemas import JobCreate, PrimerSpec, TargetSpec


def test_target_requires_email() -> None:
    with pytest.raises(ValidationError):
        TargetSpec.model_validate({"name": "t", "taxon_id": 1, "email": "not-an-email"})


def test_primer_rejects_inverted_tm() -> None:
    with pytest.raises(ValidationError):
        JobCreate.model_validate(
            {
                "target": {"name": "t", "taxon_id": 1, "email": "a@b.com"},
                "primer": {"tm_min": 70.0, "tm_max": 60.0},
            }
        )


def test_primer_rejects_inverted_gc() -> None:
    with pytest.raises(ValidationError):
        JobCreate.model_validate(
            {
                "target": {"name": "t", "taxon_id": 1, "email": "a@b.com"},
                "primer": {"gc_min": 80.0, "gc_max": 40.0},
            }
        )


def test_primer_rejects_inverted_amplicon() -> None:
    with pytest.raises(ValidationError):
        JobCreate.model_validate(
            {
                "target": {"name": "t", "taxon_id": 1, "email": "a@b.com"},
                "primer": {"amplicon_size": {"f2_b2_min": 200, "f2_b2_max": 100}},
            }
        )


def test_extra_fields_rejected() -> None:
    with pytest.raises(ValidationError):
        JobCreate.model_validate(
            {
                "target": {"name": "t", "taxon_id": 1, "email": "a@b.com"},
                "evil_extra_field": "drop tables",
            }
        )


def test_defaults_fill_in() -> None:
    job = JobCreate.model_validate(
        {"target": {"name": "t", "taxon_id": 1, "email": "a@b.com"}}
    )
    assert job.primer.tm_min == 60.0
    assert job.primer.tm_max == 65.0
    assert job.conservation.min_region_length == 250

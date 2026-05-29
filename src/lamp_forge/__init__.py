"""LAMP-Forge: reproducible LAMP primer design for microbial diagnostics.

Public re-exports keep the import surface small. Most users will only ever
touch the CLI; this matters for the notebook walkthrough and downstream
pipelines that import LAMP-Forge as a library.
"""

from lamp_forge.types import (
    ConservedRegion,
    LampPrimerSet,
    PipelineConfig,
    Primer,
    SpecificityHit,
)

__version__ = "0.1.0"

__all__ = [
    "ConservedRegion",
    "LampPrimerSet",
    "PipelineConfig",
    "Primer",
    "SpecificityHit",
    "__version__",
]

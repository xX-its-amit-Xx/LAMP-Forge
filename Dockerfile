# LAMP-Forge: reproducible LAMP primer design pipeline.
#
# Multi-stage build:
#   - "base" installs system + python deps once; reused by all runtime images
#   - "cli"  → ENTRYPOINT lamp-forge (the original single-shot path)
#   - "web"  → ENTRYPOINT lamp-forge-web (FastAPI server)
#   - "worker" → ENTRYPOINT lamp-forge-worker (arq queue consumer)
#
# Build a specific target with:
#   docker build --target web -t lamp-forge:web .

FROM python:3.11-slim-bookworm AS base

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    DEBIAN_FRONTEND=noninteractive

# System deps: MAFFT for MSA, BLAST+ for off-target screening, build tools
# for primer3 native extension, curl for healthchecks.
RUN apt-get update && apt-get install -y --no-install-recommends \
        mafft \
        ncbi-blast+ \
        build-essential \
        curl \
        ca-certificates \
        git \
        tini \
    && rm -rf /var/lib/apt/lists/*

# Non-root user for the runtime stages — never run prod containers as root.
RUN useradd --uid 1000 --create-home --shell /usr/sbin/nologin forge

WORKDIR /opt/lamp-forge

# Install Python deps separately for better layer caching.
COPY pyproject.toml README.md ./
COPY src ./src
RUN pip install --upgrade pip && pip install -e ".[dev,web]"

# Copy the rest of the project tree.
COPY Snakefile ./
COPY config ./config
COPY tests ./tests
COPY docs ./docs

# Common data layout.
RUN mkdir -p /work/data/results /work/data/cache /work/data/input /work/results /work/cache /work/input \
    && chown -R forge:forge /work /opt/lamp-forge

# Sanity check that core CLI binaries resolve.
RUN mafft --version 2>&1 | head -1 \
    && blastn -version | head -1 \
    && lamp-forge --help > /dev/null

WORKDIR /work
USER forge

# =====================================================================
# CLI image — the original single-shot pipeline runner.
# =====================================================================
FROM base AS cli
ENTRYPOINT ["/usr/bin/tini", "--", "lamp-forge"]
CMD ["--help"]

# =====================================================================
# Web image — FastAPI server.
# =====================================================================
FROM base AS web
ENV LAMP_FORGE_DATA_DIR=/work/data
EXPOSE 8000
HEALTHCHECK --interval=30s --timeout=5s --start-period=20s --retries=3 \
    CMD curl -fsS http://localhost:8000/api/v1/health || exit 1
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["lamp-forge-web", "--host", "0.0.0.0", "--port", "8000"]

# =====================================================================
# Worker image — arq queue consumer.
# =====================================================================
FROM base AS worker
ENV LAMP_FORGE_DATA_DIR=/work/data
HEALTHCHECK --interval=30s --timeout=5s --start-period=10s --retries=3 \
    CMD pgrep -f "lamp_forge.web.worker" > /dev/null || exit 1
ENTRYPOINT ["/usr/bin/tini", "--"]
CMD ["lamp-forge-worker"]

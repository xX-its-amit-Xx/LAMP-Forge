"""LAMP-Forge web service.

The `web` subpackage layers a production-grade HTTP interface on top of the
core pipeline. The split is deliberate: the CLI / library API works without
any of the web dependencies, so users who don't need a server pay nothing
for the runtime weight.

Components:
- :mod:`lamp_forge.web.app` — FastAPI application factory.
- :mod:`lamp_forge.web.worker` — arq async worker that runs pipeline jobs.
- :mod:`lamp_forge.web.db` — SQLite job-state persistence (async via aiosqlite).
- :mod:`lamp_forge.web.routes` — HTTP route modules (jobs, results, ui, health).
"""

from lamp_forge.web.app import create_app

__all__ = ["create_app"]

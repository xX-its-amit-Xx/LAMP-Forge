"""``python -m lamp_forge.web`` / ``lamp-forge-web`` entrypoint.

Boots a uvicorn server with the LAMP-Forge FastAPI app. Production
deployments will typically run uvicorn (or gunicorn+uvicorn workers)
directly with their own config, but this is the smallest defensible
local-dev path.
"""

from __future__ import annotations

import argparse

import uvicorn

from lamp_forge.web.config import get_settings


def main() -> None:
    """CLI entrypoint for the web server."""
    parser = argparse.ArgumentParser(prog="lamp-forge-web")
    parser.add_argument("--host", default="0.0.0.0")
    parser.add_argument("--port", type=int, default=8000)
    parser.add_argument("--reload", action="store_true", help="Auto-reload on edits.")
    parser.add_argument("--workers", type=int, default=1)
    args = parser.parse_args()

    s = get_settings()
    uvicorn.run(
        "lamp_forge.web.app:create_app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        workers=args.workers,
        factory=True,
        log_level=s.log_level.lower(),
        access_log=s.debug,
    )


if __name__ == "__main__":
    main()

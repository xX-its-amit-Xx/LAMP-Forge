"""API key authentication.

A simple, defensible scheme for a single-tenant scientific tool:
- Clients pass ``Authorization: Bearer <key>`` (or ``X-API-Key: <key>``).
- The server checks against a whitelist loaded from env at startup.
- Admin endpoints require the key to additionally appear in
  ``admin_api_keys``.

We hash the key on every check to derive a stable, non-secret ``api_key_id``
used in logs and the database — this means we never store or log the raw key.
For larger deployments, swap this module for OAuth/OIDC; the dependency
contract (returning an ``ApiKeyPrincipal``) is preserved.
"""

from __future__ import annotations

import hashlib
import hmac
from dataclasses import dataclass

from fastapi import Depends, Header, HTTPException, status

from lamp_forge.web.config import Settings, get_settings


@dataclass(frozen=True, slots=True)
class ApiKeyPrincipal:
    """The authenticated identity for an HTTP request.

    ``key_id`` is a deterministic hash of the raw API key — safe to log.
    ``is_admin`` is derived from whether the key appears in admin_api_keys.
    """

    key_id: str
    is_admin: bool


def _hash_key(raw: str) -> str:
    """Stable 16-hex-char hash of an API key. Used as the public key id."""
    return hashlib.sha256(raw.encode()).hexdigest()[:16]


def _extract_token(authorization: str | None, x_api_key: str | None) -> str | None:
    """Pull a raw API key from either common header form, or None."""
    if x_api_key:
        return x_api_key.strip()
    if authorization:
        parts = authorization.strip().split(None, 1)
        if len(parts) == 2 and parts[0].lower() == "bearer":
            return parts[1].strip()
        # Allow bare-token Authorization headers as a developer convenience.
        if len(parts) == 1:
            return parts[0].strip()
    return None


def _match(raw: str, accepted: list[str]) -> bool:
    """Constant-time comparison of raw key against accepted list."""
    return any(hmac.compare_digest(raw, a) for a in accepted)


async def require_api_key(
    authorization: str | None = Header(default=None),
    x_api_key: str | None = Header(default=None, alias="X-API-Key"),
    settings: Settings = Depends(get_settings),
) -> ApiKeyPrincipal:
    """FastAPI dependency that validates an API key.

    Behaviour:
    - If ``settings.auth_enabled`` is False (dev mode with no keys
      configured), returns a synthetic "anonymous" principal so endpoints
      remain callable.
    - Otherwise requires a valid key; rejects with 401 if missing/invalid.
    """
    if not settings.auth_enabled:
        return ApiKeyPrincipal(key_id="anonymous", is_admin=False)

    raw = _extract_token(authorization, x_api_key)
    if not raw:
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Missing API key. Pass Authorization: Bearer <key> or X-API-Key.",
        )
    if not _match(raw, settings.api_keys):
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED, detail="Invalid API key."
        )
    return ApiKeyPrincipal(
        key_id=_hash_key(raw), is_admin=_match(raw, settings.admin_api_keys)
    )


async def require_admin(
    principal: ApiKeyPrincipal = Depends(require_api_key),
) -> ApiKeyPrincipal:
    """Dependency that requires the caller to be an admin."""
    if not principal.is_admin:
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN, detail="Admin privileges required."
        )
    return principal

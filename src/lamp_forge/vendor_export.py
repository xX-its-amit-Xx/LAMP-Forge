"""Vendor-ready order CSV export for designed LAMP primer sets.

After ``lamp-forge run`` produces a ``primer_sets.json``, use this module (or
``lamp-forge export``) to convert the top-ranked sets into a spreadsheet that
can be pasted directly into an IDT or Twist bulk-order portal.

IDT format
----------
The IDT SDS / Custom DNA Oligos bulk-input table expects four columns::

    Name | Sequence | Scale | Purification

Scale options (IDT): ``25nm``, ``100nm``, ``250nm``, ``1umol``, ``5umol``, …
Purification options: ``STD``, ``HPLC``, ``PAGE``.

Role-based defaults applied when no override is given:

* F3, B3, LF, LB — ``25nm / STD``  (short outer and loop primers)
* FIP, BIP       — ``100nm / HPLC`` (chimeric ~40-50 nt; HPLC removes
                    incomplete-synthesis products that can cause non-specific
                    amplification in LAMP reactions)

Twist format
------------
Twist Bioscience's custom-oligo CSV requires two columns::

    Name | Sequence

LAMP-Forge appends four informational columns (Scale, Purification,
Tm_celsius, Length) so ordering parameters travel with the file when a user
opens it in a spreadsheet before uploading to the Twist portal.

References:
    IDT Bulk Input: https://www.idtdna.com/pages/tools/bulk-input
    Twist Custom Oligos: https://www.twistbioscience.com/products/oligonucleotides
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from enum import StrEnum
from pathlib import Path
from typing import Any

from lamp_forge.types import LampPrimerSet, Primer


class VendorFormat(StrEnum):
    """Supported vendor order-sheet formats."""

    IDT = "idt"
    TWIST = "twist"


# Per-role synthesis defaults.  FIP/BIP are chimeric oligos; HPLC removes
# failed-synthesis products that silence one arm and kill amplification.
_DEFAULT_SCALE: dict[str, str] = {
    "F3": "25nm",
    "B3": "25nm",
    "FIP": "100nm",
    "BIP": "100nm",
    "LF": "25nm",
    "LB": "25nm",
}
_DEFAULT_PURIFICATION: dict[str, str] = {
    "F3": "STD",
    "B3": "STD",
    "FIP": "HPLC",
    "BIP": "HPLC",
    "LF": "STD",
    "LB": "STD",
}


@dataclass(frozen=True, slots=True)
class VendorRow:
    """One primer expressed as a single vendor-order row.

    Attributes:
        name: Unique primer name for the order portal.
        sequence: Oligo sequence (5' -> 3'), uppercase DNA alphabet.
        scale: Synthesis scale string (e.g. ``"25nm"``, ``"100nm"``).
        purification: Purification level (e.g. ``"STD"``, ``"HPLC"``).
        role: LAMP primer role (``"F3"``, ``"FIP"``, etc.).
        set_id: Identifier of the parent primer set.
        tm_celsius: Melting temperature in degrees Celsius.
        length: Primer length in nucleotides.
    """

    name: str
    sequence: str
    scale: str
    purification: str
    role: str
    set_id: str
    tm_celsius: float
    length: int


# ---------------------------------------------------------------------------
# Naming
# ---------------------------------------------------------------------------


def _primer_name(set_id: str, role: str, target_label: str | None) -> str:
    """Construct a unique, informative order-portal name for one primer.

    Args:
        set_id: Parent set identifier (e.g. ``"region_01_set_001"``).
        role: LAMP primer role string (e.g. ``"FIP"``).
        target_label: Optional short target prefix (e.g. ``"dsrB_SRB"``).
            Spaces are replaced with underscores.

    Returns:
        A name such as ``"dsrB_SRB_region_01_set_001_FIP"``.
    """
    if target_label:
        prefix = target_label.replace(" ", "_")
        return f"{prefix}_{set_id}_{role}"
    return f"{set_id}_{role}"


# ---------------------------------------------------------------------------
# Core conversion helpers
# ---------------------------------------------------------------------------


def primer_to_vendor_row(
    primer: Primer,
    primer_set: LampPrimerSet,
    *,
    target_label: str | None = None,
    scale_override: str | None = None,
    purification_override: str | None = None,
) -> VendorRow:
    """Convert one :class:`~lamp_forge.types.Primer` to a :class:`VendorRow`.

    Args:
        primer: The primer to convert.
        primer_set: The set this primer belongs to (provides ``set_id``).
        target_label: Short target label prepended to the primer name.
        scale_override: If given, use this synthesis scale for all roles
            instead of the role-specific default.
        purification_override: If given, use this purification level for all
            roles instead of the role-specific default.

    Returns:
        :class:`VendorRow` ready to write to a vendor order sheet.
    """
    role = primer.role
    return VendorRow(
        name=_primer_name(primer_set.set_id, role, target_label),
        sequence=primer.sequence,
        scale=scale_override or _DEFAULT_SCALE[role],
        purification=purification_override or _DEFAULT_PURIFICATION[role],
        role=role,
        set_id=primer_set.set_id,
        tm_celsius=primer.tm,
        length=primer.length,
    )


def export_rows(
    primer_sets: list[LampPrimerSet],
    *,
    target_label: str | None = None,
    top_n: int = 1,
    scale_override: str | None = None,
    purification_override: str | None = None,
) -> list[VendorRow]:
    """Build ordered rows for the top *top_n* primer sets.

    Primers within each set are written in canonical LAMP order:
    F3 → B3 → FIP → BIP → LF → LB (loop primers omitted if not designed).

    Args:
        primer_sets: Ranked list, highest composite score first.
        target_label: Short target label prepended to primer names
            (e.g. ``"dsrB_SRB"``).
        top_n: How many sets to include. Default ``1`` exports only the
            best-scoring set, which is the typical pre-order workflow.
        scale_override: Override synthesis scale for all primers.
        purification_override: Override purification level for all primers.

    Returns:
        List of :class:`VendorRow` objects, sets in rank order, primers in
        canonical LAMP order within each set.
    """
    rows: list[VendorRow] = []
    for ps in primer_sets[:top_n]:
        for p in ps.primers:
            rows.append(
                primer_to_vendor_row(
                    p,
                    ps,
                    target_label=target_label,
                    scale_override=scale_override,
                    purification_override=purification_override,
                )
            )
    return rows


def rows_from_json_data(
    sets_data: list[dict[str, Any]],
    *,
    target_label: str | None = None,
    top_n: int = 1,
    scale_override: str | None = None,
    purification_override: str | None = None,
) -> list[VendorRow]:
    """Build :class:`VendorRow` objects directly from JSON-deserialized dicts.

    This avoids reconstructing full :class:`~lamp_forge.types.LampPrimerSet`
    dataclasses when the caller has loaded ``primer_sets.json`` from disk.
    Only the fields needed for ordering are extracted.

    Args:
        sets_data: List of primer-set dicts as returned by ``json.load()`` on
            a ``primer_sets.json`` file (the ``"primer_sets"`` array).
        target_label: Short target label prepended to primer names.
        top_n: Number of sets to include (front of list = highest ranked).
        scale_override: Override synthesis scale for all primers.
        purification_override: Override purification level for all primers.

    Returns:
        List of :class:`VendorRow`, one per primer.
    """
    rows: list[VendorRow] = []
    for set_dict in sets_data[:top_n]:
        set_id = str(set_dict["set_id"])
        for role_key in ("f3", "b3", "fip", "bip", "lf", "lb"):
            primer_dict = set_dict.get(role_key)
            if primer_dict is None:
                continue
            role = str(primer_dict["role"])
            sequence = str(primer_dict["sequence"])
            tm = float(primer_dict["tm"])
            rows.append(
                VendorRow(
                    name=_primer_name(set_id, role, target_label),
                    sequence=sequence,
                    scale=scale_override or _DEFAULT_SCALE.get(role, "100nm"),
                    purification=purification_override or _DEFAULT_PURIFICATION.get(role, "HPLC"),
                    role=role,
                    set_id=set_id,
                    tm_celsius=tm,
                    length=len(sequence),
                )
            )
    return rows


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------


def write_idt_csv(rows: list[VendorRow], path: Path) -> None:
    """Write an IDT SDS / Custom DNA Oligos bulk-input CSV.

    Column order follows the IDT Bulk Input specification:
    ``Name | Sequence | Scale | Purification``.

    Args:
        rows: Vendor rows to serialise.
        path: Destination path. Parent directories are created as needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Name", "Sequence", "Scale", "Purification"])
        for r in rows:
            writer.writerow([r.name, r.sequence, r.scale, r.purification])


def write_twist_csv(rows: list[VendorRow], path: Path) -> None:
    """Write a Twist Bioscience oligo-order CSV.

    The two columns required by the Twist upload form are ``Name`` and
    ``Sequence``. Four informational columns (Scale, Purification,
    Tm_celsius, Length) are appended so synthesis conditions travel with
    the file.

    Args:
        rows: Vendor rows to serialise.
        path: Destination path. Parent directories are created as needed.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Name", "Sequence", "Scale", "Purification", "Tm_celsius", "Length"])
        for r in rows:
            writer.writerow(
                [
                    r.name,
                    r.sequence,
                    r.scale,
                    r.purification,
                    f"{r.tm_celsius:.1f}",
                    r.length,
                ]
            )


def write_vendor_csv(
    rows: list[VendorRow],
    path: Path,
    fmt: VendorFormat,
) -> None:
    """Dispatch to the correct vendor writer based on *fmt*.

    Args:
        rows: Vendor rows to serialise.
        path: Destination path.
        fmt: :class:`VendorFormat` selecting the output schema.
    """
    if fmt is VendorFormat.IDT:
        write_idt_csv(rows, path)
    else:
        write_twist_csv(rows, path)

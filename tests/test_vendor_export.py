"""Tests for the IDT / Twist vendor order-CSV export.

All tests are pure Python — no external binaries, no network.
Primer and LampPrimerSet objects are constructed in-process using the same
frozen dataclasses the pipeline uses, so the export round-trips correctly.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from lamp_forge.vendor_export import (
    VendorFormat,
    VendorRow,
    export_rows,
    primer_to_vendor_row,
    rows_from_json_data,
    write_idt_csv,
    write_twist_csv,
    write_vendor_csv,
)

# ---------------------------------------------------------------------------
# Fixtures — minimal but valid LampPrimerSet / Primer objects
# ---------------------------------------------------------------------------


def _make_primer(
    role: str,
    sequence: str,
    tm: float = 62.0,
) -> object:
    """Build a minimal Primer suitable for export tests."""
    from lamp_forge.types import Primer

    return Primer(
        role=role,  # type: ignore[arg-type]
        sequence=sequence,
        start=0,
        end=len(sequence),
        strand="+" if role in ("F3", "FIP", "LF") else "-",
        tm=tm,
        gc_percent=50.0,
        hairpin_dg=-1.0,
        homodimer_dg=-2.0,
    )


def _make_set(set_id: str = "region_01_set_001") -> object:
    """Build a minimal LampPrimerSet with all six primers."""
    from lamp_forge.types import LampPrimerSet

    f3 = _make_primer("F3", "ACGTACGTACGTACGT", tm=61.5)
    b3 = _make_primer("B3", "TGCATGCATGCATGCA", tm=61.5)
    fip = _make_primer("FIP", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT", tm=63.0)
    bip = _make_primer("BIP", "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA", tm=63.0)
    lf = _make_primer("LF", "GGCCGGCCGGCC", tm=60.5)
    lb = _make_primer("LB", "CCGGCCGGCCGG", tm=60.5)
    return LampPrimerSet(
        set_id=set_id,
        region_id="region_01",
        f3=f3,  # type: ignore[arg-type]
        b3=b3,  # type: ignore[arg-type]
        fip=fip,  # type: ignore[arg-type]
        bip=bip,  # type: ignore[arg-type]
        lf=lf,  # type: ignore[arg-type]
        lb=lb,  # type: ignore[arg-type]
        composite_score=0.9,
    )


# ---------------------------------------------------------------------------
# VendorRow construction
# ---------------------------------------------------------------------------


class TestPrimerToVendorRow:
    def test_basic_fields(self) -> None:

        ps = _make_set()
        f3 = ps.f3  # type: ignore[attr-defined]
        row = primer_to_vendor_row(f3, ps)  # type: ignore[arg-type]
        assert row.role == "F3"
        assert row.sequence == f3.sequence
        assert row.set_id == ps.set_id  # type: ignore[attr-defined]

    def test_default_scale_outer_primer(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(ps.f3, ps)  # type: ignore[attr-defined, arg-type]
        assert row.scale == "25nm"
        assert row.purification == "STD"

    def test_default_scale_chimeric_primer(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(ps.fip, ps)  # type: ignore[attr-defined, arg-type]
        assert row.scale == "100nm"
        assert row.purification == "HPLC"

    def test_scale_override(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(
            ps.f3,
            ps,
            scale_override="250nm",  # type: ignore[attr-defined, arg-type]
        )
        assert row.scale == "250nm"

    def test_purification_override(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(
            ps.fip,
            ps,
            purification_override="PAGE",  # type: ignore[attr-defined, arg-type]
        )
        assert row.purification == "PAGE"

    def test_target_label_included_in_name(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(
            ps.f3,
            ps,
            target_label="dsrB SRB",  # type: ignore[attr-defined, arg-type]
        )
        assert row.name.startswith("dsrB_SRB_")

    def test_spaces_in_label_replaced(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(
            ps.f3,
            ps,
            target_label="my target",  # type: ignore[attr-defined, arg-type]
        )
        assert " " not in row.name

    def test_length_matches_sequence(self) -> None:
        ps = _make_set()
        fip = ps.fip  # type: ignore[attr-defined]
        row = primer_to_vendor_row(fip, ps)  # type: ignore[arg-type]
        assert row.length == len(fip.sequence)

    def test_tm_stored_correctly(self) -> None:
        ps = _make_set()
        row = primer_to_vendor_row(ps.f3, ps)  # type: ignore[attr-defined, arg-type]
        assert abs(row.tm_celsius - 61.5) < 0.01


# ---------------------------------------------------------------------------
# export_rows
# ---------------------------------------------------------------------------


class TestExportRows:
    def test_top1_exports_one_set(self) -> None:
        ps1 = _make_set("set_001")
        ps2 = _make_set("set_002")
        rows = export_rows([ps1, ps2], top_n=1)  # type: ignore[arg-type, list-item]
        set_ids = {r.set_id for r in rows}
        assert set_ids == {"set_001"}

    def test_top2_exports_two_sets(self) -> None:
        ps1 = _make_set("set_001")
        ps2 = _make_set("set_002")
        rows = export_rows([ps1, ps2], top_n=2)  # type: ignore[arg-type, list-item]
        set_ids = {r.set_id for r in rows}
        assert set_ids == {"set_001", "set_002"}

    def test_canonical_role_order_within_set(self) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        roles = [r.role for r in rows]
        assert roles == ["F3", "B3", "FIP", "BIP", "LF", "LB"]

    def test_six_primers_per_set_with_loops(self) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        assert len(rows) == 6

    def test_four_primers_per_set_without_loops(self) -> None:
        from lamp_forge.types import LampPrimerSet

        ps_base = _make_set()
        ps_no_loops = LampPrimerSet(
            set_id="no_loops",
            region_id="region_01",
            f3=ps_base.f3,  # type: ignore[attr-defined]
            b3=ps_base.b3,  # type: ignore[attr-defined]
            fip=ps_base.fip,  # type: ignore[attr-defined]
            bip=ps_base.bip,  # type: ignore[attr-defined]
        )
        rows = export_rows([ps_no_loops])
        assert len(rows) == 4

    def test_empty_input_gives_empty_rows(self) -> None:
        rows = export_rows([])
        assert rows == []

    def test_all_rows_are_vendor_row_instances(self) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        assert all(isinstance(r, VendorRow) for r in rows)


# ---------------------------------------------------------------------------
# write_idt_csv
# ---------------------------------------------------------------------------


class TestWriteIdtCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "idt.csv"
        write_idt_csv(rows, out)
        assert out.exists()

    def test_header_matches_idt_spec(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "idt.csv"
        write_idt_csv(rows, out)
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert header == ["Name", "Sequence", "Scale", "Purification"]

    def test_row_count_matches_primers(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "idt.csv"
        write_idt_csv(rows, out)
        with out.open() as fh:
            all_rows = list(csv.reader(fh))
        assert len(all_rows) == 1 + 6  # header + 6 primers

    def test_sequences_are_uppercase_dna(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "idt.csv"
        write_idt_csv(rows, out)
        with out.open() as fh:
            data_rows = list(csv.DictReader(fh))
        for row in data_rows:
            seq = row["Sequence"]
            assert seq == seq.upper()
            assert all(c in "ACGT" for c in seq)

    def test_creates_parent_directories(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "sub" / "dir" / "idt.csv"
        write_idt_csv(rows, out)
        assert out.exists()

    def test_fip_bip_get_hplc_purification(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "idt.csv"
        write_idt_csv(rows, out)
        with out.open() as fh:
            data_rows = list(csv.DictReader(fh))
        chimeric = {
            r["Purification"]
            for r in data_rows
            if r["Name"].endswith("FIP") or r["Name"].endswith("BIP")
        }
        assert chimeric == {"HPLC"}

    def test_outer_primers_get_std_purification(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "idt.csv"
        write_idt_csv(rows, out)
        with out.open() as fh:
            data_rows = list(csv.DictReader(fh))
        outer = {
            r["Purification"]
            for r in data_rows
            if r["Name"].endswith("F3") or r["Name"].endswith("B3")
        }
        assert outer == {"STD"}


# ---------------------------------------------------------------------------
# write_twist_csv
# ---------------------------------------------------------------------------


class TestWriteTwistCsv:
    def test_creates_file(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "twist.csv"
        write_twist_csv(rows, out)
        assert out.exists()

    def test_header_starts_with_name_sequence(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "twist.csv"
        write_twist_csv(rows, out)
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert header[:2] == ["Name", "Sequence"]

    def test_has_informational_columns(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "twist.csv"
        write_twist_csv(rows, out)
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert "Tm_celsius" in header
        assert "Length" in header

    def test_tm_celsius_is_numeric(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "twist.csv"
        write_twist_csv(rows, out)
        with out.open() as fh:
            data_rows = list(csv.DictReader(fh))
        for row in data_rows:
            float(row["Tm_celsius"])  # should not raise

    def test_length_column_matches_sequence(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "twist.csv"
        write_twist_csv(rows, out)
        with out.open() as fh:
            data_rows = list(csv.DictReader(fh))
        for row in data_rows:
            assert int(row["Length"]) == len(row["Sequence"])


# ---------------------------------------------------------------------------
# write_vendor_csv dispatcher
# ---------------------------------------------------------------------------


class TestWriteVendorCsv:
    def test_idt_dispatch_produces_four_column_header(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "out.csv"
        write_vendor_csv(rows, out, VendorFormat.IDT)
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert header == ["Name", "Sequence", "Scale", "Purification"]

    def test_twist_dispatch_produces_six_column_header(self, tmp_path: Path) -> None:
        ps = _make_set()
        rows = export_rows([ps])  # type: ignore[arg-type, list-item]
        out = tmp_path / "out.csv"
        write_vendor_csv(rows, out, VendorFormat.TWIST)
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert len(header) == 6


# ---------------------------------------------------------------------------
# rows_from_json_data
# ---------------------------------------------------------------------------


class TestRowsFromJsonData:
    def _json_set(self, set_id: str = "set_001") -> dict[str, object]:
        """Minimal JSON representation matching primer_sets.json schema."""
        return {
            "set_id": set_id,
            "region_id": "region_01",
            "f3": {"role": "F3", "sequence": "ACGTACGT", "tm": 61.5},
            "b3": {"role": "B3", "sequence": "TGCATGCA", "tm": 61.5},
            "fip": {"role": "FIP", "sequence": "ACGTACGTACGTACGTACGTACGT", "tm": 63.0},
            "bip": {"role": "BIP", "sequence": "TGCATGCATGCATGCATGCATGCA", "tm": 63.0},
            "lf": {"role": "LF", "sequence": "GGCCGGCC", "tm": 60.5},
            "lb": {"role": "LB", "sequence": "CCGGCCGG", "tm": 60.5},
            "composite_score": 0.9,
        }

    def test_parses_all_six_primers(self) -> None:
        rows = rows_from_json_data([self._json_set()])
        assert len(rows) == 6

    def test_role_order(self) -> None:
        rows = rows_from_json_data([self._json_set()])
        assert [r.role for r in rows] == ["F3", "B3", "FIP", "BIP", "LF", "LB"]

    def test_sequence_preserved(self) -> None:
        rows = rows_from_json_data([self._json_set()])
        f3_row = next(r for r in rows if r.role == "F3")
        assert f3_row.sequence == "ACGTACGT"

    def test_tm_preserved(self) -> None:
        rows = rows_from_json_data([self._json_set()])
        fip_row = next(r for r in rows if r.role == "FIP")
        assert abs(fip_row.tm_celsius - 63.0) < 0.01

    def test_top_n_limits_sets(self) -> None:
        sets = [self._json_set("set_001"), self._json_set("set_002")]
        rows = rows_from_json_data(sets, top_n=1)
        assert all(r.set_id == "set_001" for r in rows)

    def test_missing_loop_primers_skipped(self) -> None:
        s = self._json_set()
        s.pop("lf")  # type: ignore[union-attr]
        s.pop("lb")  # type: ignore[union-attr]
        rows = rows_from_json_data([s])
        assert len(rows) == 4
        assert all(r.role not in ("LF", "LB") for r in rows)

    def test_target_label_in_names(self) -> None:
        rows = rows_from_json_data([self._json_set()], target_label="myTarget")
        assert all(r.name.startswith("myTarget_") for r in rows)

    def test_scale_override_applied(self) -> None:
        rows = rows_from_json_data([self._json_set()], scale_override="250nm")
        assert all(r.scale == "250nm" for r in rows)

    def test_purification_override_applied(self) -> None:
        rows = rows_from_json_data([self._json_set()], purification_override="PAGE")
        assert all(r.purification == "PAGE" for r in rows)

    def test_empty_input_gives_empty_rows(self) -> None:
        rows = rows_from_json_data([])
        assert rows == []

    def test_length_computed_from_sequence(self) -> None:
        rows = rows_from_json_data([self._json_set()])
        for r in rows:
            assert r.length == len(r.sequence)


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------


class TestExportCli:
    def _write_primer_sets_json(self, path: Path) -> None:
        """Write a minimal primer_sets.json fixture."""
        payload = {
            "version": "0.1.0",
            "generated_at": "2026-01-01T00:00:00+00:00",
            "primer_sets": [
                {
                    "set_id": "region_01_set_001",
                    "region_id": "region_01",
                    "f3": {"role": "F3", "sequence": "ACGTACGTACGTACGT", "tm": 61.5},
                    "b3": {"role": "B3", "sequence": "TGCATGCATGCATGCA", "tm": 61.5},
                    "fip": {
                        "role": "FIP",
                        "sequence": "ACGTACGTACGTACGTACGTACGTACGTACGTACGT",
                        "tm": 63.0,
                    },
                    "bip": {
                        "role": "BIP",
                        "sequence": "TGCATGCATGCATGCATGCATGCATGCATGCATGCA",
                        "tm": 63.0,
                    },
                    "lf": {"role": "LF", "sequence": "GGCCGGCCGGCC", "tm": 60.5},
                    "lb": {"role": "LB", "sequence": "CCGGCCGGCCGG", "tm": 60.5},
                    "composite_score": 0.9,
                }
            ],
        }
        path.write_text(json.dumps(payload), encoding="utf-8")

    def test_idt_default_exit_zero(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        out = tmp_path / "order.csv"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["export", "--input", str(input_json), "--out", str(out)],
            obj={},
        )
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_idt_output_has_correct_header(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        out = tmp_path / "order.csv"
        runner = CliRunner()
        runner.invoke(
            cli,
            ["export", "--input", str(input_json), "--out", str(out), "--format", "idt"],
            obj={},
        )
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert header == ["Name", "Sequence", "Scale", "Purification"]

    def test_twist_output_has_name_sequence_first(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        out = tmp_path / "order.csv"
        runner = CliRunner()
        runner.invoke(
            cli,
            ["export", "--input", str(input_json), "--out", str(out), "--format", "twist"],
            obj={},
        )
        with out.open() as fh:
            header = next(csv.reader(fh))
        assert header[:2] == ["Name", "Sequence"]

    def test_default_output_path_created(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["export", "--input", str(input_json)],
            obj={},
        )
        assert result.exit_code == 0, result.output
        default_out = tmp_path / "primer_sets_idt_order.csv"
        assert default_out.exists()

    def test_target_label_appears_in_names(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        out = tmp_path / "order.csv"
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "export",
                "--input",
                str(input_json),
                "--out",
                str(out),
                "--target-label",
                "dsrB_SRB",
            ],
            obj={},
        )
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        assert all(r["Name"].startswith("dsrB_SRB_") for r in rows)

    def test_output_mentions_primer_count(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["export", "--input", str(input_json)],
            obj={},
        )
        assert "primers" in result.output.lower() or "6" in result.output

    def test_scale_override_propagates(self, tmp_path: Path) -> None:
        from click.testing import CliRunner

        from lamp_forge.cli import cli

        input_json = tmp_path / "primer_sets.json"
        self._write_primer_sets_json(input_json)
        out = tmp_path / "order.csv"
        runner = CliRunner()
        runner.invoke(
            cli,
            [
                "export",
                "--input",
                str(input_json),
                "--out",
                str(out),
                "--scale",
                "250nm",
            ],
            obj={},
        )
        with out.open() as fh:
            rows = list(csv.DictReader(fh))
        assert all(r["Scale"] == "250nm" for r in rows)

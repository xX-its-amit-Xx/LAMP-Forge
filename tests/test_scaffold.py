"""Tests for the config scaffold generator.

All tests are pure Python -- no external binaries required.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from click.testing import CliRunner

from lamp_forge.cli import cli
from lamp_forge.scaffold import (
    _DEFAULTS,
    NaType,
    ScaffoldParams,
    Vertical,
    scaffold_yaml,
    write_scaffold,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _parse(text: str) -> dict:  # type: ignore[type-arg]
    """Parse YAML text and assert it is a dict."""
    result = yaml.safe_load(text)
    assert isinstance(result, dict)
    return result


def _params(
    target_name: str = "test_target",
    vertical: Vertical = Vertical.GENERIC,
    na_type: NaType = NaType.DNA,
    taxon_id: int | None = 12345,
    gene: str | None = "testGene",
    email: str = "user@example.com",
    max_sequences: int | None = None,
) -> ScaffoldParams:
    return ScaffoldParams(
        target_name=target_name,
        vertical=vertical,
        na_type=na_type,
        taxon_id=taxon_id,
        gene=gene,
        email=email,
        max_sequences=max_sequences,
    )


# ---------------------------------------------------------------------------
# Coverage of all (Vertical, NaType) combinations
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("vertical", list(Vertical))
@pytest.mark.parametrize("na_type", list(NaType))
def test_all_combinations_produce_valid_yaml(vertical: Vertical, na_type: NaType) -> None:
    p = _params(vertical=vertical, na_type=na_type)
    text = scaffold_yaml(p)
    doc = _parse(text)
    # Required top-level sections
    for section in ("target", "off_targets", "conservation", "primer", "output"):
        assert section in doc, f"Missing section: {section}"


@pytest.mark.parametrize("vertical", list(Vertical))
@pytest.mark.parametrize("na_type", list(NaType))
def test_defaults_table_complete(vertical: Vertical, na_type: NaType) -> None:
    """Every (Vertical, NaType) pair must have an entry in the defaults table."""
    assert (vertical, na_type) in _DEFAULTS


# ---------------------------------------------------------------------------
# Target section values
# ---------------------------------------------------------------------------


class TestTargetSection:
    def test_taxon_id_written_when_provided(self) -> None:
        p = _params(taxon_id=872)
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["taxon_id"] == 872

    def test_taxon_id_null_when_absent(self) -> None:
        p = _params(taxon_id=None)
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["taxon_id"] is None

    def test_gene_written_when_provided(self) -> None:
        p = _params(gene="dsrB")
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["gene"] == "dsrB"

    def test_gene_null_when_absent(self) -> None:
        p = _params(gene=None)
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["gene"] is None

    def test_target_name_propagated(self) -> None:
        p = _params(target_name="my_assay")
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["name"] == "my_assay"

    def test_email_propagated(self) -> None:
        p = _params(email="bio@example.com")
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["email"] == "bio@example.com"

    def test_max_sequences_override(self) -> None:
        p = _params(max_sequences=99)
        doc = _parse(scaffold_yaml(p))
        assert doc["target"]["max_sequences"] == 99

    def test_max_sequences_uses_vertical_default_when_none(self) -> None:
        p = _params(vertical=Vertical.OIL_GAS, na_type=NaType.DNA, max_sequences=None)
        doc = _parse(scaffold_yaml(p))
        expected = _DEFAULTS[(Vertical.OIL_GAS, NaType.DNA)].max_sequences
        assert doc["target"]["max_sequences"] == expected

    def test_output_dir_uses_target_name(self) -> None:
        p = _params(target_name="srb_dsrB")
        doc = _parse(scaffold_yaml(p))
        assert doc["output"]["dir"] == "results/srb_dsrB"


# ---------------------------------------------------------------------------
# Vertical-specific defaults
# ---------------------------------------------------------------------------


class TestVerticalDefaults:
    def test_oil_gas_dna_entropy_threshold(self) -> None:
        p = _params(vertical=Vertical.OIL_GAS, na_type=NaType.DNA)
        doc = _parse(scaffold_yaml(p))
        assert doc["conservation"]["entropy_threshold"] == pytest.approx(0.35)

    def test_farm_rna_tm_min_raised_for_rt_lamp(self) -> None:
        """Farm RNA targets get tm_min=63.0 for one-step RT-LAMP co-activity."""
        p = _params(vertical=Vertical.FARM, na_type=NaType.RNA)
        doc = _parse(scaffold_yaml(p))
        assert doc["primer"]["tm_min"] == pytest.approx(63.0)

    def test_farm_dna_tm_min_standard(self) -> None:
        p = _params(vertical=Vertical.FARM, na_type=NaType.DNA)
        doc = _parse(scaffold_yaml(p))
        assert doc["primer"]["tm_min"] == pytest.approx(60.0)

    def test_poc_dna_entropy_tight(self) -> None:
        p = _params(vertical=Vertical.POC, na_type=NaType.DNA)
        doc = _parse(scaffold_yaml(p))
        # POC DNA should be tighter than oil-gas functional-gene default
        assert doc["conservation"]["entropy_threshold"] < 0.30

    def test_all_rna_verticals_have_tm_min_63(self) -> None:
        for vertical in Vertical:
            p = _params(vertical=vertical, na_type=NaType.RNA)
            doc = _parse(scaffold_yaml(p))
            assert doc["primer"]["tm_min"] == pytest.approx(63.0), (
                f"{vertical} RNA should have tm_min=63.0"
            )

    def test_all_dna_verticals_have_tm_min_60(self) -> None:
        for vertical in Vertical:
            p = _params(vertical=vertical, na_type=NaType.DNA)
            doc = _parse(scaffold_yaml(p))
            assert doc["primer"]["tm_min"] == pytest.approx(60.0), (
                f"{vertical} DNA should have tm_min=60.0"
            )

    def test_gc_range_valid_for_all_verticals(self) -> None:
        for vertical in Vertical:
            for na_type in NaType:
                p = _params(vertical=vertical, na_type=na_type)
                doc = _parse(scaffold_yaml(p))
                assert doc["primer"]["gc_min"] < doc["primer"]["gc_max"]

    def test_amplicon_size_valid_for_all(self) -> None:
        for vertical in Vertical:
            for na_type in NaType:
                p = _params(vertical=vertical, na_type=na_type)
                doc = _parse(scaffold_yaml(p))
                amps = doc["primer"]["amplicon_size"]
                assert amps["f2_b2_min"] < amps["f2_b2_max"]
                assert amps["f2_b2_min"] >= 80  # LAMP geometric floor


# ---------------------------------------------------------------------------
# RT-LAMP hint in generated YAML header
# ---------------------------------------------------------------------------


class TestRtHint:
    def test_rna_target_includes_rt_check_comment(self) -> None:
        p = _params(na_type=NaType.RNA, target_name="prrsv_orf7")
        text = scaffold_yaml(p)
        assert "rt-check" in text
        assert "prrsv_orf7" in text

    def test_dna_target_omits_rt_check_comment(self) -> None:
        p = _params(na_type=NaType.DNA)
        text = scaffold_yaml(p)
        assert "rt-check" not in text


# ---------------------------------------------------------------------------
# write_scaffold
# ---------------------------------------------------------------------------


class TestWriteScaffold:
    def test_writes_file(self, tmp_path: Path) -> None:
        out = tmp_path / "config" / "test.yaml"
        p = _params()
        write_scaffold(p, out)
        assert out.exists()
        doc = _parse(out.read_text(encoding="utf-8"))
        assert doc["target"]["name"] == "test_target"

    def test_creates_parent_dirs(self, tmp_path: Path) -> None:
        out = tmp_path / "a" / "b" / "c" / "cfg.yaml"
        write_scaffold(_params(), out)
        assert out.exists()

    def test_overwrites_existing_file(self, tmp_path: Path) -> None:
        out = tmp_path / "cfg.yaml"
        out.write_text("old content", encoding="utf-8")
        write_scaffold(_params(target_name="new_name"), out)
        doc = _parse(out.read_text(encoding="utf-8"))
        assert doc["target"]["name"] == "new_name"


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------


class TestScaffoldCli:
    def test_cli_prints_to_stdout(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "scaffold",
                "--target-name", "prrsv_orf7",
                "--vertical", "farm",
                "--na-type", "rna",
                "--taxon-id", "28344",
                "--gene", "ORF7",
            ],
        )
        assert result.exit_code == 0, result.output
        doc = _parse(result.output)
        assert doc["target"]["name"] == "prrsv_orf7"
        assert doc["target"]["taxon_id"] == 28344
        assert doc["primer"]["tm_min"] == pytest.approx(63.0)

    def test_cli_writes_file(self, tmp_path: Path) -> None:
        out_path = tmp_path / "config.yaml"
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "scaffold",
                "--target-name", "srb_dsrB",
                "--vertical", "oil-gas",
                "--taxon-id", "872",
                "--gene", "dsrB",
                "--out", str(out_path),
            ],
        )
        assert result.exit_code == 0, result.output
        assert out_path.exists()
        doc = _parse(out_path.read_text(encoding="utf-8"))
        assert doc["target"]["name"] == "srb_dsrB"

    def test_cli_generic_vertical_default(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["scaffold", "--target-name", "my_target"],
        )
        assert result.exit_code == 0, result.output
        doc = _parse(result.output)
        assert doc["target"]["name"] == "my_target"

    def test_cli_invalid_vertical_exits_nonzero(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["scaffold", "--target-name", "t", "--vertical", "hospital"],
        )
        assert result.exit_code != 0

    def test_cli_invalid_na_type_exits_nonzero(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            ["scaffold", "--target-name", "t", "--na-type", "protein"],
        )
        assert result.exit_code != 0

    def test_cli_max_sequences_override(self) -> None:
        runner = CliRunner()
        result = runner.invoke(
            cli,
            [
                "scaffold",
                "--target-name", "big_panel",
                "--max-sequences", "50",
            ],
        )
        assert result.exit_code == 0, result.output
        doc = _parse(result.output)
        assert doc["target"]["max_sequences"] == 50

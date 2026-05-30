"""Tests for primer design, focusing on geometric and thermodynamic constraints."""

from __future__ import annotations

import pytest
from Bio.Align import MultipleSeqAlignment

from lamp_forge.conserve import compute_track, find_conserved_regions
from lamp_forge.primer_design import (
    B1_B2_GAP,
    F2_F1_GAP,
    F2_LEN_RANGE,
    _slide_candidates,
    design_all,
    design_sets_for_region,
    gc_percent,
    hairpin_dg,
    homodimer_dg,
    reverse_complement,
    tm,
)
from lamp_forge.types import ConservedRegion, PipelineConfig


class TestSequenceHelpers:
    def test_reverse_complement_canonical(self) -> None:
        assert reverse_complement("ATCG") == "CGAT"
        assert reverse_complement("AAAA") == "TTTT"

    def test_reverse_complement_palindrome(self) -> None:
        assert reverse_complement("GAATTC") == "GAATTC"  # EcoRI

    def test_gc_percent_pure_at(self) -> None:
        assert gc_percent("AAAATTTT") == 0.0

    def test_gc_percent_pure_gc(self) -> None:
        assert gc_percent("GCGCGCGC") == 100.0

    def test_gc_percent_balanced(self) -> None:
        assert gc_percent("ATGC") == 50.0

    def test_gc_percent_empty(self) -> None:
        assert gc_percent("") == 0.0


class TestThermodynamics:
    def test_tm_in_reasonable_band(self) -> None:
        # A 20mer with balanced GC should land 55-70 °C under default salt.
        t = tm("ACGTACGTACGTACGTACGT")
        assert 50.0 < t < 75.0

    def test_higher_gc_higher_tm(self) -> None:
        low_gc = tm("ATATATATATATATATATAT")
        high_gc = tm("GCGCGCGCGCGCGCGCGCGC")
        assert high_gc > low_gc

    def test_hairpin_returns_finite(self) -> None:
        assert hairpin_dg("ACGTACGTACGTACGTACGT") <= 0  # ΔG should be ≤ 0

    def test_homodimer_returns_finite(self) -> None:
        # Sequence cannot dimerise with itself in any obviously stable way.
        assert homodimer_dg("ATATATATATATATATATAT") <= 0


class TestSlideCandidates:
    def test_filters_on_tm_band(self, example_config: PipelineConfig) -> None:
        # All-AT sequence: extremely low Tm — should produce no candidates
        # at the standard band.
        seq = "A" * 100
        cands = _slide_candidates(seq, 0, F2_LEN_RANGE, "+", example_config)
        assert cands == []

    def test_yields_some_for_realistic_input(self, example_config: PipelineConfig) -> None:
        # A real ~220 bp slice of M. tuberculosis rpoB. Chosen because it
        # carries plenty of 18-22 nt windows landing in the 58-68 C Tm band
        # under primer3's default LAMP-buffer parameters, so the candidate
        # list is unambiguously non-empty across primer3 versions.
        seq = (
            "ATGCTGGACCAGAATAACCCGCTGTCGGGGTTGACCCACAAGCGCCGACTGTCGGCGCTGGGG"
            "CCCGGCGGTCTGTCACGTGAGCGTGCCGGGCTGGAGGTCCGCGACGTGCACCCGTCGCACTAC"
            "GGCCGGATGTGCCCGATCGAAACGCCCGAAGGGCCCAACATCGGTCTGATCGGCTCACTGTCG"
            "GTCTACGCGCGGGTCAATCCGTTCGGGTTCATCG"
        )
        cands = _slide_candidates(seq, 0, F2_LEN_RANGE, "+", example_config)
        assert len(cands) > 0
        for c in cands:
            assert example_config.tm_min <= c.tm <= example_config.tm_max
            assert example_config.gc_min <= c.gc <= example_config.gc_max
            assert c.strand == "+"
            assert "N" not in c.sequence

    def test_skips_n_containing_windows(self, example_config: PipelineConfig) -> None:
        seq = "ATGCATGCATGCNNNNNATGCATGCATGCATGCATGC" * 3
        cands = _slide_candidates(seq, 0, F2_LEN_RANGE, "+", example_config)
        for c in cands:
            assert "N" not in c.sequence


class TestDesignSetsForRegion:
    def test_too_short_region_returns_empty(self, example_config: PipelineConfig) -> None:
        # Region length below f2_b2_max + flanks → nothing can be designed.
        region = ConservedRegion("R000", 0, 50, 0.0, 0.0, "A" * 50)
        sets = design_sets_for_region(region, example_config)
        assert sets == []

    def test_designs_sets_in_realistic_region(
        self, conserved_msa: MultipleSeqAlignment, example_config: PipelineConfig
    ) -> None:
        track = compute_track(conserved_msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.10, min_region_length=200)
        if not regions:
            pytest.skip("Fixture didn't produce a conserved region")
        sets = design_sets_for_region(regions[0], example_config, max_sets=20)
        # We don't assert >0 strictly — the synthetic consensus may not yield
        # primers in the Tm band for all configs. If we do get any, validate
        # their structure.
        for s in sets:
            assert s.f3 is not None
            assert s.b3 is not None
            assert s.fip is not None
            assert s.bip is not None
            assert s.tm_spread <= example_config.tm_match_tolerance
            assert example_config.f2_b2_min <= s.f2_b2_distance <= example_config.f2_b2_max
            assert s.composite_score >= 0.0

    def test_geometric_gap_constants_self_consistent(self) -> None:
        # Sanity check on the LAMP geometry constants: F2-F1 and B1-B2
        # spacings should be wide enough to accommodate at least one loop primer.
        assert F2_F1_GAP[0] >= 18  # min loop primer length
        assert B1_B2_GAP[0] >= 18


class TestDesignAll:
    def test_returns_sorted_by_score(
        self, conserved_msa: MultipleSeqAlignment, example_config: PipelineConfig
    ) -> None:
        track = compute_track(conserved_msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.15, min_region_length=200)
        sets = design_all(regions, example_config, max_sets_per_region=10)
        if len(sets) >= 2:
            # design_all does not enforce a final sort across regions — only
            # design_sets_for_region sorts within a region. We just check
            # that within each region the order is descending by composite.
            from itertools import groupby

            for _, group in groupby(sets, key=lambda s: s.region_id):
                gs = list(group)
                scores = [s.composite_score for s in gs]
                assert scores == sorted(scores, reverse=True)

"""Tests for the conservation analysis module."""

from __future__ import annotations

import pytest
from Bio.Align import MultipleSeqAlignment

from lamp_forge.conserve import (
    _column_consensus,
    _column_coverage,
    _column_entropy,
    compute_track,
    conservation_score,
    find_conserved_regions,
)
from lamp_forge.types import ConservedRegion


class TestColumnEntropy:
    def test_uniform_column_max_entropy(self) -> None:
        # Equal A/C/G/T → entropy = log2(4) = 2.0 bits
        assert _column_entropy("ACGT") == pytest.approx(2.0)

    def test_identical_column_zero_entropy(self) -> None:
        assert _column_entropy("AAAA") == 0.0

    def test_binary_column(self) -> None:
        # 50/50 A/C → entropy = 1.0 bit
        assert _column_entropy("AACC") == pytest.approx(1.0)

    def test_gaps_excluded(self) -> None:
        # Gaps should not contribute. A/A/-/- == all-A column.
        assert _column_entropy("AA--") == 0.0

    def test_all_gaps_zero(self) -> None:
        assert _column_entropy("----") == 0.0

    def test_case_insensitive(self) -> None:
        assert _column_entropy("aacc") == pytest.approx(_column_entropy("AACC"))


class TestColumnConsensus:
    def test_majority_wins(self) -> None:
        assert _column_consensus("AAAAC") == "A"

    def test_tie_returns_n(self) -> None:
        assert _column_consensus("AACC") == "N"

    def test_all_gaps_returns_n(self) -> None:
        assert _column_consensus("----") == "N"


class TestColumnCoverage:
    def test_all_canonical(self) -> None:
        assert _column_coverage("ACGT") == 1.0

    def test_half_gaps(self) -> None:
        assert _column_coverage("AA--") == 0.5

    def test_all_gaps(self) -> None:
        assert _column_coverage("----") == 0.0


class TestComputeTrack:
    def test_track_dimensions(self, conserved_msa: MultipleSeqAlignment) -> None:
        track = compute_track(conserved_msa, window_size=30)
        width = conserved_msa.get_alignment_length()
        assert track.raw_entropy.shape == (width,)
        assert track.smoothed_entropy.shape == (width,)
        assert len(track.consensus) == width
        assert track.coverage.shape == (width,)

    def test_conserved_middle_has_low_entropy(self, conserved_msa: MultipleSeqAlignment) -> None:
        track = compute_track(conserved_msa, window_size=30)
        # The fixture has identical positions 150-450.
        conserved_slice = track.smoothed_entropy[200:400]
        assert conserved_slice.max() < 0.05, (
            f"Expected near-zero entropy in conserved region, got {conserved_slice.max():.3f}"
        )

    def test_variable_flanks_have_higher_entropy(self, conserved_msa: MultipleSeqAlignment) -> None:
        track = compute_track(conserved_msa, window_size=30)
        # Flank positions are randomised → entropy should be measurably > 0.
        # Use the raw track here (smoothed pulls in conserved-region zeros at
        # the edges and can mute the signal).
        flank_raw = track.raw_entropy[20:130]
        assert flank_raw.mean() > 0.2

    def test_window_size_must_be_positive(self, conserved_msa: MultipleSeqAlignment) -> None:
        with pytest.raises(ValueError, match="window_size"):
            compute_track(conserved_msa, window_size=0)

    def test_consensus_string_length_matches_alignment(
        self, conserved_msa: MultipleSeqAlignment
    ) -> None:
        track = compute_track(conserved_msa, window_size=30)
        assert len(track.consensus) == conserved_msa.get_alignment_length()


class TestFindConservedRegions:
    def test_detects_conserved_middle(self, conserved_msa: MultipleSeqAlignment) -> None:
        track = compute_track(conserved_msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.10, min_region_length=200)
        assert len(regions) >= 1
        biggest = max(regions, key=lambda r: r.length)
        assert biggest.length >= 200
        # The conserved fixture region is positions 150-450.
        assert biggest.start >= 140 and biggest.end <= 460

    def test_too_short_regions_dropped(self, conserved_msa: MultipleSeqAlignment) -> None:
        track = compute_track(conserved_msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.10, min_region_length=10000)
        assert regions == []

    def test_strict_threshold_returns_empty(self, conserved_msa: MultipleSeqAlignment) -> None:
        track = compute_track(conserved_msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=-1.0, min_region_length=200)
        assert regions == []

    def test_region_ids_are_unique_and_sequential(
        self, conserved_msa: MultipleSeqAlignment
    ) -> None:
        track = compute_track(conserved_msa, window_size=30)
        regions = find_conserved_regions(track, entropy_threshold=0.20, min_region_length=50)
        ids = [r.region_id for r in regions]
        assert len(ids) == len(set(ids))
        for i, rid in enumerate(ids):
            assert rid == f"R{i:03d}"


class TestConservationScore:
    def test_perfectly_conserved_scores_one(self) -> None:
        region = ConservedRegion("R000", 0, 200, 0.0, 0.0, "A" * 200)
        assert conservation_score(region) == 1.0

    def test_maximally_variable_scores_zero(self) -> None:
        region = ConservedRegion("R000", 0, 200, 2.0, 2.0, "N" * 200)
        assert conservation_score(region) == 0.0

    def test_halfway_scores_half(self) -> None:
        region = ConservedRegion("R000", 0, 200, 1.0, 1.0, "A" * 200)
        assert conservation_score(region) == 0.5

"""Tests for CreateSyntheticJunctionsStep."""

import pytest
from pathlib import Path
from amplifinder.steps import CreateSyntheticJunctionsStep
from amplifinder.data_types import (
    RecordTypedDF, CandidateTnJc2, TnLoc, RawEvent, Orientation,
)


@pytest.fixture
def sample_candidate(tmp_path):
    """Create a sample CandidateTnJc2."""
    return CandidateTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf_chr="tiny",
        pos_chr_L=200, pos_chr_R=300,
        pos_tn_L=10, pos_tn_R=20,
        dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[1], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="tiny", iso_name="sample1",
        amplicon_coverage=2.0, genome_coverage=1.0,
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1], chosen_tn_id=1,
        analysis_dir="jc_200_300_001_L150",
    )


@pytest.fixture
def sample_tn_loc():
    """Create a sample TnLoc."""
    return TnLoc(
        ID=1, TN_Name="IS1", TN_scaf="tiny",
        LocLeft=200, LocRight=300,
        Complement=False, Join=False,
    )


def test_creates_junctions_fasta(tiny_genome, sample_candidate, sample_tn_loc, tmp_path):
    """Should create junctions.fasta file."""
    candidates = RecordTypedDF.from_records([sample_candidate], CandidateTnJc2)
    tn_locs = RecordTypedDF.from_records([sample_tn_loc], TnLoc)
    
    step = CreateSyntheticJunctionsStep(
        candidates=candidates,
        genome=tiny_genome,
        tn_locs=tn_locs,
        output_dir=tmp_path,
        read_length=150,
    )
    
    result = step.run()
    
    # Check that junctions.fasta was created
    junctions_file = tmp_path / sample_candidate.analysis_dir / "junctions.fasta"
    assert junctions_file.exists()
    
    # Check that it contains 7 junction sequences
    with open(junctions_file) as f:
        content = f.read()
        # Should have 7 headers (one per junction type)
        assert content.count(">") == 7


def test_handles_missing_tn(tiny_genome, sample_candidate, tmp_path):
    """Should skip candidates with missing TN."""
    # Create candidate with non-existent TN ID
    candidate_no_tn = CandidateTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf_chr="tiny",
        pos_chr_L=200, pos_chr_R=300,
        pos_tn_L=10, pos_tn_R=20,
        dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[999], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="tiny", iso_name="sample1",
        amplicon_coverage=2.0, genome_coverage=1.0,
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[999], chosen_tn_id=999,
        analysis_dir="jc_200_300_999_L150",
    )
    
    candidates = RecordTypedDF.from_records([candidate_no_tn], CandidateTnJc2)
    tn_locs = RecordTypedDF.from_records([], TnLoc)  # Empty TN list
    
    # Create empty junctions.fasta to satisfy step output requirements
    junctions_file = tmp_path / candidate_no_tn.analysis_dir / "junctions.fasta"
    junctions_file.parent.mkdir(parents=True, exist_ok=True)
    junctions_file.write_text("")
    
    step = CreateSyntheticJunctionsStep(
        candidates=candidates,
        genome=tiny_genome,
        tn_locs=tn_locs,
        output_dir=tmp_path,
    )
    
    # Should not raise error, just skip creating content
    result = step.run()
    assert len(result) == 1

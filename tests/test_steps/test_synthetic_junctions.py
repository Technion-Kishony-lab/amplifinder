"""Tests for CreateSyntheticJunctionsStep."""

import pytest
from pathlib import Path
from amplifinder.steps import CreateSyntheticJunctionsStep
from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, RefTnLoc, RawEvent, Orientation, Coverage,
)
from amplifinder.utils.file_utils import ensure_parent_dir


@pytest.fixture
def sample_candidate(tmp_path):
    """Create a sample FilteredTnJc2."""
    return FilteredTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf="tiny",
        pos_scaf_L=200, pos_scaf_R=300,
        pos_tn_L=10, pos_tn_R=20,
        dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[1], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="tiny", iso_name="sample1",
        amplicon_coverage=2.0, scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1], chosen_tn_id=1,
        analysis_dir="jc_200_300_001_L150",
    )


@pytest.fixture
def sample_tn_loc():
    """Create a sample RefTnLoc."""
    return RefTnLoc(
        tn_id=1, tn_name="IS1", tn_scaf="tiny",
        loc_left=200, loc_right=300,
        complement=False, join=False,
    )


def test_creates_junctions_fasta(tiny_genome, sample_candidate, sample_tn_loc, tmp_path):
    """Should create junctions.fasta file."""
    filtered_tnjc2s = RecordTypedDf.from_records([sample_candidate], FilteredTnJc2)
    tn_locs = RecordTypedDf.from_records([sample_tn_loc], RefTnLoc)
    
    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
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
    candidate_no_tn = FilteredTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf="tiny",
        pos_scaf_L=200, pos_scaf_R=300,
        pos_tn_L=10, pos_tn_R=20,
        dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[999], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="tiny", iso_name="sample1",
        amplicon_coverage=2.0, scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[999], chosen_tn_id=999,
        analysis_dir="jc_200_300_999_L150",
    )
    
    filtered_tnjc2s = RecordTypedDf.from_records([candidate_no_tn], FilteredTnJc2)
    tn_locs = RecordTypedDf.from_records([], RefTnLoc)  # Empty TN list
    
    # Create empty junctions.fasta to satisfy step output requirements
    junctions_file = tmp_path / candidate_no_tn.analysis_dir / "junctions.fasta"
    ensure_parent_dir(junctions_file)
    junctions_file.write_text("")
    
    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        tn_locs=tn_locs,
        output_dir=tmp_path,
    )
    
    # Should not raise error, just skip creating content
    result = step.run()
    assert len(result) == 1

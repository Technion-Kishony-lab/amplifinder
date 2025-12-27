"""Tests for AnalyzeTnJc2AlignmentsStep."""

import pytest
from pathlib import Path
from amplifinder.steps import AnalyzeTnJc2AlignmentsStep
from amplifinder.data_types import RecordTypedDf, FilteredTnJc2, RawEvent, Orientation, Coverage
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


def test_step_initialization(sample_candidate, tmp_path):
    """Should initialize step correctly."""
    filtered_tnjc2s = RecordTypedDf.from_records([sample_candidate], FilteredTnJc2)
    
    # Create dummy BAM file (empty, just for initialization test)
    bam_file = tmp_path / sample_candidate.analysis_dir / "iso.sorted.bam"
    ensure_parent_dir(bam_file)
    bam_file.write_text("dummy")
    
    step = AnalyzeTnJc2AlignmentsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        output_dir=tmp_path,
        read_length=150,
        has_ancestor=False,
    )
    
    assert step.filtered_tnjc2s == filtered_tnjc2s
    assert step.read_length == 150
    assert step.has_ancestor is False

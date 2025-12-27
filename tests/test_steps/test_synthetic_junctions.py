"""Tests for CreateSyntheticJunctionsStep."""

from amplifinder.steps import CreateSyntheticJunctionsStep
from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, RefTnLoc, Orientation,
)
from amplifinder.utils.file_utils import ensure_parent_dir


def test_creates_junctions_fasta(tiny_genome, filtered_tnjc2_record, ref_tn_loc_record, tmp_path):
    """Should create junctions.fasta file."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], FilteredTnJc2)
    tn_locs = RecordTypedDf.from_records([ref_tn_loc_record], RefTnLoc)

    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        tn_locs=tn_locs,
        output_dir=tmp_path,
        read_length=150,
    )

    step.run()

    # Check that junctions.fasta was created
    junctions_file = tmp_path / filtered_tnjc2_record.analysis_dir / "junctions.fasta"
    assert junctions_file.exists()

    # Check that it contains 7 junction sequences
    with open(junctions_file) as f:
        content = f.read()
        # Should have 7 headers (one per junction type)
        assert content.count(">") == 7


def test_handles_missing_tn(tiny_genome, filtered_tnjc2_record, tmp_path):
    """Should skip candidates with missing TN."""
    candidate_no_tn = FilteredTnJc2.from_other(
        filtered_tnjc2_record,
        tn_ids=[999],
        tn_orientations=[Orientation.FORWARD],
        shared_tn_ids=[999],
        chosen_tn_id=999,
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

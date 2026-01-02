"""Tests for CreateSyntheticJunctionsStep."""

from amplifinder.steps import CreateSyntheticJunctionsStep
from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, Orientation,
)
from amplifinder.utils.file_utils import ensure_parent_dir


def test_creates_junctions_fasta(tiny_genome, filtered_tnjc2_record, tmp_path):
    """Should create junctions.fasta file."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], FilteredTnJc2)

    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        output_dir=tmp_path,
        read_length=150,
    )

    step.run()

    # Check that junctions.fasta was created
    junctions_file = tmp_path / "junctions" / filtered_tnjc2_record.analysis_dir / "junctions.fasta"
    assert junctions_file.exists()

    # Check that it contains 7 junction sequences
    with open(junctions_file) as f:
        content = f.read()
        # Should have 7 headers (one per junction type)
        assert content.count(">") == 7


def test_handles_missing_tn(tiny_genome, filtered_tnjc2_record, tmp_path):
    """Should skip candidates with missing TN (chosen_tn_id=None)."""
    candidate_no_tn = FilteredTnJc2.from_other(
        filtered_tnjc2_record,
        chosen_tn_id=None,  # No chosen TN
        analysis_dir="jc_200_300_no_tn",
    )

    filtered_tnjc2s = RecordTypedDf.from_records([candidate_no_tn], FilteredTnJc2)

    # Create empty junctions.fasta to satisfy step output requirements
    junctions_file = tmp_path / "junctions" / candidate_no_tn.analysis_dir / "junctions.fasta"
    ensure_parent_dir(junctions_file)
    junctions_file.write_text("")

    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        output_dir=tmp_path,
    )

    # Should not raise error, just skip creating content
    result = step.run()
    assert len(result) == 1

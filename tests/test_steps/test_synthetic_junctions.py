"""Tests for CreateSyntheticJunctionsStep."""

from amplifinder.steps import CreateSyntheticJunctionsStep
from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2


def test_creates_junctions_fasta(tiny_genome, filtered_tnjc2_record, ref_tns_indexed, tmp_path):
    """Should create junctions.fasta file."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], SynJctsTnJc2)

    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        ref_tns=ref_tns_indexed,
        output_dir=tmp_path,
        junction_length=150,
    )

    result = step.run()

    # Check that result is returned
    assert len(result) == 1
    first_record = result[0]
    assert first_record.analysis_dir is not None

    # Check that junctions.fasta was created
    junctions_file = tmp_path / "junctions" / first_record.analysis_dir / "junctions.fasta"
    assert junctions_file.exists()

    # Check that it contains 7 junction sequences
    with open(junctions_file) as f:
        content = f.read()
        # Should have 7 headers (one per junction type)
        assert content.count(">") == 7

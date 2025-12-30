"""Tests for AlignReadsToJunctionsStep."""

from amplifinder.steps import AlignReadsToJunctionsStep
from amplifinder.data_types import RecordTypedDf, FilteredTnJc2
from amplifinder.utils.file_utils import ensure_parent_dir


def test_step_initialization(filtered_tnjc2_record, tmp_path):
    """Should initialize step correctly."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], FilteredTnJc2)

    # Create dummy FASTQ file
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

    # Create dummy junctions.fasta
    junctions_file = tmp_path / "junctions" / filtered_tnjc2_record.analysis_dir / "junctions.fasta"
    ensure_parent_dir(junctions_file)
    junctions_file.write_text(">1\nACGTACGT\n")

    step = AlignReadsToJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        output_dir=tmp_path,
        iso_fastq_path=fastq_file,
    )

    assert step.filtered_tnjc2s == filtered_tnjc2s
    assert step.iso_fastq_path == fastq_file
    assert step.anc_fastq_path is None

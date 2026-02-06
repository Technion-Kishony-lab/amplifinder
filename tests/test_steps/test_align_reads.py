"""Tests for AlignReadsToJunctionsStep."""

from amplifinder.steps import AlignReadsToJunctionsStep
from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2
from amplifinder.utils.file_utils import ensure_parent_dir


def test_step_initialization(filtered_tnjc2_record, tmp_path):
    """Should initialize step correctly."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], SynJctsTnJc2)

    # Create dummy FASTQ file
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text("@read1\nACGT\n+\nIIII\n")

    # Create dummy junctions.fasta
    junctions_file = tmp_path / "junctions" / filtered_tnjc2_record.analysis_dir / "junctions.fasta"
    ensure_parent_dir(junctions_file)
    junctions_file.write_text(">1\nACGTACGT\n")

    step = AlignReadsToJunctionsStep(
        synjcs_tnjc2s=filtered_tnjc2s,
        output_dir=tmp_path,
        fastq_path=fastq_file,
        read_length=100,
    )

    assert step.synjcs_tnjc2s == filtered_tnjc2s
    assert step.fastq_path == fastq_file
    assert step.is_ancestor is False

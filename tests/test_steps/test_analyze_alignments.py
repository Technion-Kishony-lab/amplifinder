"""Tests for AnalyzeTnJc2AlignmentsStep."""

from amplifinder.steps import AnalyzeTnJc2AlignmentsStep
from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2
from amplifinder.utils.file_utils import ensure_parent_dir


def test_step_initialization(filtered_tnjc2_record, tmp_path):
    """Should initialize step correctly."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], SynJctsTnJc2)

    # Create dummy BAM file (empty, just for initialization test)
    bam_file = tmp_path / "junctions" / filtered_tnjc2_record.analysis_dir / "iso.sorted.bam"
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

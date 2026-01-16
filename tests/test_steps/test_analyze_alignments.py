"""Tests for AnalyzeTnJc2AlignmentsStep."""

from amplifinder.steps import AnalyzeTnJc2AlignmentsStep
from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2
from amplifinder.utils.file_utils import ensure_parent_dir


def test_step_initialization(filtered_tnjc2_record, tiny_genome, tmp_path):
    """Should initialize step correctly."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], SynJctsTnJc2)

    # Create dummy BAM file (empty, just for initialization test)
    bam_file = tmp_path / "junctions" / filtered_tnjc2_record.analysis_dir / "sorted.bam"
    ensure_parent_dir(bam_file)
    bam_file.write_text("dummy")

    # Create dummy breseq path
    iso_breseq_path = tmp_path / "breseq"
    iso_breseq_path.mkdir(exist_ok=True)

    step = AnalyzeTnJc2AlignmentsStep(
        synjct_tnjc2s=filtered_tnjc2s,
        output_dir=tmp_path,
        genome=tiny_genome,
        iso_breseq_path=iso_breseq_path,
        iso_read_length=150,
        anc_output_dir=None,
        anc_breseq_path=None,
        create_plots=False,
    )

    assert step.synjct_tnjc2s == filtered_tnjc2s
    assert step.iso_read_length == 150
    assert step.has_ancestor is False

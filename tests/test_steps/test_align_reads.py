"""Tests for AlignReadsToJunctionsStep."""

import pytest
from pathlib import Path
from amplifinder.steps import AlignReadsToJunctionsStep
from amplifinder.data_types import RecordTypedDf, CandidateTnJc2, RawEvent, Orientation


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


def test_step_initialization(sample_candidate, tmp_path):
    """Should initialize step correctly."""
    candidates = RecordTypedDf.from_records([sample_candidate], CandidateTnJc2)
    
    # Create dummy FASTQ file
    fastq_file = tmp_path / "test.fastq"
    fastq_file.write_text("@read1\nACGT\n+\nIIII\n")
    
    # Create dummy junctions.fasta
    junctions_file = tmp_path / sample_candidate.analysis_dir / "junctions.fasta"
    junctions_file.parent.mkdir(parents=True)
    junctions_file.write_text(">1\nACGTACGT\n")
    
    step = AlignReadsToJunctionsStep(
        candidates=candidates,
        output_dir=tmp_path,
        iso_fastq_path=fastq_file,
    )
    
    assert step.candidates == candidates
    assert step.iso_fastq_path == fastq_file
    assert step.anc_fastq_path is None

"""Tests for Bowtie2 alignment tools."""

import pytest

from amplifinder.tools.bowtie2 import align_reads_to_fasta
from tests.env import RUN_BOWTIE2_TESTS
from tests.conftest import make_random_seq, write_fasta


skip_no_bowtie2 = pytest.mark.skipif(not RUN_BOWTIE2_TESTS, reason="Bowtie2 tests disabled")


def write_fastq(path, name, seq, quality=None):
    """Write single-read FASTQ file."""
    if quality is None:
        quality = "I" * len(seq)  # High quality scores
    with open(path, "w") as f:
        f.write(f"@{name}\n{seq}\n+\n{quality}\n")


@pytest.fixture
def ref_fasta(tmp_path):
    """Create reference FASTA."""
    path = tmp_path / "ref.fasta"
    seq = make_random_seq(500, seed=42)
    write_fasta(path, "test_ref", seq)
    return path, seq


@pytest.fixture
def reads_fastq(tmp_path, ref_fasta):
    """Create FASTQ with reads from reference."""
    path = tmp_path / "reads.fastq"
    _, ref_seq = ref_fasta

    # Extract reads from reference (simulating perfect alignments)
    reads = [
        ("read1", str(ref_seq[50:150])),
        ("read2", str(ref_seq[100:200])),
        ("read3", str(ref_seq[200:300])),
    ]

    with open(path, "w") as f:
        for name, seq in reads:
            f.write(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n")

    return path


@skip_no_bowtie2
def test_align_reads_to_fasta_creates_bam(tmp_path, ref_fasta, reads_fastq):
    """Should create BAM and BAI files."""
    ref_path, _ = ref_fasta
    output_bam = tmp_path / "output" / "alignments.bam"

    align_reads_to_fasta(
        ref_fasta=ref_path,
        fastq_path=reads_fastq,
        output_bam=output_bam,
        threads=1,
    )

    # Check BAM created
    assert output_bam.exists(), "BAM file not created"
    assert output_bam.stat().st_size > 0, "BAM file is empty"

    # Check BAI index created
    bai_path = output_bam.parent / f"{output_bam.name}.bai"
    assert bai_path.exists(), "BAI index not created"


@skip_no_bowtie2
def test_align_reads_to_fasta_cleanup(tmp_path, ref_fasta, reads_fastq):
    """Should cleanup intermediate files (SAM, index)."""
    ref_path, _ = ref_fasta
    output_bam = tmp_path / "output" / "alignments.bam"

    align_reads_to_fasta(
        ref_fasta=ref_path,
        fastq_path=reads_fastq,
        output_bam=output_bam,
        keep_sam=False,
    )

    work_dir = output_bam.parent

    # SAM should be cleaned up
    assert not (work_dir / "alignments.sam").exists(), "SAM not cleaned up"

    # Index files should be cleaned up
    index_files = list(work_dir.glob("bowtie2_index.*"))
    assert len(index_files) == 0, f"Index files not cleaned up: {index_files}"

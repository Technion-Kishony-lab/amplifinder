"""Tests for BreseqStep."""

import pytest
from unittest.mock import patch

from amplifinder.steps import BreseqStep
from tests.conftest import FIXTURES_DIR


@pytest.fixture
def breseq_step(tmp_path, tiny_ref_fasta):
    """Create BreseqStep with mock inputs."""
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")
    
    output = tmp_path / "breseq_out"
    output.mkdir()
    
    return BreseqStep(
        fastq_path=fastq,
        ref_file=tiny_ref_fasta,
        output_path=output,
        docker=False,
        threads=1,
    )


def test_breseq_step_calls_run_breseq(breseq_step):
    """Should call run_breseq with correct params."""
    def mock_run_fn(*args, **kwargs):
        # Create output file when "running"
        (breseq_step.output_path / "output").mkdir(exist_ok=True)
        (breseq_step.output_path / "output" / "output.gd").touch()
    
    with patch("amplifinder.steps.run_breseq.run_breseq", side_effect=mock_run_fn) as mock_run:
        breseq_step.run()
        
        mock_run.assert_called_once_with(
            ref_paths=[breseq_step.ref_file],
            fastq_path=breseq_step.fastq_path,
            output_path=breseq_step.output_path,
            docker=False,
            threads=1,
        )


def test_breseq_step_skips_if_output_exists(breseq_step):
    """Should skip if output.gd already exists."""
    (breseq_step.output_path / "output").mkdir(exist_ok=True)
    (breseq_step.output_path / "output" / "output.gd").touch()
    
    assert breseq_step.run() is False


def test_breseq_step_read_outputs(tmp_path):
    """Should parse breseq output correctly."""
    step = BreseqStep(
        fastq_path=tmp_path / "dummy.fastq",
        ref_file=tmp_path / "dummy.fasta",
        output_path=FIXTURES_DIR / "breseq",
        docker=False,
        threads=1,
    )
    
    # Create dummy inputs
    (tmp_path / "dummy.fastq").touch()
    (tmp_path / "dummy.fasta").touch()
    
    result = step.read_outputs()
    
    assert isinstance(result, dict)
    assert "JC" in result

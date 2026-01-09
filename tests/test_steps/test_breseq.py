"""Tests for BreseqStep."""

import pytest
from unittest.mock import patch

from amplifinder.steps import BreseqStep
from amplifinder.data_types import RecordTypedDf, BreseqJunction
from amplifinder.utils.file_utils import ensure_dir
from tests.conftest import FIXTURES_DIR


@pytest.fixture
def breseq_step(tmp_path, tiny_ref_fasta):
    """Create BreseqStep with mock inputs."""
    fastq = tmp_path / "reads.fastq"
    fastq.write_text("@read1\nACGT\n+\nIIII\n")

    output = ensure_dir(tmp_path / "breseq_out")

    return BreseqStep(
        breseq_path=output,
        fastq_path=fastq,
        ref_file=tiny_ref_fasta,
        docker=False,
        threads=1,
    )


def test_breseq_step_calls_run_breseq(breseq_step):
    """Should call run_breseq with correct params."""
    def mock_run_fn(*args, **kwargs):
        # Create output file when "running"
        ensure_dir(breseq_step.breseq_output_path / "output")
        (breseq_step.breseq_output_path / "output" / "output.gd").touch()

    with patch("amplifinder.steps.run_breseq.run_breseq", side_effect=mock_run_fn) as mock_run:
        breseq_step.run()

        mock_run.assert_called_once_with(
            ref_paths=[breseq_step.ref_file],
            fastq_path=breseq_step.fastq_path,
            output_path=breseq_step.breseq_output_path,
            docker=False,
            threads=1,
        )


def test_breseq_step_skips_if_output_exists(breseq_step):
    """Should skip if output.gd already exists."""
    ensure_dir(breseq_step.breseq_output_path / "output")
    (breseq_step.breseq_output_path / "output" / "output.gd").touch()

    breseq_step.run()
    assert breseq_step.run_count == 0  # skipped


def test_breseq_step_read_outputs():
    """Should parse breseq output correctly when output exists."""
    # When output exists, fastq_path and ref_file are optional
    step = BreseqStep(
        breseq_path=FIXTURES_DIR / "breseq",
    )

    result = step.load_outputs()

    assert isinstance(result, RecordTypedDf)
    assert len(result.df) > 0
    assert all(col in result.df.columns for col in ["num", "scaf1", "pos1", "scaf2", "pos2"])


def test_breseq_step_requires_inputs_when_no_output(tmp_path):
    """Should raise error when output doesn't exist and inputs not provided."""
    output = ensure_dir(tmp_path / "breseq_out")

    with pytest.raises(ValueError, match="fastq_path and ref_file are required"):
        BreseqStep(breseq_path=output)


# =============================================================================
# Integration tests with real breseq output
# =============================================================================

@pytest.mark.integration
class TestBreseqIntegration:
    """Test breseq output parsing using existing breseq results."""

    def test_parse_existing_breseq_output(self, tmp_path, isolate_srr25242877):
        """Parse pre-computed breseq output for SRR25242877."""
        breseq_path = isolate_srr25242877["breseq_path"]

        if not breseq_path.exists():
            pytest.skip(f"Breseq output not found: {breseq_path}")

        # Parse the output.gd file directly
        from amplifinder.tools.breseq import parse_breseq_output

        output_gd = breseq_path / "output" / "output.gd"
        if not output_gd.exists():
            pytest.skip(f"output.gd not found: {output_gd}")

        # parse_breseq_output expects the breseq output directory, not the .gd file
        result = parse_breseq_output(breseq_path)

        # Check JC (junction) output
        assert "JC" in result
        jc_df = result["JC"]
        assert len(jc_df) > 0

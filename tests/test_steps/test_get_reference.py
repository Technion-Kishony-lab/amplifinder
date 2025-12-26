"""Tests for GetRefGenomeStep."""

import pytest

from amplifinder.steps import GetRefGenomeStep
from amplifinder.utils.tools import ensure_dir


@pytest.mark.integration
class TestGetReference:
    """Test reference genome loading for U00096."""

    def test_load_u00096_from_ncbi(self, tmp_path, isolate_srr25242877):
        """Load E. coli U00096 reference from NCBI."""
        ref_path = ensure_dir(tmp_path / "genomesDB")

        step = GetRefGenomeStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        )
        genome = step.run()

        assert genome.name == "U00096"
        assert genome.genbank_path.exists()
        assert genome.fasta_path.exists()
        assert genome.length > 4_000_000  # E. coli ~4.6M bp


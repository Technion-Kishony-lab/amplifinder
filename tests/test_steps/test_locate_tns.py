"""Tests for TN element location steps (Genbank and ISfinder)."""

import pandas as pd
import pytest

from amplifinder.steps import LocateTNsUsingGenbankStep, LocateTNsUsingISfinderStep
from amplifinder.data_types import RecordTypedDF
from tests.env import RUN_BLAST_TESTS


skip_no_blast = pytest.mark.skipif(not RUN_BLAST_TESTS, reason="BLAST tests disabled")


@pytest.fixture(params=[
    pytest.param("genbank", marks=[]),
    pytest.param("isfinder", marks=[skip_no_blast]),
])
def step_factory(request, tmp_output, tiny_genome, tiny_tn_db):
    """Factory to create TN location steps - parametrized for both backends."""
    output_dir = tmp_output / "tn_loc" / tiny_genome.name
    def make_step(force=False):
        if request.param == "genbank":
            return LocateTNsUsingGenbankStep(
                genome=tiny_genome,
                output_dir=output_dir,
                force=force,
            )
        else:
            return LocateTNsUsingISfinderStep(
                genome=tiny_genome,
                output_dir=output_dir,
                isdb_path=tiny_tn_db,
                evalue=1e-4,
                critical_coverage=0.9,
                force=force,
            )
    return make_step


@pytest.fixture
def locate_tns_step(step_factory):
    """Create TN location step with default params."""
    return step_factory()


def test_extracts_tn_elements(locate_tns_step):
    """Should extract TN elements correctly."""
    tn_loc = locate_tns_step.run()

    assert isinstance(tn_loc, RecordTypedDF)

    expected = pd.DataFrame({
        "ID": [1, 2],
        "TN_Name": ["IS_test1", "IS_test2"],
        "TN_scaf": ["tiny", "tiny"],
        "LocLeft": [501, 1601],
        "LocRight": [1200, 2200],
        "Complement": [False, True],
        "Join": [False, False],
    })
    pd.testing.assert_frame_equal(tn_loc.df, expected)


def test_skips_if_output_exists(locate_tns_step):
    """Should skip execution if output already exists."""
    locate_tns_step.run()
    assert locate_tns_step.run_count == 1  # first run
    locate_tns_step.run()
    assert locate_tns_step.run_count == 1  # second run skips


def test_force_reruns(step_factory):
    """Force=True should re-run even if output exists."""
    step_factory().run()
    step = step_factory(force=True)
    step.run()
    assert step.run_count == 1


# =============================================================================
# Integration tests with real NCBI data
# =============================================================================

@pytest.mark.integration
class TestLocateTNsIntegration:
    """Test TN element location using GenBank annotations with real data."""

    @pytest.fixture
    def u00096_genome(self, tmp_path):
        """Load U00096 reference (cached)."""
        from amplifinder.steps import GetRefGenomeStep
        
        ref_path = tmp_path / "genomesDB"
        ref_path.mkdir()

        step = GetRefGenomeStep(
            ref_name="U00096",
            ref_path=ref_path,
            ncbi=True,
        )
        return step.run()

    def test_locate_tns_genbank(self, tmp_path, u00096_genome, isolate_srr25242877):
        """Locate TN elements from GenBank annotations."""
        output_dir = tmp_path / "tn_loc" / u00096_genome.name
        step = LocateTNsUsingGenbankStep(
            genome=u00096_genome,
            output_dir=output_dir,
        )
        tn_loc = step.run()

        # MATLAB found IS elements in U00096
        assert tn_loc is not None
        assert len(tn_loc) > 0

        # Check expected columns
        assert "TN_Name" in tn_loc.df.columns
        assert "LocLeft" in tn_loc.df.columns
        assert "LocRight" in tn_loc.df.columns

    def test_locate_tns_isfinder(self, tmp_path, u00096_genome):
        """Locate TN elements using ISfinder database."""
        from amplifinder.data import get_builtin_isfinder_db_path

        output_dir = tmp_path / "tn_loc" / u00096_genome.name
        step = LocateTNsUsingISfinderStep(
            genome=u00096_genome,
            output_dir=output_dir,
            isdb_path=get_builtin_isfinder_db_path(),
            evalue=1e-4,
            critical_coverage=0.9,
        )
        tn_loc = step.run()

        assert len(tn_loc) > 0

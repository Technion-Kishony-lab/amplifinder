"""Tests for TN element location steps (Genbank and ISfinder)."""

import pytest

from amplifinder.steps import LocateTNsUsingGenbankStep, LocateTNsUsingISfinderStep
from amplifinder.data_types import Orientation, RecordTypedDf
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

    assert isinstance(tn_loc, RecordTypedDf)

    # Check that we found the right number of TNs
    assert len(tn_loc) == 2

    # Check specific fields
    assert list(tn_loc.df["tn_name"]) == ["IS_test1", "IS_test2"]
    assert list(tn_loc.df["scaf"]) == ["tiny", "tiny"]
    # Coordinates from conftest: IS1=902bp, IS5=1195bp
    # IS_test1 (forward): start=501, end=1402
    # IS_test2 (reverse): start=2997 (higher), end=1803 (lower)
    assert tn_loc.df["start"].iloc[0] in [500, 501]  # IS_test1 start (forward)
    assert tn_loc.df["start"].iloc[1] in [2997, 2998]  # IS_test2 start (reverse - higher coord)
    assert tn_loc.df["end"].iloc[0] in [1402, 1403]  # IS_test1 end (forward)
    assert tn_loc.df["end"].iloc[1] in [1803, 1804]  # IS_test2 end (reverse - lower coord)
    assert list(tn_loc.df["orientation"]) == [Orientation.FORWARD, Orientation.REVERSE]
    assert list(tn_loc.df["join"]) == [False, False]


def test_skips_if_output_exists(locate_tns_step, step_factory):
    """Should skip execution if output already exists."""
    locate_tns_step.run()
    assert locate_tns_step._artifacts_generated is True  # first run

    # New instance should skip due to cached artifacts
    step2 = step_factory()
    step2.run()
    assert step2._artifacts_generated is False  # skipped


def test_force_reruns(step_factory):
    """Force=True should re-run even if output exists."""
    step_factory().run()
    step = step_factory(force=True)
    step.run()
    assert step._artifacts_generated is True


# =============================================================================
# Integration tests with real NCBI data
# =============================================================================

@pytest.mark.integration
class TestLocateTNsIntegration:
    """Test TN element location using GenBank annotations with real data."""

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
        assert "tn_name" in tn_loc.df.columns
        assert "start" in tn_loc.df.columns
        assert "end" in tn_loc.df.columns

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

"""Tests for IS element location steps (Genbank and ISfinder)."""

import pandas as pd
import pytest

from amplifinder.steps import LocateISsUsingGenbank, LocateISsUsingISfinder
from tests.env import RUN_BLAST_TESTS


skip_no_blast = pytest.mark.skipif(not RUN_BLAST_TESTS, reason="BLAST tests disabled")


@pytest.fixture(params=[
    pytest.param("genbank", marks=[]),
    pytest.param("isfinder", marks=[skip_no_blast]),
])
def step_factory(request, tmp_output, tiny_ref_gbk, tiny_ref_fasta, tiny_is_db):
    """Factory to create IS location steps - parametrized for both backends."""
    def make_step(force=False):
        if request.param == "genbank":
            return LocateISsUsingGenbank(
                genbank_path=tiny_ref_gbk,
                ref_name="tiny",
                ref_path=tmp_output,
                force=force,
            )
        else:
            return LocateISsUsingISfinder(
                ref_fasta=tiny_ref_fasta,
                ref_name="tiny",
                ref_path=tmp_output,
                isdb_path=tiny_is_db,
                force=force,
            )
    return make_step


@pytest.fixture
def locate_iss_step(step_factory):
    """Create IS location step with default params."""
    return step_factory()


def test_extracts_is_elements(locate_iss_step):
    """Should extract IS elements correctly."""
    is_loc = locate_iss_step.run_and_read_outputs()
    
    expected = pd.DataFrame({
        "ID": [1, 2],
        "IS_Name": ["IS_test1", "IS_test2"],
        "IS_scaf": ["tiny", "tiny"],
        "LocLeft": [501, 1601],
        "LocRight": [1200, 2200],
        "Complement": [False, True],
        "Join": [False, False],
    })
    pd.testing.assert_frame_equal(is_loc, expected)


def test_skips_if_output_exists(locate_iss_step):
    """Should skip execution if output already exists."""
    assert locate_iss_step.run() is True   # first run
    assert locate_iss_step.run() is False  # second run skips


def test_force_reruns(step_factory):
    """Force=True should re-run even if output exists."""
    step_factory().run()
    assert step_factory(force=True).run() is True


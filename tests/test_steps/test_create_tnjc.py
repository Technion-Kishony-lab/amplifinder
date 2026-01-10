"""Tests for CreateTnJcStep."""

import pytest

from amplifinder.steps import CreateRefTnJcStep, CreateTnJcStep
from amplifinder.data_types import RecordTypedDf, Junction


@pytest.fixture
def ref_tnjcs(locate_tns_step_factory, tiny_genome, tmp_output):
    """Create reference TN junctions."""
    tn_loc = locate_tns_step_factory().run()

    return CreateRefTnJcStep(
        ref_tn_locs=tn_loc,
        genome=tiny_genome,
        output_dir=tmp_output,
        source="isfinder",
        reference_IS_out_span=50,
    ).run()


@pytest.fixture
def mock_junctions(locate_tns_step_factory, tiny_genome, tmp_output):
    """Create mock junctions (using ref TN junctions as input)."""
    tn_loc = locate_tns_step_factory().run()

    ref_jc = CreateRefTnJcStep(
        ref_tn_locs=tn_loc,
        genome=tiny_genome,
        output_dir=tmp_output,
        source="isfinder",
        reference_IS_out_span=50,
    ).run()

    # Convert to base Junction type
    junctions = [
        Junction(
            scaf1=jc.scaf1, pos1=jc.pos1, dir1=jc.dir1,
            scaf2=jc.scaf2, pos2=jc.pos2, dir2=jc.dir2,
            flanking1=jc.flanking1, flanking2=jc.flanking2,
        )
        for jc in ref_jc
    ]
    return RecordTypedDf.from_records(junctions, Junction)


@pytest.fixture
def tnjc_step_factory(mock_junctions, ref_tnjcs, tiny_genome, tmp_output):
    """Factory to create new TnJc steps."""

    def make():
        return CreateTnJcStep(
            junctions=mock_junctions,
            ref_tnjcs=ref_tnjcs,
            genome=tiny_genome,
            output_dir=tmp_output,
            max_dist_to_tn=20,
            trim_jc_flanking=5,
        )

    return make


def test_runs_without_error(tnjc_step_factory):
    """Should run and produce output file."""
    tnjc_step = tnjc_step_factory()
    tnjcs = tnjc_step.run()

    assert isinstance(tnjcs, RecordTypedDf)
    assert tnjc_step.output_file.exists()


def test_output_has_correct_columns(tnjc_step_factory):
    """TnJc output should have expected columns."""
    tnjc_step = tnjc_step_factory()
    tnjcs = tnjc_step.run()

    expected_cols = {"num", "scaf1", "pos1", "dir1", "scaf2", "pos2", "dir2",
                     "flanking1", "flanking2", "ref_tn_sides", "swapped"}
    assert expected_cols.issubset(set(tnjcs.df.columns))


def test_skips_if_exists(tnjc_step_factory):
    """Should skip if output exists."""
    step1 = tnjc_step_factory()
    step1.run()
    assert step1._artifacts_generated is True

    # New instance sees cached artifacts and skips generation
    step2 = tnjc_step_factory()
    step2.run()
    assert step2._artifacts_generated is False

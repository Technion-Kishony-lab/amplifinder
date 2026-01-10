"""Tests for CreateRefTnJcStep."""

import pytest

from amplifinder.steps import CreateRefTnJcStep
from amplifinder.data_types import RecordTypedDf, Side


@pytest.fixture
def ref_jc_factory(locate_tns_step, tiny_genome, tmp_output):
    """Factory to create new CreateRefTnJcStep instances."""
    tn_loc = locate_tns_step.run()

    def make():
        return CreateRefTnJcStep(
            ref_tn_locs=tn_loc,
            genome=tiny_genome,
            output_dir=tmp_output,
            source="isfinder",
            reference_IS_out_span=50,
        )

    return make


def test_creates_junction_records(ref_jc_factory):
    """Should create 2 junctions per TN (start + end)."""
    ref_jc = ref_jc_factory().run()

    assert isinstance(ref_jc, RecordTypedDf)
    assert len(ref_jc.df) == 4  # 2 TNs * 2 sides
    assert set(jc.ref_tn_side.side for jc in ref_jc) == {Side.START, Side.END}


def test_junction_positions_correct(ref_jc_factory, locate_tns_step_factory):
    """Junction positions should match TN boundaries."""
    ref_jc = ref_jc_factory().run()
    tn_loc = locate_tns_step_factory().run()  # run() returns output

    tn1 = tn_loc.df.iloc[0]
    start_jc = next(jc for jc in ref_jc if jc.ref_tn_side.tn_id == tn1["tn_id"] and jc.ref_tn_side.side == Side.START)
    end_jc = next(jc for jc in ref_jc if jc.ref_tn_side.tn_id == tn1["tn_id"] and jc.ref_tn_side.side == Side.END)

    assert start_jc.pos1 == tn1["start"]
    assert end_jc.pos1 == tn1["end"]


def test_skips_if_exists(ref_jc_factory):
    """Should skip if output exists."""
    step1 = ref_jc_factory()
    step1.run()
    assert step1._artifacts_generated is True

    # New instance sees cached artifacts and skips generation
    step2 = ref_jc_factory()
    step2.run()
    assert step2._artifacts_generated is False

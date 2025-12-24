"""Tests for CreateRefTnJcStep and CreateRefTnEndSeqsStep."""

import pytest

from amplifinder.steps import CreateRefTnJcStep, CreateRefTnEndSeqsStep
from amplifinder.data_types import RecordTypedDF, Side


@pytest.fixture
def ref_jc_step(locate_tns_step, tiny_genome, tmp_output):
    """Create reference junctions step."""
    tn_loc = locate_tns_step.run()
    return CreateRefTnJcStep(
        tn_loc=tn_loc,
        genome=tiny_genome,
        output_dir=tmp_output,
        source="isfinder",
        reference_tn_out_span=50,
    )


def test_creates_junction_records(ref_jc_step):
    """Should create 2 junctions per TN (left + right)."""
    ref_jc = ref_jc_step.run()

    assert isinstance(ref_jc, RecordTypedDF)
    assert len(ref_jc.df) == 4  # 2 TNs * 2 sides
    assert set(jc.ref_tn_side.side for jc in ref_jc) == {Side.LEFT, Side.RIGHT}


def test_junction_positions_correct(ref_jc_step, locate_tns_step):
    """Junction positions should match TN boundaries."""
    ref_jc = ref_jc_step.run()
    tn_loc = locate_tns_step.load_outputs()

    tn1 = tn_loc.df.iloc[0]
    left_jc = next(jc for jc in ref_jc if jc.ref_tn_side.tn_id == tn1["tn_id"] and jc.ref_tn_side.side == Side.LEFT)
    right_jc = next(jc for jc in ref_jc if jc.ref_tn_side.tn_id == tn1["tn_id"] and jc.ref_tn_side.side == Side.RIGHT)

    assert left_jc.pos1 == tn1["loc_left"]
    assert right_jc.pos1 == tn1["loc_right"]


def test_skips_if_exists(ref_jc_step):
    """Should skip if output exists."""
    ref_jc_step.run()
    assert ref_jc_step.run_count == 1
    ref_jc_step.run()
    assert ref_jc_step.run_count == 1  # didn't run again


# --- CreateRefTnEndSeqsStep ---

@pytest.fixture
def end_seqs_step(ref_jc_step, locate_tns_step, tiny_genome, tmp_output):
    """Create TN end sequences step."""
    ref_jc = ref_jc_step.run()
    tn_loc = locate_tns_step.run()
    return CreateRefTnEndSeqsStep(
        ref_tn_jc=ref_jc,
        tn_loc=tn_loc,
        genome=tiny_genome,
        output_dir=tmp_output,
        source="isfinder",
        max_dist_to_tn=20,
    )


def test_creates_end_sequences(end_seqs_step):
    """Should create end sequences for each junction."""
    end_seqs = end_seqs_step.run()

    assert isinstance(end_seqs, RecordTypedDF)
    assert len(end_seqs.df) == 4  # 4 junctions
    assert all(len(seq.seq_fwd) > 0 for seq in end_seqs)
    assert all(len(seq.seq_rc) > 0 for seq in end_seqs)

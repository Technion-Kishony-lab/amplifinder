"""Tests for CreateReferenceTnJunctionsStep and CreateRefTnEndSeqsStep."""

import pytest

from amplifinder.steps import CreateReferenceTnJunctionsStep, CreateRefTnEndSeqsStep
from amplifinder.data_types import RecordTypedDF, Side


@pytest.fixture
def ref_jc_step(locate_tns_step, tmp_output):
    """Create reference junctions step."""
    tn_loc = locate_tns_step.run_and_read_outputs()
    return CreateReferenceTnJunctionsStep(
        tn_loc=tn_loc,
        ref_name="tiny",
        output_dir=tmp_output,
        reference_tn_out_span=50,
    )


def test_creates_junction_records(ref_jc_step):
    """Should create 2 junctions per TN (left + right)."""
    ref_jc = ref_jc_step.run_and_read_outputs()

    assert isinstance(ref_jc, RecordTypedDF)
    assert len(ref_jc.df) == 4  # 2 TNs * 2 sides
    assert set(ref_jc.df["tn_side"]) == {Side.LEFT, Side.RIGHT}


def test_junction_positions_correct(ref_jc_step, locate_tns_step):
    """Junction positions should match TN boundaries."""
    ref_jc = ref_jc_step.run_and_read_outputs()
    tn_loc = locate_tns_step.read_outputs()

    tn1 = tn_loc.df.iloc[0]
    left_jc = ref_jc.df[(ref_jc.df["refTN"] == tn1["ID"]) & (ref_jc.df["tn_side"] == Side.LEFT)].iloc[0]
    right_jc = ref_jc.df[(ref_jc.df["refTN"] == tn1["ID"]) & (ref_jc.df["tn_side"] == Side.RIGHT)].iloc[0]

    assert left_jc["pos1"] == tn1["LocLeft"]
    assert right_jc["pos1"] == tn1["LocRight"]


def test_skips_if_exists(ref_jc_step):
    """Should skip if output exists."""
    assert ref_jc_step.run() is True
    assert ref_jc_step.run() is False


# --- CreateRefTnEndSeqsStep ---

@pytest.fixture
def end_seqs_step(ref_jc_step, tiny_genome, tmp_output):
    """Create TN end sequences step."""
    ref_jc = ref_jc_step.run_and_read_outputs()
    return CreateRefTnEndSeqsStep(
        ref_tn_jc=ref_jc,
        genome=tiny_genome,
        ref_path=tmp_output,
        max_dist_to_tn=20,
    )


def test_creates_end_sequences(end_seqs_step):
    """Should create end sequences for each junction."""
    end_seqs = end_seqs_step.run_and_read_outputs()

    assert isinstance(end_seqs, RecordTypedDF)
    assert len(end_seqs.df) == 4  # 4 junctions
    assert all(len(s) > 0 for s in end_seqs.df["seq_fwd"])
    assert all(len(s) > 0 for s in end_seqs.df["seq_rc"])

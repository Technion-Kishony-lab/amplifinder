"""Tests for CreateTNJC2Step."""

import pytest
from typing import List

from amplifinder.steps import (
    CreateReferenceTnJunctionsStep,
    CreateRefTnEndSeqsStep,
    CreateTNJCStep,
    CreateTNJC2Step,
)
from amplifinder.data_types import RecordTypedDF, TnJunction, TnJunctionPair, Junction, TnMatch, Side


# =============================================================================
# Helpers
# =============================================================================

def make_tnjc(
    num: int,
    pos2: int,
    dir2: int,
    tn_side: Side,
    scaf: str = "tiny",
    tn_id: int = 1,
) -> TnJunction:
    """Create a TnJunction with sensible defaults for testing."""
    return TnJunction(
        num=num,
        scaf1=scaf, pos1=100 * num, dir1=1 if tn_side == Side.LEFT else -1,
        scaf2=scaf, pos2=pos2, dir2=dir2,
        flanking_left=50, flanking_right=50,
        matches=[TnMatch(tn_id=tn_id, side=tn_side, distance=0)],
        switched=False,
    )


def run_tnjc2(tnjc_records: List[TnJunction], genome, output_dir) -> RecordTypedDF[TnJunctionPair]:
    """Create TNJC2 from junction records."""
    tnjc = RecordTypedDF.from_records(tnjc_records, TnJunction)
    return CreateTNJC2Step(
        tnjc=tnjc,
        genome=genome,
        output_dir=output_dir,
    ).run_and_read_outputs()


@pytest.fixture
def tnjc(locate_tns_step, tiny_genome, tmp_output):
    """Create TNJC (TN-associated junctions) from test data."""
    tn_loc = locate_tns_step.run_and_read_outputs()

    ref_jc = CreateReferenceTnJunctionsStep(
        tn_loc=tn_loc,
        ref_name="tiny",
        output_dir=tmp_output,
        reference_tn_out_span=50,
    ).run_and_read_outputs()

    ref_tn_end_seqs = CreateRefTnEndSeqsStep(
        ref_tn_jc=ref_jc,
        genome=tiny_genome,
        ref_path=tmp_output,
        max_dist_to_tn=20,
    ).run_and_read_outputs()

    # Convert to base Junction type
    junctions = [
        Junction(
            num=jc.num, scaf1=jc.scaf1, pos1=jc.pos1, dir1=jc.dir1,
            scaf2=jc.scaf2, pos2=jc.pos2, dir2=jc.dir2,
            flanking_left=jc.flanking_left, flanking_right=jc.flanking_right,
        )
        for jc in ref_jc
    ]
    mock_junctions = RecordTypedDF.from_records(junctions, Junction)

    return CreateTNJCStep(
        jc_df=mock_junctions,
        ref_tn_end_seqs=ref_tn_end_seqs,
        genome=tiny_genome,
        output_dir=tmp_output,
        max_dist_to_tn=20,
        trim_jc_flanking=5,
    ).run_and_read_outputs()


@pytest.fixture
def tnjc2_step(tnjc, tiny_genome, tmp_output):
    """Create TNJC2 step."""
    return CreateTNJC2Step(
        tnjc=tnjc,
        genome=tiny_genome,
        output_dir=tmp_output,
    )


def test_runs_without_error(tnjc2_step):
    """Should run and produce output file."""
    tnjc2 = tnjc2_step.run_and_read_outputs()

    assert isinstance(tnjc2, RecordTypedDF)
    assert tnjc2_step.output_file.exists()


def test_output_has_correct_columns(tnjc2_step):
    """TNJC2 output should have expected columns."""
    tnjc2 = tnjc2_step.run_and_read_outputs()

    expected_cols = {
        "jc_num_L", "jc_num_R", "scaf_chr",
        "pos_chr_L", "pos_chr_R", "pos_tn_L", "pos_tn_R",
        "dir_chr_L", "dir_chr_R", "dir_tn_L", "dir_tn_R",
        "tn_ids", "tn_orientation", "span_origin",
        "amplicon_length", "complementary_length",
    }
    assert expected_cols.issubset(set(tnjc2.df.columns))


def test_skips_if_exists(tnjc2_step):
    """Should skip if output exists."""
    assert tnjc2_step.run() is True
    assert tnjc2_step.run() is False


def test_pairs_opposing_junctions(tiny_genome, tmp_output):
    """Should pair junctions on same scaffold with opposite directions."""
    tnjc2 = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=1, tn_side=Side.LEFT),   # points right
        make_tnjc(num=2, pos2=1000, dir2=-1, tn_side=Side.RIGHT),  # points left
    ], tiny_genome, tmp_output)

    assert len(tnjc2) == 1
    pair = list(tnjc2)[0]
    assert pair.scaf_chr == "tiny"
    assert pair.pos_chr_L == 500
    assert pair.pos_chr_R == 1000
    assert 1 in pair.tn_ids


def test_no_pair_same_direction(tiny_genome, tmp_output):
    """Should not pair junctions with same chromosome direction."""
    tnjc2 = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=1, tn_side=Side.LEFT),   # points right
        make_tnjc(num=2, pos2=1000, dir2=1, tn_side=Side.RIGHT),  # also right - no pair!
    ], tiny_genome, tmp_output)

    assert len(tnjc2) == 0


def test_no_pair_different_scaffold(tiny_genome, tmp_output):
    """Should not pair junctions on different scaffolds."""
    tnjc2 = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=1, tn_side=Side.LEFT, scaf="tiny"),
        make_tnjc(num=2, pos2=1000, dir2=-1, tn_side=Side.RIGHT, scaf="other"),
    ], tiny_genome, tmp_output)

    assert len(tnjc2) == 0


def test_no_pair_same_tn_side(tiny_genome, tmp_output):
    """Should not pair junctions matching same TN side."""
    tnjc2 = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=1, tn_side=Side.LEFT),
        make_tnjc(num=2, pos2=1000, dir2=-1, tn_side=Side.LEFT),  # also left - no pair!
    ], tiny_genome, tmp_output)

    assert len(tnjc2) == 0


def test_calculates_amplicon_length(tiny_genome, tmp_output):
    """Should calculate amplicon length correctly."""
    tnjc2 = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=1, tn_side=Side.LEFT),
        make_tnjc(num=2, pos2=1000, dir2=-1, tn_side=Side.RIGHT),
    ], tiny_genome, tmp_output)

    pair = list(tnjc2)[0]
    # Length should be 1000 - 500 + 1 = 501
    assert pair.amplicon_length == 501


def test_handles_span_origin(tiny_genome, tmp_output):
    """Should detect span_origin when left junction points left."""
    tnjc2 = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=-1, tn_side=Side.LEFT),   # points left -> span_origin
        make_tnjc(num=2, pos2=1000, dir2=1, tn_side=Side.RIGHT),  # points right
    ], tiny_genome, tmp_output)

    pair = list(tnjc2)[0]
    assert pair.span_origin is True


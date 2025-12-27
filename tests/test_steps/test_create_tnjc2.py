"""Tests for PairTnJcToRawTnJc2Step."""

import pytest
from typing import List

from amplifinder.steps import (
    CreateRefTnJcStep,
    CreateRefTnEndSeqsStep,
    CreateTnJcStep,
    PairTnJcToRawTnJc2Step,
)
from amplifinder.data_types import RecordTypedDf, TnJunction, RawTnJc2, Junction, RefTnSide, Side, Orientation


# =============================================================================
# Helpers
# =============================================================================

def make_tnjc(
    num: int,
    pos2: int,
    dir2: Orientation,
    tn_side: Side,
    scaf: str = "tiny",
    tn_id: int = 1,
) -> TnJunction:
    """Create a TnJunction with sensible defaults for testing."""
    return TnJunction(
        num=num,
        scaf1=scaf, pos1=100 * num, dir1=Orientation.FORWARD if tn_side == Side.LEFT else Orientation.REVERSE,
        scaf2=scaf, pos2=pos2, dir2=dir2,
        flanking_left=50, flanking_right=50,
        ref_tn_sides=[RefTnSide(tn_id=tn_id, side=tn_side, distance=0)],
        switched=False,
    )


def run_tnjc2(tnjc_records: List[TnJunction], genome, output_dir) -> RecordTypedDf[RawTnJc2]:
    """Create RawTnJc2 from junction records."""
    tnjcs = RecordTypedDf.from_records(tnjc_records, TnJunction)
    return PairTnJcToRawTnJc2Step(
        tnjcs=tnjcs,
        genome=genome,
        output_dir=output_dir,
    ).run()


@pytest.fixture
def tnjc(locate_tns_step, tiny_genome, tmp_output):
    """Create TnJc (TN-associated junctions) from test data."""
    tn_loc = locate_tns_step.run()

    ref_jc = CreateRefTnJcStep(
        ref_tn_locs=tn_loc,
        genome=tiny_genome,
        output_dir=tmp_output,
        source="isfinder",
        reference_tn_out_span=50,
    ).run()

    ref_tn_end_seqs = CreateRefTnEndSeqsStep(
        ref_tn_jcs=ref_jc,
        ref_tn_locs=tn_loc,
        genome=tiny_genome,
        output_dir=tmp_output,
        source="isfinder",
        max_dist_to_tn=20,
    ).run()

    # Convert to base Junction type
    junctions = [
        Junction(
            num=jc.num, scaf1=jc.scaf1, pos1=jc.pos1, dir1=jc.dir1,
            scaf2=jc.scaf2, pos2=jc.pos2, dir2=jc.dir2,
            flanking_left=jc.flanking_left, flanking_right=jc.flanking_right,
        )
        for jc in ref_jc
    ]
    mock_junctions = RecordTypedDf.from_records(junctions, Junction)

    return CreateTnJcStep(
        junctions=mock_junctions,
        ref_tn_end_seqs=ref_tn_end_seqs,
        genome=tiny_genome,
        output_dir=tmp_output,
        max_dist_to_tn=20,
        trim_jc_flanking=5,
    ).run()


@pytest.fixture
def tnjc2_step(tnjc, tiny_genome, tmp_output):
    """Create TnJc2 step."""
    return PairTnJcToRawTnJc2Step(
        tnjcs=tnjc,
        genome=tiny_genome,
        output_dir=tmp_output,
    )


def test_runs_without_error(tnjc2_step):
    """Should run and produce output file."""
    raw_tnjc2s = tnjc2_step.run()

    assert isinstance(raw_tnjc2s, RecordTypedDf)
    assert tnjc2_step.output_file.exists()


def test_output_has_correct_columns(tnjc2_step):
    """RawTnJc2 output should have expected columns."""
    raw_tnjc2s = tnjc2_step.run()

    expected_cols = {
        "jc_num_L", "jc_num_R", "scaf",
        "pos_scaf_L", "pos_scaf_R", "pos_tn_L", "pos_tn_R",
        "dir_scaf_L", "dir_scaf_R", "dir_tn_L", "dir_tn_R",
        "tn_ids", "tn_orientations", "span_origin",
        "amplicon_length", "complementary_length",
    }
    assert expected_cols.issubset(set(raw_tnjc2s.df.columns))


def test_skips_if_exists(tnjc2_step):
    """Should skip if output exists."""
    tnjc2_step.run()
    assert tnjc2_step.run_count == 1
    tnjc2_step.run()
    assert tnjc2_step.run_count == 1  # didn't run again


def test_pairs_opposing_junctions(tiny_genome, tmp_output):
    """Should pair junctions on same scaffold with opposite directions."""
    raw_tnjc2s = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=Orientation.FORWARD, tn_side=Side.LEFT),   # points right
        make_tnjc(num=2, pos2=1000, dir2=Orientation.REVERSE, tn_side=Side.RIGHT),  # points left
    ], tiny_genome, tmp_output)

    assert len(raw_tnjc2s) == 1
    pair = list(raw_tnjc2s)[0]
    assert pair.scaf == "tiny"
    assert pair.pos_scaf_L == 500
    assert pair.pos_scaf_R == 1000
    assert 1 in pair.tn_ids


def test_no_pair_same_direction(tiny_genome, tmp_output):
    """Should not pair junctions with same chromosome direction."""
    raw_tnjc2s = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=Orientation.FORWARD, tn_side=Side.LEFT),   # points right
        make_tnjc(num=2, pos2=1000, dir2=Orientation.FORWARD, tn_side=Side.RIGHT),  # also right - no pair!
    ], tiny_genome, tmp_output)

    assert len(raw_tnjc2s) == 0


def test_no_pair_different_scaffold(tiny_genome, tmp_output):
    """Should not pair junctions on different scaffolds."""
    raw_tnjc2s = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=Orientation.FORWARD, tn_side=Side.LEFT, scaf="tiny"),
        make_tnjc(num=2, pos2=1000, dir2=Orientation.REVERSE, tn_side=Side.RIGHT, scaf="other"),
    ], tiny_genome, tmp_output)

    assert len(raw_tnjc2s) == 0


def test_no_pair_same_tn_side(tiny_genome, tmp_output):
    """Should not pair junctions matching same TN side."""
    raw_tnjc2s = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=Orientation.FORWARD, tn_side=Side.LEFT),
        make_tnjc(num=2, pos2=1000, dir2=Orientation.REVERSE, tn_side=Side.LEFT),  # also left - no pair!
    ], tiny_genome, tmp_output)

    assert len(raw_tnjc2s) == 0


def test_calculates_amplicon_length(tiny_genome, tmp_output):
    """Should calculate amplicon length correctly."""
    raw_tnjc2s = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=Orientation.FORWARD, tn_side=Side.LEFT),
        make_tnjc(num=2, pos2=1000, dir2=Orientation.REVERSE, tn_side=Side.RIGHT),
    ], tiny_genome, tmp_output)

    pair = list(raw_tnjc2s)[0]
    # Length should be 1000 - 500 + 1 = 501
    assert pair.amplicon_length == 501


def test_handles_span_origin(tiny_genome, tmp_output):
    """Should detect span_origin when left junction points left."""
    raw_tnjc2s = run_tnjc2([
        make_tnjc(num=1, pos2=500, dir2=Orientation.REVERSE, tn_side=Side.LEFT),   # points left -> span_origin
        make_tnjc(num=2, pos2=1000, dir2=Orientation.FORWARD, tn_side=Side.RIGHT),  # points right
    ], tiny_genome, tmp_output)

    pair = list(raw_tnjc2s)[0]
    assert pair.span_origin is True

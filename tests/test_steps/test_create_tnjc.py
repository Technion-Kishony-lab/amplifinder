"""Tests for CreateTNJCStep."""

import pytest

from amplifinder.steps import CreateRefTnJcsStep, CreateRefTnEndSeqsStep, CreateTNJCStep
from amplifinder.data_types import RecordTypedDF, Junction


@pytest.fixture
def ref_tn_end_seqs(locate_tns_step, tiny_genome, tmp_output):
    """Create TN end sequences."""
    tn_loc = locate_tns_step.run()

    ref_jc = CreateRefTnJcsStep(
        tn_loc=tn_loc,
        ref_name="tiny",
        output_dir=tmp_output,
        reference_tn_out_span=50,
    ).run()

    return CreateRefTnEndSeqsStep(
        ref_tn_jc=ref_jc,
        genome=tiny_genome,
        ref_path=tmp_output,
        max_dist_to_tn=20,
    ).run()


@pytest.fixture
def mock_junctions(locate_tns_step, tmp_output):
    """Create mock junctions (using ref TN junctions as input)."""
    tn_loc = locate_tns_step.run()

    ref_jc = CreateRefTnJcsStep(
        tn_loc=tn_loc,
        ref_name="tiny",
        output_dir=tmp_output,
        reference_tn_out_span=50,
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
    return RecordTypedDF.from_records(junctions, Junction)


@pytest.fixture
def tnjc_step(mock_junctions, ref_tn_end_seqs, tiny_genome, tmp_output):
    """Create TNJC step."""
    return CreateTNJCStep(
        jc_df=mock_junctions,
        ref_tn_end_seqs=ref_tn_end_seqs,
        genome=tiny_genome,
        output_dir=tmp_output,
        max_dist_to_tn=20,
        trim_jc_flanking=5,
    )


def test_runs_without_error(tnjc_step):
    """Should run and produce output file."""
    tnjc = tnjc_step.run()

    assert isinstance(tnjc, RecordTypedDF)
    assert tnjc_step.output_file.exists()


def test_output_has_correct_columns(tnjc_step):
    """TNJC output should have expected columns."""
    tnjc = tnjc_step.run()

    expected_cols = {"num", "scaf1", "pos1", "dir1", "scaf2", "pos2", "dir2",
                     "flanking_left", "flanking_right", "matches", "switched"}
    assert expected_cols.issubset(set(tnjc.df.columns))


def test_skips_if_exists(tnjc_step):
    """Should skip if output exists."""
    tnjc_step.run()
    assert tnjc_step.run_count == 1
    tnjc_step.run()
    assert tnjc_step.run_count == 1  # didn't run again

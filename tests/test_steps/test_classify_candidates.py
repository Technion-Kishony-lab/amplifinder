"""Tests for ClassifyTnJc2CandidatesStep."""

import pytest
from amplifinder.steps import ClassifyTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, AnalyzedTnJc2, RawEvent, Orientation, Coverage,
)


@pytest.fixture
def sample_analyzed(tmp_path):
    """Create a sample AnalyzedTnJc2."""
    return AnalyzedTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf="tiny",
        pos_scaf_L=200, pos_scaf_R=300,
        pos_tn_L=10, pos_tn_R=20,
        dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[1], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="tiny", iso_name="sample1",
        amplicon_coverage=2.0,
        scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        iso_amplicon_coverage=Coverage(mean=2.0, median=2.0, mode=2.0),
        iso_scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1], chosen_tn_id=1,
        analysis_dir="jc_200_300_001_L150",
        jc_cov_left=[0] * 7, jc_cov_right=[0] * 7, jc_cov_spanning=[0] * 7,
        isolate_architecture=RawEvent.FLANKED,
        ancestor_architecture=None,
        event="flanked",
        event_modifiers=[],
    )


def test_step_initialization(sample_analyzed, tmp_path):
    """Should initialize step correctly."""
    analyzed_tnjc2s = RecordTypedDf.from_records([sample_analyzed], AnalyzedTnJc2)

    step = ClassifyTnJc2CandidatesStep(
        analyzed_tnjc2s=analyzed_tnjc2s,
        output_dir=tmp_path,
        has_ancestor=False,
    )

    assert step.analyzed_tnjc2s == analyzed_tnjc2s
    assert step.has_ancestor is False
    assert step.output_file == tmp_path / "tnjc2_analyzed.csv"

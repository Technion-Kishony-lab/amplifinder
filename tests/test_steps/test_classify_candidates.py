"""Tests for ClassifyTnJc2CandidatesStep."""

import pytest
from pathlib import Path
from amplifinder.steps import ClassifyTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, AnalyzedTnJc2, RawEvent, Orientation, EventModifier,
)


@pytest.fixture
def sample_analyzed(tmp_path):
    """Create a sample AnalyzedTnJc2."""
    return AnalyzedTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf_chr="tiny",
        pos_chr_L=200, pos_chr_R=300,
        pos_tn_L=10, pos_tn_R=20,
        dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[1], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="tiny", iso_name="sample1",
        amplicon_coverage=2.0, genome_coverage=1.0,
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1], chosen_tn_id=1,
        analysis_dir="jc_200_300_001_L150",
        jc_cov_left=[0]*7, jc_cov_right=[0]*7, jc_cov_spanning=[0]*7,
        isolate_architecture=RawEvent.FLANKED,
        ancestor_architecture=None,
        event="flanked",
        event_modifiers=[],
    )


def test_step_initialization(sample_analyzed, tmp_path):
    """Should initialize step correctly."""
    analyzed = RecordTypedDf.from_records([sample_analyzed], AnalyzedTnJc2)
    
    step = ClassifyTnJc2CandidatesStep(
        analyzed=analyzed,
        output_dir=tmp_path,
        has_ancestor=False,
    )
    
    assert step.analyzed == analyzed
    assert step.has_ancestor is False
    assert step.output_file == tmp_path / "tn_jc2_analyzed.csv"

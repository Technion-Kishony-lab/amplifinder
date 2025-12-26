"""Tests for FilterTnJc2CandidatesStep."""

import pytest
from pathlib import Path
from amplifinder.steps import FilterTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, RawEvent, Orientation,
)


@pytest.fixture
def sample_classified_tnjc2(tmp_path):
    """Create sample ClassifiedTnJc2 records."""
    records = [
        ClassifiedTnJc2(
            jc_num_L=1, jc_num_R=2,
            scaf_chr="chr1",
            pos_chr_L=100, pos_chr_R=200,
            pos_tn_L=10, pos_tn_R=20,
            dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[1], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=100, complementary_length=900,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=2.0, genome_coverage=1.0,
            copy_number=2.0, amplicon_coverage_mode=2.0,
            raw_event=RawEvent.FLANKED,
            shared_tn_ids=[1], chosen_tn_id=1,
        ),
        ClassifiedTnJc2(
            jc_num_L=3, jc_num_R=4,
            scaf_chr="chr1",
            pos_chr_L=300, pos_chr_R=320,
            pos_tn_L=30, pos_tn_R=40,
            dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[2], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=20,  # Too short
            complementary_length=980,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=1.0, genome_coverage=1.0,
            copy_number=1.0, amplicon_coverage_mode=1.0,
            raw_event=RawEvent.UNFLANKED,
            shared_tn_ids=[2], chosen_tn_id=2,
        ),
    ]
    return RecordTypedDf.from_records(records, ClassifiedTnJc2)


def test_filters_by_length(sample_classified_tnjc2, tmp_path):
    """Should filter candidates by amplicon length."""
    step = FilterTnJc2CandidatesStep(
        classified_tnjc2s=sample_classified_tnjc2,
        output_dir=tmp_path,
        min_amplicon_length=30,
        max_amplicon_length=1_000_000,
    )
    
    result = step.run()
    
    # Should filter out the short one (length 20)
    assert len(result) == 1
    result_list = list(result)
    assert result_list[0].amplicon_length == 100
    assert result_list[0].analysis_dir.startswith("jc_100_200")


def test_assigns_analysis_dir(sample_classified_tnjc2, tmp_path):
    """Should assign analysis directory names."""
    step = FilterTnJc2CandidatesStep(
        classified_tnjc2s=sample_classified_tnjc2,
        output_dir=tmp_path,
        min_amplicon_length=30,
    )
    
    result = step.run()
    
    assert len(result) == 1
    result_list = list(result)
    assert result_list[0].analysis_dir.startswith("jc_")
    assert "_L150" in result_list[0].analysis_dir

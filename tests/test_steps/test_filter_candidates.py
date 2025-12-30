"""Tests for FilterTnJc2CandidatesStep."""

import pytest
from amplifinder.steps import FilterTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, RawEvent,
)


@pytest.fixture
def sample_classified_tnjc2(classified_tnjc2_record):
    """Create sample ClassifiedTnJc2 records."""
    first = ClassifiedTnJc2.from_other(
        classified_tnjc2_record,
        scaf="chr1",
        start=100,
        end=200,
        pos_tn_S=10,
        pos_tn_E=20,
        amplicon_length=100,
        ref_name="U00096",
        iso_name="sample1",
        iso_amplicon_coverage=2.0,
        copy_number=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1],
        chosen_tn_id=1,
    )

    second = ClassifiedTnJc2.from_other(
        classified_tnjc2_record,
        jc_num_S=3,
        jc_num_E=4,
        scaf="chr1",
        start=300,
        end=320,
        pos_tn_S=30,
        pos_tn_E=40,
        tn_ids=[2],
        amplicon_length=20,  # Too short
        iso_amplicon_coverage=1.0,
        copy_number=1.0,
        raw_event=RawEvent.UNFLANKED,
        shared_tn_ids=[2],
        chosen_tn_id=2,
    )

    return RecordTypedDf.from_records([first, second], ClassifiedTnJc2)


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

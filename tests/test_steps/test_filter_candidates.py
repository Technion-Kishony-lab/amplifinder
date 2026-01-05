"""Tests for FilterTnJc2CandidatesStep."""

from amplifinder.data_types.enums import BaseRawEvent
import pytest
from amplifinder.steps import FilterTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, RawEvent,
)


@pytest.fixture
def sample_classified_tnjc2(classified_tnjc2_record):
    """Create sample ClassifiedTnJc2 records."""
    # Use different base_raw_event values to differentiate
    first = ClassifiedTnJc2.from_other(
        classified_tnjc2_record,
        base_raw_event=BaseRawEvent.LOCUS_JOINING,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )

    second = ClassifiedTnJc2.from_other(
        classified_tnjc2_record,
        base_raw_event=BaseRawEvent.TRANSPOSITION,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=1.0,
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

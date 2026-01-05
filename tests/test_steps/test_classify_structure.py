"""Tests for ClassifyTnJc2StructureStep."""

import pytest
from amplifinder.steps import ClassifyTnJc2StructureStep
from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, RefTn, RawEvent,
)


@pytest.fixture
def sample_covered_tnjc2(covered_tnjc2_record):
    """Create sample CoveredTnJc2 records - one of each major type."""
    # FLANKED: both TN sides match
    flanked = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
        tn_id_start=1,
        tn_id_end=1,
    )

    # TRANSPOSITION: reference location junction
    transposition = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        iso_scaf_avg=10.0,
        iso_amplicon_avg=1.0,
        tn_id_start=1,
        tn_id_end=2,
    )

    # UNFLANKED: no TN sides match
    unflanked = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=1.0,
        tn_id_start=None,
        tn_id_end=None,
    )

    return RecordTypedDf.from_records([flanked, transposition, unflanked], CoveredTnJc2)


@pytest.fixture
def sample_tn_locs(ref_tn_record):
    """Create sample RefTn records."""
    first = RefTn.from_other(ref_tn_record, scaf="chr1", start=10, end=20, tn_id=1, tn_name="IS1")
    second = RefTn.from_other(ref_tn_record, scaf="chr1", start=30, end=40, tn_id=2, tn_name="IS2")
    return RecordTypedDf.from_records([first, second], RefTn)


def test_classify_structure(sample_covered_tnjc2, sample_tn_locs, tmp_path, tiny_genome):
    """Should classify junction pairs."""
    step = ClassifyTnJc2StructureStep(
        covered_tnjc2s=sample_covered_tnjc2,
        genome=tiny_genome,
        tn_locs=sample_tn_locs,
        output_dir=tmp_path,
    )

    result = step.run()

    assert len(result) == 3
    result_list = list(result)
    # Test each of the 3 types
    assert result_list[0].raw_event == RawEvent.FLANKED
    assert result_list[1].raw_event == RawEvent.TRANSPOSITION
    assert result_list[2].raw_event == RawEvent.UNFLANKED


def test_filters_by_length(sample_covered_tnjc2, sample_tn_locs, covered_tnjc2_record, tmp_path, tiny_genome):
    """Should filter candidates by amplicon length."""
    # Create a short amplicon
    short_record = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=1.0,
        tn_id_start=1,
        tn_id_end=1,
    )

    all_records = RecordTypedDf.from_records(
        list(sample_covered_tnjc2) + [short_record],
        CoveredTnJc2
    )

    step = ClassifyTnJc2StructureStep(
        covered_tnjc2s=all_records,
        genome=tiny_genome,
        tn_locs=sample_tn_locs,
        output_dir=tmp_path,
        transposition_threshold=30,
    )

    result = step.run()
    # Now expecting 4 records (3 from fixture + 1 new)
    assert len(result) == 4

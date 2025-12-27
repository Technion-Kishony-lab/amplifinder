"""Tests for ClassifyTnJc2StructureStep."""

import pytest
from amplifinder.steps import ClassifyTnJc2StructureStep
from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, RefTnLoc, RawEvent, Coverage,
)


@pytest.fixture
def sample_covered_tnjc2(covered_tnjc2_record):
    """Create sample CoveredTnJc2 records."""
    first = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        scaf="chr1",
        pos_scaf_L=100,
        pos_scaf_R=200,
        pos_tn_L=10,
        pos_tn_R=20,
        ref_name="U00096",
        iso_name="sample1",
        amplicon_coverage=2.0,
        iso_amplicon_coverage=Coverage(mean=2.0, median=2.0, mode=2.0),
        copy_number=2.0,
        amplicon_coverage_mode=2.0,
    )

    second = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        jc_num_L=0,
        jc_num_R=0,  # reference junctions
        scaf="chr1",
        pos_scaf_L=300,
        pos_scaf_R=400,
        pos_tn_L=30,
        pos_tn_R=40,
        tn_ids=[2],
        ref_name="U00096",
        iso_name="sample1",
        amplicon_coverage=1.0,
        iso_amplicon_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        copy_number=1.0,
        amplicon_coverage_mode=1.0,
    )

    return RecordTypedDf.from_records([first, second], CoveredTnJc2)


@pytest.fixture
def sample_tn_locs(ref_tn_loc_record):
    """Create sample RefTnLoc records."""
    first = RefTnLoc.from_other(ref_tn_loc_record, tn_scaf="chr1", loc_left=10, loc_right=20, tn_id=1, tn_name="IS1")
    second = RefTnLoc.from_other(ref_tn_loc_record, tn_scaf="chr1", loc_left=30, loc_right=40, tn_id=2, tn_name="IS2")
    return RecordTypedDf.from_records([first, second], RefTnLoc)


def test_classify_structure(sample_covered_tnjc2, sample_tn_locs, tmp_path):
    """Should classify junction pairs."""
    step = ClassifyTnJc2StructureStep(
        covered_tnjc2s=sample_covered_tnjc2,
        tn_locs=sample_tn_locs,
        output_dir=tmp_path,
    )

    result = step.run()

    assert len(result) == 2
    result_list = list(result)
    # First should be unflanked (no shared IS)
    assert result_list[0].raw_event == RawEvent.UNFLANKED
    # Second should be reference (both junctions are reference)
    assert result_list[1].raw_event == RawEvent.REFERENCE


def test_filters_by_length(sample_covered_tnjc2, sample_tn_locs, covered_tnjc2_record, tmp_path):
    """Should filter candidates by amplicon length."""
    # Create a short amplicon
    short_record = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        jc_num_L=3,
        jc_num_R=4,
        scaf="chr1",
        pos_scaf_L=500,
        pos_scaf_R=520,
        pos_tn_L=50,
        pos_tn_R=60,
        tn_ids=[1],
        amplicon_length=20,  # Too short
        complementary_length=980,
        ref_name="U00096",
        iso_name="sample1",
        amplicon_coverage=1.0,
        iso_amplicon_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        copy_number=1.0,
        amplicon_coverage_mode=1.0,
    )

    all_records = RecordTypedDf.from_records(
        list(sample_covered_tnjc2) + [short_record],
        CoveredTnJc2
    )

    step = ClassifyTnJc2StructureStep(
        covered_tnjc2s=all_records,
        tn_locs=sample_tn_locs,
        output_dir=tmp_path,
        min_amplicon_length=30,
    )

    result = step.run()
    # Short amplicon should be classified as transposition
    transpositions = [r for r in result if r.raw_event == RawEvent.TRANSPOSITION]
    assert len(transpositions) == 1

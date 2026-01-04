"""Tests for ClassifyTnJc2StructureStep."""

import pytest
from amplifinder.steps import ClassifyTnJc2StructureStep
from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, RefTn, RawEvent,
)


@pytest.fixture
def sample_covered_tnjc2(covered_tnjc2_record):
    """Create sample CoveredTnJc2 records."""
    first = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        scaf="chr1",
        start=100,
        end=200,
        pos_tn_S=10,
        pos_tn_E=20,
        ref_name="U00096",
        iso_name="sample1",
        iso_amplicon_coverage=2.0,
        copy_number=2.0,
    )

    second = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        jc_num_S=0,
        jc_num_E=0,  # reference junctions
        scaf="chr1",
        start=300,
        end=400,
        pos_tn_S=30,
        pos_tn_E=40,
        tn_ids=[2],
        ref_name="U00096",
        iso_name="sample1",
        iso_amplicon_coverage=1.0,
        copy_number=1.0,
    )

    return RecordTypedDf.from_records([first, second], CoveredTnJc2)


@pytest.fixture
def sample_tn_locs(ref_tn_loc_record):
    """Create sample RefTn records."""
    first = RefTn.from_other(ref_tn_loc_record, scaf="chr1", start=10, end=20, tn_id=1, tn_name="IS1")
    second = RefTn.from_other(ref_tn_loc_record, scaf="chr1", start=30, end=40, tn_id=2, tn_name="IS2")
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

    assert len(result) == 2
    result_list = list(result)
    # First should be unflanked (no shared IS)
    assert result_list[0].raw_event == RawEvent.UNFLANKED
    # Second should be reference (both junctions are reference)
    assert result_list[1].raw_event == RawEvent.REFERENCE


def test_filters_by_length(sample_covered_tnjc2, sample_tn_locs, covered_tnjc2_record, tmp_path, tiny_genome):
    """Should filter candidates by amplicon length."""
    # Create a short amplicon
    short_record = CoveredTnJc2.from_other(
        covered_tnjc2_record,
        jc_num_S=3,
        jc_num_E=4,
        scaf="chr1",
        start=500,
        end=520,
        pos_tn_S=50,
        pos_tn_E=60,
        tn_ids=[1],
        amplicon_length=20,  # Too short
        ref_name="U00096",
        iso_name="sample1",
        iso_amplicon_coverage=1.0,
        copy_number=1.0,
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
    # Short amplicon should be classified as transposition
    transpositions = [r for r in result if r.raw_event == RawEvent.TRANSPOSITION]
    assert len(transpositions) == 1

"""Tests for ClassifyTnJc2StructureStep."""

import pytest
from pathlib import Path
from amplifinder.steps import ClassifyTnJc2StructureStep
from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, RefTnLoc, RawEvent, Orientation, Coverage,
)


@pytest.fixture
def sample_covered_tnjc2(tmp_path):
    """Create sample CoveredTnJc2 records."""
    records = [
        CoveredTnJc2(
            jc_num_L=1, jc_num_R=2,
            scaf="chr1",
            pos_scaf_L=100, pos_scaf_R=200,
            pos_tn_L=10, pos_tn_R=20,
            dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[1], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=100, complementary_length=900,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=2.0, scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
            copy_number=2.0, amplicon_coverage_mode=2.0,
        ),
        CoveredTnJc2(
            jc_num_L=0, jc_num_R=0,  # reference junctions
            scaf="chr1",
            pos_scaf_L=300, pos_scaf_R=400,
            pos_tn_L=30, pos_tn_R=40,
            dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[2], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=100, complementary_length=900,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=1.0, scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
            copy_number=1.0, amplicon_coverage_mode=1.0,
        ),
    ]
    return RecordTypedDf.from_records(records, CoveredTnJc2)


@pytest.fixture
def sample_tn_locs():
    """Create sample RefTnLoc records."""
    records = [
        RefTnLoc(tn_id=1, tn_name="IS1", tn_scaf="chr1", loc_left=10, loc_right=20, complement=False, join=False),
        RefTnLoc(tn_id=2, tn_name="IS2", tn_scaf="chr1", loc_left=30, loc_right=40, complement=False, join=False),
    ]
    return RecordTypedDf.from_records(records, RefTnLoc)


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


def test_filters_by_length(sample_covered_tnjc2, sample_tn_locs, tmp_path):
    """Should filter candidates by amplicon length."""
    # Create a short amplicon
    short_record = CoveredTnJc2(
        jc_num_L=3, jc_num_R=4,
        scaf="chr1",
        pos_scaf_L=500, pos_scaf_R=520,
        pos_tn_L=50, pos_tn_R=60,
        dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[1], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=20,  # Too short
        complementary_length=980,
        ref_name="U00096", iso_name="sample1",
        amplicon_coverage=1.0, genome_coverage=1.0,
        copy_number=1.0, amplicon_coverage_mode=1.0,
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

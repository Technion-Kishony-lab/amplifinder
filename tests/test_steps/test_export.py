"""Tests for ExportTnJc2Step."""

import pytest
from amplifinder.steps import ExportTnJc2Step
from amplifinder.data_types import (
    RecordTypedDf, AnalyzedTnJc2, RawEvent, Orientation, Coverage,
)


@pytest.fixture
def sample_analyzed(tmp_path):
    """Create sample AnalyzedTnJc2 records."""
    records = [
        AnalyzedTnJc2(
            jc_num_L=1, jc_num_R=2,
            scaf="chr1",
            pos_scaf_L=100, pos_scaf_R=200,
            pos_tn_L=10, pos_tn_R=20,
            dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[1], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=150, complementary_length=850,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=2.5,
            scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
            iso_amplicon_coverage=Coverage(mean=2.5, median=2.5, mode=2.5),
            iso_scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
            copy_number=2.5, amplicon_coverage_mode=2.5,
            raw_event=RawEvent.FLANKED,
            shared_tn_ids=[1], chosen_tn_id=1,
            analysis_dir="jc_100_200_001_L150",
            jc_cov_left=[0] * 7, jc_cov_right=[0] * 7, jc_cov_spanning=[0] * 7,
            isolate_architecture=RawEvent.FLANKED,
            ancestor_architecture=None,
            event="flanked",
            event_modifiers=[],
        ),
        AnalyzedTnJc2(
            jc_num_L=3, jc_num_R=4,
            scaf="chr1",
            pos_scaf_L=300, pos_scaf_R=400,
            pos_tn_L=30, pos_tn_R=40,
            dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[2], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=50,  # Too short for filtering
            complementary_length=950,
            ref_name="U00096", iso_name="sample1",
        amplicon_coverage=0.2,
        scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        iso_amplicon_coverage=Coverage(mean=0.2, median=0.2, mode=0.2),
        iso_scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
            copy_number=0.2, amplicon_coverage_mode=0.2,
            raw_event=RawEvent.UNFLANKED,
            shared_tn_ids=[2], chosen_tn_id=2,
            analysis_dir="jc_300_400_002_L150",
            jc_cov_left=[0] * 7, jc_cov_right=[0] * 7, jc_cov_spanning=[0] * 7,
            isolate_architecture=RawEvent.UNFLANKED,
            ancestor_architecture=None,
            event="unflanked",
            event_modifiers=[],
        ),
    ]
    return RecordTypedDf.from_records(records, AnalyzedTnJc2)


def test_export_creates_files(sample_analyzed, tmp_path):
    """Should create tnjc2_exported.csv and candidate_amplifications.csv."""
    step = ExportTnJc2Step(
        analyzed_tnjc2s=sample_analyzed,
        output_dir=tmp_path,
        ref_name="U00096",
        iso_name="sample1",
        copy_number_threshold=1.5,
        del_copy_number_threshold=0.3,
        filter_amplicon_length=100,
    )

    step.run()

    # Check that files were created
    exported_file = tmp_path / "tnjc2_exported.csv"
    candidates_file = tmp_path / "candidate_amplifications.csv"

    assert exported_file.exists()
    assert candidates_file.exists()

    # Check that tnjc2_exported.csv contains all candidates
    import pandas as pd
    df_all = pd.read_csv(exported_file)
    assert len(df_all) == 2

    # Check that candidate_amplifications.csv contains filtered candidates
    df_filtered = pd.read_csv(candidates_file)
    # First candidate should be included (copy_number > threshold, length > filter)
    # Second should be excluded (length too short)
    assert len(df_filtered) == 1
    assert df_filtered.iloc[0]['amplicon_length'] == 150

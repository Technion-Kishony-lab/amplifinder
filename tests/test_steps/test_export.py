"""Tests for ExportStep."""

import pytest
from pathlib import Path
from amplifinder.steps import ExportStep
from amplifinder.data_types import (
    RecordTypedDF, AnalyzedTnJc2, RawEvent, Orientation, EventModifier,
)


@pytest.fixture
def sample_analyzed(tmp_path):
    """Create sample AnalyzedTnJc2 records."""
    records = [
        AnalyzedTnJc2(
            jc_num_L=1, jc_num_R=2,
            scaf_chr="chr1",
            pos_chr_L=100, pos_chr_R=200,
            pos_tn_L=10, pos_tn_R=20,
            dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[1], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=150, complementary_length=850,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=2.5, genome_coverage=1.0,
            copy_number=2.5, amplicon_coverage_mode=2.5,
            raw_event=RawEvent.FLANKED,
            shared_tn_ids=[1], chosen_tn_id=1,
            analysis_dir="jc_100_200_001_L150",
            jc_cov_left=[0]*7, jc_cov_right=[0]*7, jc_cov_spanning=[0]*7,
            isolate_architecture=RawEvent.FLANKED,
            ancestor_architecture=None,
            event="flanked",
            event_modifiers=[],
        ),
        AnalyzedTnJc2(
            jc_num_L=3, jc_num_R=4,
            scaf_chr="chr1",
            pos_chr_L=300, pos_chr_R=400,
            pos_tn_L=30, pos_tn_R=40,
            dir_chr_L=Orientation.FORWARD, dir_chr_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[2], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=50,  # Too short for filtering
            complementary_length=950,
            ref_name="U00096", iso_name="sample1",
            amplicon_coverage=0.2, genome_coverage=1.0,
            copy_number=0.2, amplicon_coverage_mode=0.2,
            raw_event=RawEvent.UNFLANKED,
            shared_tn_ids=[2], chosen_tn_id=2,
            analysis_dir="jc_300_400_002_L150",
            jc_cov_left=[0]*7, jc_cov_right=[0]*7, jc_cov_spanning=[0]*7,
            isolate_architecture=RawEvent.UNFLANKED,
            ancestor_architecture=None,
            event="unflanked",
            event_modifiers=[],
        ),
    ]
    return RecordTypedDF.from_records(records, AnalyzedTnJc2)


def test_export_creates_files(sample_analyzed, tmp_path):
    """Should create ISJC2.csv and candidate_amplifications.csv."""
    step = ExportStep(
        analyzed_candidates=sample_analyzed,
        output_dir=tmp_path,
        copy_number_threshold=1.5,
        del_copy_number_threshold=0.3,
        filter_amplicon_length=100,
    )
    
    step.run()
    
    # Check that files were created
    isjc2_file = tmp_path / "ISJC2.csv"
    candidates_file = tmp_path / "candidate_amplifications.csv"
    
    assert isjc2_file.exists()
    assert candidates_file.exists()
    
    # Check that ISJC2.csv contains all candidates
    import pandas as pd
    df_all = pd.read_csv(isjc2_file)
    assert len(df_all) == 2
    
    # Check that candidate_amplifications.csv contains filtered candidates
    df_filtered = pd.read_csv(candidates_file)
    # First candidate should be included (copy_number > threshold, length > filter)
    # Second should be excluded (length too short)
    assert len(df_filtered) == 1
    assert df_filtered.iloc[0]['amplicon_length'] == 150

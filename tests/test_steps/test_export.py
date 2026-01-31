"""Tests for ExportTnJc2Step."""

import pytest
from amplifinder.steps import ExportTnJc2Step
from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, Architecture, EventDescriptor,
)


@pytest.fixture
def sample_analyzed(analyzed_tnjc2_record):
    """Create sample ClassifiedTnJc2 records."""
    first = ClassifiedTnJc2.from_other(
        analyzed_tnjc2_record,
        iso_architecture=Architecture.FLANKED,
        event_descriptors=[EventDescriptor.DENOVO_LEFT],
    )

    second = ClassifiedTnJc2.from_other(
        analyzed_tnjc2_record,
        iso_architecture=Architecture.UNFLANKED,
        event_descriptors=[EventDescriptor.DENOVO_LEFT],
    )

    return RecordTypedDf.from_records([first, second], ClassifiedTnJc2)


def test_export_creates_files(sample_analyzed, tmp_path, tiny_genome):
    """Should create tnjc2_exported.csv and candidate_amplifications.csv."""
    step = ExportTnJc2Step(
        classified_tnjc2s=sample_analyzed,
        genome=tiny_genome,
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
    # Both candidates should be included (both have length > filter_amplicon_length=100)
    # First: length=101, copy_number=2.5 > 1.5 threshold -> included
    # Second: length=101, copy_number=0.2 < 0.3 del_threshold -> included (deletion)
    assert len(df_filtered) == 2
    assert df_filtered.iloc[0]['amplicon_length'] == 101

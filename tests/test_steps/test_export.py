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


def test_export_creates_files(sample_analyzed, tmp_path, tiny_genome, read_lengths_fixture):
    """Should create amplifications.yaml file."""
    step = ExportTnJc2Step(
        classified_tnjc2s=sample_analyzed,
        output_dir=tmp_path,
        ref_name="U00096",
        iso_name="sample1",
        read_lengths=read_lengths_fixture,
    )

    result = step.run()

    # Check that YAML file was created
    yaml_file = tmp_path / "amplifications.yaml"
    assert yaml_file.exists()
    
    # Check export result
    assert 'sample' in result
    assert 'amplicons' in result
    assert len(result['amplicons']) == 2

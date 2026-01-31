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


def test_export_creates_files(sample_analyzed, filtered_tnjc2_record, tmp_path, tiny_genome, ref_tns_indexed):
    """Should create summary.yaml and summary.json files."""
    from amplifinder.steps.read_length import ReadLengths

    # Create a linked_tnjc2s RecordTypedDf from filtered_tnjc2_record
    linked_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], type(filtered_tnjc2_record))

    step = ExportTnJc2Step(
        classified_tnjc2s=sample_analyzed,
        linked_tnjc2s=linked_tnjc2s,
        output_dir=tmp_path,
        ref_name="U00096",
        iso_name="sample1",
        read_lengths=ReadLengths(read_len_iso=150, read_len_anc=None, jc_arm_len_iso=300, jc_arm_len_anc=None),
        ref_tns=ref_tns_indexed,
    )

    result = step.run()

    # Check that YAML and JSON files were created
    yaml_file = tmp_path / "summary.yaml"
    json_file = tmp_path / "summary.json"
    assert yaml_file.exists()
    assert json_file.exists()

    # Check export result
    assert 'sample' in result
    assert 'amplicons' in result
    assert len(result['amplicons']) == 2

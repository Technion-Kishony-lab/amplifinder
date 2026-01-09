"""Tests for ClassifyTnJc2CandidatesStep."""

import pytest
from amplifinder.steps import ClassifyTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, AnalyzedTnJc2,
)


@pytest.fixture
def sample_analyzed(analyzed_tnjc2_record):
    """Create a sample AnalyzedTnJc2."""
    return AnalyzedTnJc2.from_other(analyzed_tnjc2_record)


def test_step_initialization(sample_analyzed, tmp_path):
    """Should initialize step correctly."""
    analyzed_tnjc2s = RecordTypedDf.from_records([sample_analyzed], AnalyzedTnJc2)

    step = ClassifyTnJc2CandidatesStep(
        analyzed_tnjc2s=analyzed_tnjc2s,
        output_dir=tmp_path,
        has_ancestor=False,
    )

    assert step.analyzed_tnjc2s == analyzed_tnjc2s
    assert step.has_ancestor is False
    assert step.output_file == tmp_path / "tnjc2_classified_final.csv"

"""Tests for breseq output parsing using real breseq output."""

import pytest

from amplifinder.tools.breseq import parse_breseq_output, _convert_value, RECORD_TYPES
from tests.conftest import FIXTURES_DIR


@pytest.fixture
def breseq_output():
    """Parse real breseq output."""
    return parse_breseq_output(FIXTURES_DIR / "breseq")


def test_parse_output_gd(breseq_output):
    """Should parse output.gd file into DataFrames."""
    assert len(breseq_output["SNP"]) > 0
    assert len(breseq_output["JC"]) > 0
    assert len(breseq_output["DEL"]) > 0
    assert len(breseq_output["MOB"]) > 0


def test_parsed_columns_match_schema(breseq_output):
    """Parsed DataFrames should have columns matching required schema fields."""
    for record_type, schema in RECORD_TYPES.items():
        df = breseq_output[record_type]
        if df.empty:
            continue
        # Schema is tuple of Column(name, dtype, optional)
        required_names = [c.name for c in schema if not c.optional]
        for field in required_names:
            assert field in df.columns, f"{record_type} missing required field: {field}"


def test_convert_int():
    assert _convert_value("123", int) == 123
    assert _convert_value("123.0", int) == 123


def test_convert_float():
    assert _convert_value("1.5", float) == 1.5
    assert _convert_value("1e-10", float) == 1e-10


def test_convert_string():
    assert _convert_value("hello", str) == "hello"


def test_convert_na():
    assert _convert_value("NA", int) is None
    assert _convert_value("", str) is None


def test_missing_output_gd(tmp_path):
    """Should raise FileNotFoundError for missing output."""
    with pytest.raises(FileNotFoundError):
        parse_breseq_output(tmp_path / "nonexistent")

import pytest
import pandas as pd
import numpy as np
from typing import Optional, List, Dict, Tuple

from amplifinder.data_types.records import Column, Schema
from amplifinder.data_types.validate_and_cast_df import validate_and_cast_df, parse_compound


@pytest.fixture
def simple_schema():
    return Schema([
        Column("name", str),
        Column("value", int),
    ])


@pytest.fixture
def optional_schema():
    return Schema([
        Column("name", str),
        Column("notes", Optional[str], optional=True),
    ])


@pytest.mark.parametrize("value,expected_type,expected", [
    ("[1, 2, 3]", List[int], [1, 2, 3]),
    ("{'a': 1}", Dict[str, int], {'a': 1}),
    ("(1, 'x')", Tuple[int, str], (1, 'x')),
])
def test_parse_compund(value, expected_type, expected):
    assert parse_compound(value, expected_type) == expected


def test_parse_compound_na_passthrough_optional():
    assert pd.isna(parse_compound(np.nan, Optional[List[int]]))


def test_parse_compound_na_raises_for_non_optional():
    with pytest.raises(TypeError, match="not optional"):
        parse_compound(np.nan, List[int])


@pytest.mark.parametrize("value,expected_type,match", [
    ({'a': 1}, List[int], r"Input should be a valid list"),
    (['not_int'], List[int], r"Input should be a valid integer"),
])
def test_parse_compound_type_mismatch(value, expected_type, match):
    from pydantic import ValidationError
    with pytest.raises(ValidationError, match=match):
        parse_compound(value, expected_type)


def test_valid_df(simple_schema):
    df = pd.DataFrame({"name": ["a", "b"], "value": [1, 2]})
    result = validate_and_cast_df(df, simple_schema)
    assert list(result.columns) == ["name", "value"]


def test_missing_required(simple_schema):
    df = pd.DataFrame({"name": ["a"]})
    with pytest.raises(ValueError, match="Missing required"):
        validate_and_cast_df(df, simple_schema)


def test_missing_optional_ok(optional_schema):
    df = pd.DataFrame({"name": ["a"]})
    result = validate_and_cast_df(df, optional_schema)
    assert "name" in result.columns


def test_extra_columns_error(simple_schema):
    df = pd.DataFrame({"name": ["a"], "value": [1], "extra": [0]})
    with pytest.raises(ValueError, match="Extra columns"):
        validate_and_cast_df(df, simple_schema, check_extra=True)


def test_extra_columns_ok(simple_schema):
    df = pd.DataFrame({"name": ["a"], "value": [1], "extra": [0]})
    result = validate_and_cast_df(df, simple_schema, check_extra=False)
    assert "extra" in result.columns


def test_cast_types(simple_schema):
    df = pd.DataFrame({"name": ["a"], "value": ["1"]})  # value as string
    result = validate_and_cast_df(df, simple_schema, cast=True)
    assert result["value"].dtype == np.int64


def test_no_cast_type_error(simple_schema):
    df = pd.DataFrame({"name": ["a"], "value": ["1"]})  # wrong dtype
    with pytest.raises(TypeError, match="expected"):
        validate_and_cast_df(df, simple_schema, cast=False)


def test_compound_type_parsing():
    schema = Schema([Column("items", List[int])])
    df = pd.DataFrame({"items": ["[1, 2, 3]"]})
    result = validate_and_cast_df(df, schema)
    assert result["items"].iloc[0] == [1, 2, 3]

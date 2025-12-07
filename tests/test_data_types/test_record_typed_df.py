"""Test RecordTypedDF save/load with various field types."""

import pytest
from enum import Enum
from typing import List, NamedTuple, Tuple

from amplifinder.data_types import RecordTypedDF
from amplifinder.data_types.records import Record


class Color(str, Enum):
    RED = "red"
    BLUE = "blue"


class Direction(int, Enum):
    UP = 1
    DOWN = -1


class Match(NamedTuple):
    """NamedTuple with enum field (like TnMatch)."""
    id: int
    color: Color
    value: int


class SampleRecord(Record):
    """Sample record with various field types."""
    name: str
    count: int
    color: Color
    direction: Direction
    tags: List[int]
    pair: Tuple[int, int]


class RecordWithMatches(Record):
    """Record with List of NamedTuples containing enums."""
    name: str
    matches: List[Match]


@pytest.fixture
def sample_records():
    return [
        SampleRecord(name="a", count=1, color=Color.RED, direction=Direction.UP, tags=[1, 2], pair=(10, 20)),
        SampleRecord(name="b", count=2, color=Color.BLUE, direction=Direction.DOWN, tags=[3], pair=(30, 40)),
    ]


def test_save_load_csv(sample_records, tmp_path):
    """Test round-trip save/load to CSV."""
    csv_path = tmp_path / "test.csv"

    # Save
    df = RecordTypedDF.from_records(sample_records, SampleRecord)
    df.to_csv(csv_path)

    # Load
    df2 = RecordTypedDF.from_csv(csv_path, SampleRecord)

    # Check
    assert len(df2) == 2
    records = list(df2)

    r0 = records[0]
    assert r0.name == "a"
    assert r0.count == 1
    assert r0.color == Color.RED
    assert r0.direction == Direction.UP
    assert r0.tags == [1, 2]
    assert r0.pair == (10, 20)

    r1 = records[1]
    assert r1.name == "b"
    assert r1.count == 2
    assert r1.color == Color.BLUE
    assert r1.direction == Direction.DOWN
    assert r1.tags == [3]
    assert r1.pair == (30, 40)


def test_save_load_namedtuple_with_enum(tmp_path):
    """Test round-trip for List[NamedTuple] where NamedTuple contains Enum."""
    csv_path = tmp_path / "matches.csv"

    records = [
        RecordWithMatches(name="x", matches=[Match(1, Color.RED, 10), Match(2, Color.BLUE, 20)]),
        RecordWithMatches(name="y", matches=[Match(3, Color.RED, 30)]),
    ]

    # Save
    df = RecordTypedDF.from_records(records, RecordWithMatches)
    df.to_csv(csv_path)

    # Load
    df2 = RecordTypedDF.from_csv(csv_path, RecordWithMatches)

    # Check
    assert len(df2) == 2
    loaded = list(df2)

    r0 = loaded[0]
    assert r0.name == "x"
    assert len(r0.matches) == 2
    assert r0.matches[0] == (1, "red", 10)  # Tuple, enum as value
    assert r0.matches[1] == (2, "blue", 20)


def test_list_int_type_validation(tmp_path):
    """Test that List[int] validates element types on load."""
    from amplifinder.data_types.validate_and_cast_df import parse_compound
    from pydantic import ValidationError

    # This should raise ValidationError - list contains string, not int
    with pytest.raises(ValidationError, match=r"Input should be a valid integer"):
        parse_compound(["not_an_int"], List[int])

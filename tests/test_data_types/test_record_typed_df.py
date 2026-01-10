"""Test RecordTypedDF save/load with various field types."""

import pytest
from enum import Enum
from typing import List, NamedTuple, Tuple, ClassVar

from amplifinder.data_types import RecordTypedDf
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
    df = RecordTypedDf.from_records(sample_records, SampleRecord)
    df.to_csv(csv_path)

    # Load
    df2 = RecordTypedDf.from_csv(csv_path, SampleRecord)

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
    df = RecordTypedDf.from_records(records, RecordWithMatches)
    df.to_csv(csv_path)

    # Load
    df2 = RecordTypedDf.from_csv(csv_path, RecordWithMatches)

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


def test_save_load_with_index(sample_records, tmp_path):
    """Test round-trip save/load with custom index (from dict)."""
    csv_path = tmp_path / "test_indexed.csv"

    # Create from dict with custom keys
    records_dict = {
        "key1": sample_records[0],
        "key2": sample_records[1],
    }
    df = RecordTypedDf.from_dict(records_dict, SampleRecord)

    # Verify index before save
    assert list(df.keys()) == ["key1", "key2"]
    assert df["key1"].name == "a"
    assert df["key2"].name == "b"

    # Save (should auto-save index since it's meaningful)
    df.to_csv(csv_path)

    # Load with index
    df2 = RecordTypedDf.from_csv(csv_path, SampleRecord, index_col=0)

    # Verify index is preserved
    assert list(df2.keys()) == ["key1", "key2"]
    assert len(df2) == 2

    # Verify records match
    assert df2["key1"].name == "a"
    assert df2["key1"].count == 1
    assert df2["key1"].color == Color.RED

    assert df2["key2"].name == "b"
    assert df2["key2"].count == 2
    assert df2["key2"].color == Color.BLUE

    # Verify items() works
    items = list(df2.items())
    assert len(items) == 2
    assert items[0][0] == "key1"
    assert items[0][1].name == "a"
    assert items[1][0] == "key2"
    assert items[1][1].name == "b"


def test_to_records_round_trip(sample_records):
    """Test round-trip: records → from_records() → to_records() → same records."""
    # Convert to DataFrame
    df = RecordTypedDf.from_records(sample_records, SampleRecord)

    # Convert back to records
    records2 = df.to_records()

    # Verify count
    assert len(records2) == len(sample_records)

    # Verify each record matches field-by-field
    for orig, loaded in zip(sample_records, records2):
        assert orig.name == loaded.name
        assert orig.count == loaded.count
        assert orig.color == loaded.color
        assert orig.direction == loaded.direction
        assert orig.tags == loaded.tags
        assert orig.pair == loaded.pair


def test_csv_export_fields(tmp_path):
    """Test CSV_EXPORT_FIELDS functionality."""

    class InnerRecord(Record):
        """Inner record for testing nested objects."""
        x: int
        y: str

    class RecordWithComplexFields(Record):
        """Record with complex objects and CSV_EXPORT_FIELDS."""
        name: str
        inner: InnerRecord
        value: int

        CSV_EXPORT_FIELDS: ClassVar[List[str]] = ['name', 'value']  # Exclude 'inner'

    # Create records with complex objects
    records = [
        RecordWithComplexFields(name="a", inner=InnerRecord(x=1, y="foo"), value=10),
        RecordWithComplexFields(name="b", inner=InnerRecord(x=2, y="bar"), value=20),
    ]

    # Convert to DataFrame
    df = RecordTypedDf.from_records(records, RecordWithComplexFields)

    # Verify DataFrame contains all fields (including complex objects)
    assert 'name' in df.df.columns
    assert 'inner' in df.df.columns
    assert 'value' in df.df.columns

    # Verify inner is stored as object (not serialized)
    assert df.df['inner'].dtype == object
    assert isinstance(df.df.iloc[0]['inner'], InnerRecord)

    # Save to CSV
    csv_path = tmp_path / "test_csv_export.csv"
    df.to_csv(csv_path)

    # Verify CSV contains only CSV_EXPORT_FIELDS columns
    import pandas as pd
    csv_df = pd.read_csv(csv_path)
    assert set(csv_df.columns) == {'name', 'value'}  # No 'inner' column
    assert len(csv_df) == 2
    assert csv_df.iloc[0]['name'] == 'a'
    assert csv_df.iloc[0]['value'] == 10

    # Verify step can detect artifact files exist
    from amplifinder.steps.base import RecordTypedDfStep

    # Create a minimal step instance to test has_artifact_files
    class TestStep(RecordTypedDfStep[RecordWithComplexFields]):
        def _calculate_output(self):
            return RecordTypedDf.from_records([], RecordWithComplexFields)

    step = TestStep(output_file=csv_path)
    # has_artifact_files returns True when output file exists
    assert step.has_artifact_files()


def test_csv_export_fields_with_property(tmp_path):
    """Test CSV_EXPORT_FIELDS with properties."""
    from typing import ClassVar, List

    class RecordWithProperty(Record):
        """Record with property included in CSV_EXPORT_FIELDS."""
        base_value: int

        CSV_EXPORT_FIELDS: ClassVar[List[str]] = ['base_value', 'computed_value']

        @property
        def computed_value(self) -> int:
            """Computed property."""
            return self.base_value * 2

    records = [
        RecordWithProperty(base_value=5),
        RecordWithProperty(base_value=10),
    ]

    # Convert to DataFrame
    df = RecordTypedDf.from_records(records, RecordWithProperty)

    # Save to CSV
    csv_path = tmp_path / "test_property_export.csv"
    df.to_csv(csv_path)

    # Verify CSV contains both field and property
    import pandas as pd
    csv_df = pd.read_csv(csv_path)
    assert set(csv_df.columns) == {'base_value', 'computed_value'}
    assert csv_df.iloc[0]['base_value'] == 5
    assert csv_df.iloc[0]['computed_value'] == 10  # 5 * 2
    assert csv_df.iloc[1]['base_value'] == 10
    assert csv_df.iloc[1]['computed_value'] == 20  # 10 * 2


def test_backward_compatibility(sample_records, tmp_path):
    """Test that existing Record classes without CSV_EXPORT_FIELDS still work."""
    # SampleRecord doesn't have CSV_EXPORT_FIELDS set (defaults to None)
    # Verify it works as before

    # Convert to DataFrame
    df = RecordTypedDf.from_records(sample_records, SampleRecord)

    # Verify round-trip works
    records2 = df.to_records()
    assert len(records2) == len(sample_records)
    assert records2[0].name == sample_records[0].name

    # Verify CSV save/load still works
    csv_path = tmp_path / "test_backward.csv"
    df.to_csv(csv_path)
    df2 = RecordTypedDf.from_csv(csv_path, SampleRecord)
    assert len(df2) == 2
    assert df2.to_records()[0].name == "a"

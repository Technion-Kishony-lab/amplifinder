"""Typed DataFrame - schema from Pydantic model or explicit column definitions."""
from __future__ import annotations

import pandas as pd

from enum import Enum
from pathlib import Path
from typing import Any, Callable, Dict, Generic, Iterator, List, Optional, Type, TypeVar, get_args

from amplifinder.data_types.records import Record, Schema
from amplifinder.data_types.validate_and_cast_df import validate_and_cast_df

T = TypeVar("T", bound=Record)
SelfTypedDF = TypeVar("SelfTypedDF", bound="TypedDF")
SelfRecordTypedDF = TypeVar("SelfRecordTypedDF", bound="RecordTypedDf")


def _clean_nan(v: Any) -> Any:
    """Convert NaN to None. Handles scalar values only; lists/arrays pass through."""
    try:
        return None if pd.isna(v) else v
    except (ValueError, TypeError):
        # pd.isna() fails on lists/arrays - just return as-is
        return v


def _serialize_for_csv(v: Any) -> Any:
    """Convert enums to values for CSV serialization (handles nested in lists/tuples/Records)."""
    from amplifinder.data_types.records import Record
    if isinstance(v, Enum):
        return v.value
    if isinstance(v, Record):
        return v.model_dump(mode='json')  # Convert Record to dict for CSV
    if isinstance(v, tuple):
        return list(_serialize_for_csv(x) for x in v)  # tuples as lists for CSV
    if isinstance(v, list):
        return [_serialize_for_csv(x) for x in v]
    if isinstance(v, dict):
        return {k: _serialize_for_csv(val) for k, val in v.items()}  # Recursively serialize dict values
    return v


class TypedDF:
    """Base class for typed DataFrames."""

    __slots__ = ("df", "schema", "headers", "_csv_field_formats")

    def __init__(self, df: pd.DataFrame, schema: Schema, headers: bool = True):
        self.df = df
        self.schema = schema
        self.headers = headers
        self._csv_field_formats = {}  # Dict[field_name, format_spec]

    def __len__(self) -> int:
        return len(self.df)

    def _has_meaningful_index(self) -> bool:
        """Check if DataFrame has a meaningful (non-default) index."""
        idx = self.df.index
        # Default RangeIndex starting at 0 is not meaningful
        if isinstance(idx, pd.RangeIndex) and idx.start == 0 and idx.step == 1:
            return False
        # Any other index (named, custom values, etc.) is meaningful
        return True

    def assert_matches(self, df: Optional[pd.DataFrame] = None, check_missing=True, check_extra=True) -> None:
        """Assert DataFrame matches schema exactly."""
        df = df if df is not None else self.df
        validate_and_cast_df(df, self.schema, check_missing=check_missing, check_extra=check_extra, cast=False)

    def to_csv(self, path: Path, index: Optional[bool] = None) -> None:
        """Save DataFrame to CSV. Auto-detects meaningful indices if index=None.

        Only serializes columns defined in schema.
        """
        from amplifinder.data_types.validate_and_cast_df import _is_optional

        # Filter DataFrame to schema columns only
        df = self.df[list(self.schema.column_names)].copy()

        # Convert Optional[int] columns to Int64
        dtypes = self.schema.dtypes
        for col in df.columns:
            if col in dtypes:
                col_dtype = dtypes[col].dtype
                # Check if it's Optional[int]
                if _is_optional(col_dtype):
                    inner_type = next((a for a in get_args(col_dtype) if a is not type(None)), None)
                    if inner_type is int:
                        # Convert to Int64 (nullable integer) to preserve integer format
                        df[col] = df[col].astype("Int64")

        # Serialize object columns (enums, Records, etc.)
        for col in df.columns:
            if df[col].dtype == object:
                df[col] = df[col].apply(_serialize_for_csv)

        # Format float columns: field-specific formats first, then default to 1 decimal
        for col in df.columns:
            if df[col].dtype in [float, 'float64', 'float32']:
                fmt = self._csv_field_formats.get(col, '.1f')  # Default: 1 decimal place
                df[col] = df[col].apply(lambda x: f"{x:{fmt}}" if pd.notna(x) else x)

        # Auto-detect if index should be saved

        if index is None:
            index = self._has_meaningful_index()

        df.to_csv(path, index=index, header=self.headers)

    @staticmethod
    def _read_csv_df(path: Path, schema: Schema, headers: bool, index_col: Optional[int] = None) -> pd.DataFrame:
        """Read CSV and validate against schema."""
        path = Path(path)
        if headers:
            df = pd.read_csv(path, index_col=index_col)
        else:
            df = pd.read_csv(path, header=None, names=schema.column_names, index_col=index_col)
        return validate_and_cast_df(df, schema, check_missing=True, check_extra=True, cast=True)

    @classmethod
    def from_csv(
        cls: Type[SelfTypedDF], path: Path, schema: Schema,
        headers: bool = True, index_col: Optional[int] = None
    ) -> SelfTypedDF:
        df = cls._read_csv_df(path, schema, headers, index_col=index_col)
        return cls(df, schema, headers)

    @classmethod
    def empty(cls: Type[SelfTypedDF], schema: Schema, headers: bool = True) -> SelfTypedDF:
        return cls(pd.DataFrame(columns=schema.column_names), schema, headers)

    def _construct(self: SelfTypedDF, df: pd.DataFrame) -> SelfTypedDF:
        """Construct new instance with same config but different df."""
        return type(self)(df, self.schema, self.headers)

    def pipe(self: SelfTypedDF, func: Callable[[pd.DataFrame], pd.DataFrame], inplace: bool = False) -> SelfTypedDF:
        """Pipe DataFrame through a function. If inplace, mutates self.df; else returns new instance."""
        new_df = validate_and_cast_df(func(self.df), self.schema, check_missing=True, check_extra=True, cast=True)
        if inplace:
            self.df = new_df
            return self
        return self._construct(new_df)

    def _row_to_dict(self, row: pd.Series) -> Dict[str, Any]:
        """Convert pandas Series row to cleaned dict (NaN → None)."""
        return {k: _clean_nan(v) for k, v in row.items()}

    def __getitem__(self, key: Any) -> Dict[str, Any]:
        """Access row by index key as dict."""
        row_series = self.df.loc[key]
        return self._row_to_dict(row_series)

    def items(self) -> Iterator[tuple[Any, Dict[str, Any]]]:
        """Iterate (key, row_dict) pairs."""
        for key, row_series in self.df.iterrows():
            yield key, self._row_to_dict(row_series)

    def keys(self) -> Iterator[Any]:
        """Iterate index keys."""
        return iter(self.df.index)

    def values(self) -> Iterator[Dict[str, Any]]:
        """Iterate rows as dicts (same as __iter__)."""
        return iter(self)

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        """Iterate rows as dicts with NaN → None."""
        for _, row in self.df.iterrows():
            yield self._row_to_dict(row)


class RecordTypedDf(TypedDF, Generic[T]):
    """DataFrame with Record schema. Iteration yields typed instances."""

    __slots__ = ("_record_type",)

    def __init__(self, df: pd.DataFrame, record_type: Type[T], headers: bool = True):
        self._record_type = record_type
        schema = record_type.schema()
        super().__init__(df, schema, headers)

    def _construct(self: SelfRecordTypedDF, df: pd.DataFrame) -> SelfRecordTypedDF:
        return type(self)(df, self._record_type, self.headers)

    @classmethod
    def from_records(cls: Type[SelfRecordTypedDF], records: List[T], record_type: Type[T]) -> SelfRecordTypedDF:
        schema = record_type.schema()
        if not records:
            return cls(pd.DataFrame(columns=schema.column_names), record_type)
        data = [{name: getattr(r, name) for name in record_type.model_fields} for r in records]
        return cls(pd.DataFrame(data), record_type)

    @classmethod
    def from_dict(cls: Type[SelfRecordTypedDF], records_dict: Dict[Any, T], record_type: Type[T]) -> SelfRecordTypedDF:
        """Create RecordTypedDf from dict[Any:Record] with keys as DataFrame index."""
        schema = record_type.schema()
        if not records_dict:
            return cls(pd.DataFrame(columns=schema.column_names), record_type)
        data = [r.model_dump() for r in records_dict.values()]
        df = pd.DataFrame(data, index=list(records_dict.keys()))
        return cls(df, record_type)

    @classmethod
    def from_csv(cls: Type[SelfRecordTypedDF], path: Path, record_type: Type[T],
                 headers: bool = True, index_col: Optional[int] = None) -> SelfRecordTypedDF:
        """Load from CSV. If index_col=0, reads first column as index. None means no index."""
        df = cls._read_csv_df(path, record_type.schema(), headers, index_col=index_col)
        return cls(df, record_type, headers)

    @classmethod
    def empty(cls: Type[SelfRecordTypedDF], record_type: Type[T], headers: bool = True) -> SelfRecordTypedDF:
        return cls(pd.DataFrame(columns=record_type.schema().column_names), record_type, headers)

    def to_records(self) -> List[T]:
        """Convert DataFrame to list of Record objects.

        Provides explicit round-trip: from_records() → to_records() should recreate same records.
        Uses existing __iter__() which calls model_validate() on each row.
        """
        return list(self)

    def to_dict(self) -> Dict[Any, T]:
        """Convert DataFrame to dict of Record objects using index as keys.

        Provides explicit round-trip: from_dict() → to_dict() should recreate same dict.
        Useful when DataFrame has meaningful index (e.g., tn_id).
        """
        return dict(self.items())

    def to_csv(self, path: Path, index: Optional[bool] = None) -> None:
        """Save DataFrame to CSV.

        If schema includes properties not in DataFrame, computes them first.
        Otherwise delegates to parent (which filters to schema columns).
        Applies field-specific formatting from record_type.CSV_FIELD_FORMATS if defined.
        """
        # Check if all schema columns are in DataFrame
        if not set(self.schema.column_names).issubset(set(self.df.columns)):
            # Some schema columns are properties - compute them via records
            records = self.to_records()
            filtered_data = [{col: getattr(record, col) for col in self.schema.column_names} for record in records]
            temp_df = TypedDF(pd.DataFrame(filtered_data), self.schema, self.headers)
            # Pass record type's formats to parent
            if hasattr(self._record_type, 'CSV_FIELD_FORMATS'):
                temp_df._csv_field_formats = self._record_type.CSV_FIELD_FORMATS
            temp_df.to_csv(path, index)
        else:
            # All schema columns are in DataFrame - parent filters to schema columns
            # Pass record type's formats to parent
            if hasattr(self._record_type, 'CSV_FIELD_FORMATS'):
                self._csv_field_formats = self._record_type.CSV_FIELD_FORMATS
            super().to_csv(path, index)

    def __getitem__(self, key: Any) -> T:
        """Access record by index key."""
        row_dict = super().__getitem__(key)
        return self._record_type.model_validate(row_dict)

    def items(self) -> Iterator[tuple[Any, T]]:
        """Iterate (key, record) pairs."""
        for key, row_dict in super().items():
            yield key, self._record_type.model_validate(row_dict)

    def __iter__(self) -> Iterator[T]:
        for row in super().__iter__():
            yield self._record_type.model_validate(row)

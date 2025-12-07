"""Typed DataFrame - schema from dataclass or explicit column definitions."""
from __future__ import annotations

import pandas as pd

from pathlib import Path
from typing import Any, Callable, Dict, Generic, Iterator, List, Mapping, Optional, Type, TypeVar

from amplifinder.data_types.records import Record, Schema
from amplifinder.data_types.validate_and_cast_df import validate_and_cast_df

T = TypeVar("T", bound=Record)
SelfTypedDF = TypeVar("SelfTypedDF", bound="TypedDF")
SelfRecordTypedDF = TypeVar("SelfRecordTypedDF", bound="RecordTypedDF")


def _clean_nan(v: Any) -> Any:
    """Convert NaN to None. Handles scalar values only; lists/arrays pass through."""
    try:
        return None if pd.isna(v) else v
    except (ValueError, TypeError):
        # pd.isna() fails on lists/arrays - just return as-is
        return v


class TypedDF:
    """Base class for typed DataFrames."""
    
    __slots__ = ("df", "schema", "headers")
    
    def __init__(self, df: pd.DataFrame, schema: Schema, headers: bool = True):
        self.df = df
        self.schema = schema
        self.headers = headers
    
    def __len__(self) -> int:
        return len(self.df)
    
    @property
    def is_empty(self) -> bool:
        return len(self.df) == 0
    
    def assert_matches(self, df: Optional[pd.DataFrame] = None, check_missing=True, check_extra=True) -> None:
        """Assert DataFrame matches schema exactly."""
        df = df if df is not None else self.df
        validate_and_cast_df(df, self.schema, check_missing=check_missing, check_extra=check_extra, cast=False)

    def replace_df(self, df: pd.DataFrame) -> None:
        """Validate and cast the DataFrame to the schema."""
        self.df = validate_and_cast_df(df, self.schema, check_missing=True, check_extra=True, cast=True)
    
    def to_csv(self, path: Path) -> None:
        self.df.to_csv(path, index=False, header=self.headers)
    
    @staticmethod
    def _read_csv_df(path: Path, schema: Schema, headers: bool) -> pd.DataFrame:
        """Read CSV and validate against schema."""
        path = Path(path)
        if headers:
            df = pd.read_csv(path)
        else:
            df = pd.read_csv(path, header=None, names=schema.column_names)
        return validate_and_cast_df(df, schema, check_missing=True, check_extra=True, cast=True)

    @classmethod
    def from_csv(cls: Type[SelfTypedDF], path: Path, schema: Schema, headers: bool = True) -> SelfTypedDF:
        df = cls._read_csv_df(path, schema, headers)
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
    
    def __iter__(self) -> Iterator[Dict[str, Any]]:
        """Iterate rows as dicts with NaN → None."""
        for _, row in self.df.iterrows():
            yield {k: _clean_nan(v) for k, v in row.items()}
        

class RecordTypedDF(TypedDF, Generic[T]):
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
        data = []
        for record in records:
            d = {col.name: getattr(record, col.name) for col in schema}
            d.update(record.extra)
            data.append(d)
        return cls(pd.DataFrame(data), record_type)
    
    @classmethod
    def from_csv(cls: Type[SelfRecordTypedDF], path: Path, record_type: Type[T], headers: bool = True) -> SelfRecordTypedDF:
        df = cls._read_csv_df(path, record_type.schema(), headers)
        return cls(df, record_type, headers)

    @classmethod
    def empty(cls: Type[SelfRecordTypedDF], record_type: Type[T], headers: bool = True) -> SelfRecordTypedDF:
        return cls(pd.DataFrame(columns=record_type.schema().column_names), record_type, headers)

    def __iter__(self) -> Iterator[T]:
        for row in super().__iter__():
            yield self._record_type(**row)

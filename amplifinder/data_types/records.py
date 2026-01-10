"""Base Record class using Pydantic."""
from __future__ import annotations

import pandas as pd

from pathlib import Path
from pydantic import BaseModel
from typing import Any, ClassVar, Dict, List, NamedTuple, Optional, Tuple, Type, TypeVar, get_origin, get_args, Union

T = TypeVar("T", bound="Record")


_TYPE_MAP = {
    "string": str,
    "double": float,
    "int32": int,
    "bool": bool,
    "cell": str,
}


class Column(NamedTuple):
    """Schema column definition."""
    name: str
    dtype: Type  # Full type hint (may be Optional[T])
    optional: bool = False

    @property
    def is_nullable(self) -> bool:
        """True if dtype is Optional[T] (Union with None)."""
        origin = get_origin(self.dtype)
        return origin is Union and type(None) in get_args(self.dtype)

    @property
    def base_dtype(self) -> Type:
        """Inner type for casting (unwraps Optional[T] → T)."""
        if self.is_nullable:
            inner = [t for t in get_args(self.dtype) if t is not type(None)]
            return inner[0] if inner else str
        return self.dtype


class Schema(NamedTuple):
    """Schema: Column definitions with allow_extra flag."""
    columns: Tuple[Column, ...]
    allow_extra: bool = True

    def __iter__(self):
        return iter(self.columns)

    def __len__(self):
        return len(self.columns)

    @property
    def column_names(self) -> Tuple[str, ...]:
        """All column names."""
        return tuple(col.name for col in self.columns)

    @property
    def required_columns(self) -> Tuple[str, ...]:
        """Names of required columns (optional=False)."""
        return tuple(col.name for col in self.columns if not col.optional)

    @property
    def dtypes(self) -> Dict[str, Column]:
        """Dictionary of column names to Column objects."""
        return {col.name: col for col in self.columns}

    @classmethod
    def from_csv(cls, path: Path, allow_extra: bool = True, *,
                 fields_col: str = "fields", types_col: str = "types",
                 optional_col: str = "optional") -> Schema:
        """Load schema from CSV file with columns: fields, types, optional."""
        df = pd.read_csv(path, skipinitialspace=True)
        return Schema(
            columns=tuple(
                Column(row[fields_col], _TYPE_MAP[row[types_col]], bool(row[optional_col]))
                for _, row in df.iterrows()
            ),
            allow_extra=allow_extra,
        )


class Record(BaseModel):
    """Base class for record models with schema support.

    Features:
        schema(): Returns Column definitions derived from model fields.
        from_other(): Create new record copying shared fields from another.
        extra: Dict for unknown fields (when model_config extra='allow').

    Optional columns:
        If a field's type is Optional[T], the column is marked as optional
        (may not exist in DataFrame) and nullable (can be None).

    Usage:
        class MyRecord(Record):
            name: str
            value: int
            notes: Optional[str] = None

        MyRecord.schema()
        # → (Column('name', str, False),
        #    Column('value', int, False),
        #    Column('notes', Optional[str], True))
    """

    NAME: ClassVar[str] = "records"  # Display name for logging

    # Class variable to track if extra fields allowed (for schema generation)
    ALLOW_EXTRA: ClassVar[bool] = False

    # CSV export control: None = export all fields, List[str] = export only specified fields/properties
    CSV_EXPORT_FIELDS: ClassVar[Optional[List[str]]] = None

    @classmethod
    def schema(cls) -> Schema:
        """Extract schema (Column definitions) for CSV read/write.

        If CSV_EXPORT_FIELDS is set, schema includes only those fields/properties.
        Otherwise includes all model fields.
        """
        # Determine which fields to include in schema
        if hasattr(cls, 'CSV_EXPORT_FIELDS') and cls.CSV_EXPORT_FIELDS is not None:
            field_names = cls.CSV_EXPORT_FIELDS
        else:
            field_names = list(cls.model_fields.keys())

        columns = []
        for name in field_names:
            if name in cls.model_fields:
                # It's a field - get type from model_fields
                field_info = cls.model_fields[name]
                dtype = field_info.annotation
                is_optional = (get_origin(dtype) is Union and type(None) in get_args(dtype))
            else:
                # It's a property - infer type from property return annotation
                from amplifinder.data_types.type_inference import infer_property_type
                dtype, is_optional = infer_property_type(cls, name)

            columns.append(Column(name=name, dtype=dtype, optional=is_optional))

        return Schema(columns=tuple(columns), allow_extra=cls.ALLOW_EXTRA)

    @classmethod
    def from_other(cls: type[T], obj: "Record", **kwargs) -> T:
        """Create instance from another record, copying shared fields."""
        shared = set(cls.model_fields.keys()) & set(type(obj).model_fields.keys())
        data = {name: getattr(obj, name) for name in shared}
        data.update(kwargs)
        return cls.model_validate(data)

    @property
    def extra(self) -> Dict[str, Any]:
        """Return extra fields not in schema."""
        return self.__pydantic_extra__ or {}

    def __repr__(self) -> str:
        """Custom repr that excludes sequence fields."""
        fields = []
        for field_name in self.model_fields.keys():
            # Skip sequence fields
            value = getattr(self, field_name, None)
            if ('sequence' in field_name.lower() or 'seq' in field_name.lower()) and isinstance(value, str):
                s = f"{field_name}='<{len(value)} chars>'"
            else:
                s = f"{field_name}={repr(value)}"
            fields.append(s)

        return f"{self.__class__.__name__}({', '.join(fields)})"

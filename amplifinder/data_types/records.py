"""Base Record class and schema types."""
from __future__ import annotations

import pandas as pd

from dataclasses import dataclass, field, fields as dataclass_fields
from pathlib import Path
from typing import Any, ClassVar, Dict, NamedTuple, Tuple, Type, TypeVar, Union, get_type_hints, get_origin, get_args

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


@dataclass(frozen=True)
class Schema:
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
    def optional_columns(self) -> Tuple[str, ...]:
        """Names of optional columns (optional=True)."""
        return tuple(col.name for col in self.columns if col.optional)

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


@dataclass(kw_only=True)
class Record:
    """Base class for record dataclasses with schema support.

    Features:
        schema(): Returns Column definitions derived from dataclass fields.
        from_instance(): Create new record copying shared fields from another.
        extra: Dict for unknown columns (when ALLOW_EXTRA=True).
        ALLOW_EXTRA: ClassVar controlling whether unknown columns are allowed.

    Optional columns:
        If a field's type is Optional[T], the column is marked as optional
        (may not exist in DataFrame) and nullable (can be None).

    Usage:
        @dataclass(kw_only=True)
        class MyRecord(Record):
            name: str
            value: int
            notes: Optional[str] = None

        MyRecord.schema()
        # → (Column('name', str, False),
        #    Column('value', int, False),
        #    Column('notes', Optional[str], True))

        # Column properties:
        # col.dtype       → Optional[str]  (full type)
        # col.base_dtype  → str            (for casting)
        # col.is_nullable → True           (can be None)
    """

    ALLOW_EXTRA: ClassVar[bool] = True
    extra: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def schema(cls) -> Schema:
        """Extract schema (Column definitions) from this dataclass."""
        hints = get_type_hints(cls)
        columns = []

        for f in dataclass_fields(cls):
            if f.name == "extra":
                continue

            dtype = hints.get(f.name, str)
            # Check if type is Optional (Union with None)
            is_optional = (
                get_origin(dtype) is Union and
                type(None) in get_args(dtype)
            )

            columns.append(Column(
                name=f.name,
                dtype=dtype,
                optional=is_optional,
                )
            )

        return Schema(columns=tuple(columns), allow_extra=cls.ALLOW_EXTRA)

    @classmethod
    def from_dict(cls: type[T], data: Dict[str, Any]) -> T:
        """Create instance from dict, separating schema fields from extras."""
        schema_cols = {col.name for col in cls.schema()}
        schema_fields = {k: v for k, v in data.items() if k in schema_cols}
        extra_fields = {k: v for k, v in data.items() if k not in schema_cols}
        if extra_fields and cls.ALLOW_EXTRA:
            schema_fields["extra"] = extra_fields
        return cls(**schema_fields)

    @classmethod
    def from_other(cls: type[T], obj: Record, **kwargs) -> T:
        """Create instance from another record, copying shared fields."""
        shared = {f.name for f in dataclass_fields(cls)} & {f.name for f in dataclass_fields(obj)}
        data = {name: getattr(obj, name) for name in shared if name != "extra"}
        data.update(kwargs)
        return cls(**data)

"""Validate and cast DataFrame columns using Pydantic TypeAdapter."""
import ast
import pandas as pd
import numpy as np
from enum import Enum
from typing import Any, Type, Union, get_origin, get_args, NamedTuple

from pydantic import TypeAdapter

from amplifinder.data_types.records import Schema


# Cache TypeAdapters for performance
_adapter_cache: dict[Type, TypeAdapter] = {}


def _get_adapter(dtype: Type) -> TypeAdapter:
    """Get or create a TypeAdapter for a type."""
    if dtype not in _adapter_cache:
        _adapter_cache[dtype] = TypeAdapter(dtype)
    return _adapter_cache[dtype]


def _is_optional(t: Type) -> bool:
    return get_origin(t) is Union and type(None) in get_args(t)


def _unwrap_optional(t: Type) -> Type:
    """Unwrap Optional[T] -> T."""
    if _is_optional(t):
        return next(a for a in get_args(t) if a is not type(None))
    return t


def _is_compound(t: Type) -> bool:
    """Check if type is a compound type (list, dict, tuple)."""
    origin = get_origin(_unwrap_optional(t))
    return origin in (list, dict, tuple)


def _is_enum(t: Type) -> bool:
    """Check if type is an Enum subclass."""
    try:
        return isinstance(t, type) and issubclass(t, Enum)
    except TypeError:
        return False


def _is_namedtuple(t: Type) -> bool:
    """Check if type is a NamedTuple subclass."""
    try:
        return isinstance(t, type) and issubclass(t, tuple) and hasattr(t, '_fields')
    except TypeError:
        return False


def parse_compound(value: Any, expected_type: Type) -> Any:
    """Parse compound type using Pydantic TypeAdapter."""
    # Handle NaN
    if isinstance(value, (list, dict, tuple)):
        pass  # Already parsed
    elif pd.isna(value):
        if not _is_optional(expected_type):
            raise TypeError(f"NaN value for type {expected_type} which is not optional")
        return value
    elif isinstance(value, str):
        value = ast.literal_eval(value)

    # Validate/coerce with Pydantic
    adapter = _get_adapter(expected_type)
    return adapter.validate_python(value)


def validate_and_cast_df(
    df: pd.DataFrame,
    columns: Schema,
    check_missing: bool = True,
    check_extra: bool = False,
    cast: bool = True,
) -> pd.DataFrame:
    """Validate and cast DataFrame against schema using Pydantic."""

    if check_missing:
        missing = set(columns.required_columns) - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

    if check_extra:
        extra = set(df.columns) - set(columns.column_names)
        if extra:
            raise ValueError(f"Extra columns: {extra}")

    dtypes = columns.dtypes
    scalar_cols = {}

    for col_name in df.columns:
        if col_name not in dtypes:
            continue
        _dtype = dtypes[col_name]
        base_type = _dtype.base_dtype  # Unwraps Optional[T] → T

        if cast:
            if _is_compound(base_type) or _is_compound(_dtype.dtype) or _is_namedtuple(base_type):
                # Compound types and NamedTuple: parse with Pydantic
                # (NamedTuple isn't detected as compound by get_origin, so check separately)
                df = df.copy()
                df[col_name] = df[col_name].apply(lambda x, et=_dtype.dtype: parse_compound(x, et))
            elif _is_enum(base_type):
                # Enum: convert values via Pydantic
                df = df.copy()
                adapter = _get_adapter(base_type)
                df[col_name] = df[col_name].apply(
                    lambda x, a=adapter: a.validate_python(x) if pd.notna(x) else x
                )
            else:
                scalar_cols[col_name] = base_type
        else:
            expected = _dtype.dtype
            actual = df[col_name].dtype
            expected_dtype = np.dtype('O') if expected is str else np.dtype(expected)
            if actual != expected_dtype:
                raise TypeError(f"Column {col_name}: expected {expected_dtype}, got {actual}")

    if cast and scalar_cols:
        df = df.astype(scalar_cols)

    return df

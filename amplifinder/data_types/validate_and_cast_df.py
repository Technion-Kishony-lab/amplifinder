import pandas as pd
import numpy as np
from typing import Any, Type, Union, get_origin, get_args
import ast

from amplifinder.data_types.records import Schema


def _is_optional(t: Type) -> bool:
    return get_origin(t) is Union and type(None) in get_args(t)


def _unwrap_optional(t: Type) -> Type:
    """Unwrap Optional[T] -> T."""
    if _is_optional(t):
        return next(a for a in get_args(t) if a is not type(None))
    return t


def parse_compound(value: Any, expected_type: Type) -> Any:
    """Parse compound type (list/dict/tuple) from string."""
    if pd.isna(value):
        if not _is_optional(expected_type):
            raise TypeError(f"NaN value for type {expected_type} which is not optional")
        return value
    if isinstance(value, str):
        value = ast.literal_eval(value)
    
    inner_type = _unwrap_optional(expected_type)
    origin = get_origin(inner_type)
    args = get_args(inner_type)
    
    if origin is list and not isinstance(value, list):
        raise TypeError(f"Expected list, got {type(value).__name__}")
    elif origin is dict and not isinstance(value, dict):
        raise TypeError(f"Expected dict, got {type(value).__name__}")
    elif origin is tuple and not isinstance(value, tuple):
        raise TypeError(f"Expected tuple, got {type(value).__name__}")
    
    if args and value:
        if origin is list:
            inner = args[0]
            for i, item in enumerate(value):
                if not isinstance(item, inner):
                    raise TypeError(f"List[{i}]: expected {inner.__name__}, got {type(item).__name__}")
        elif origin is tuple:
            for i, (item, inner) in enumerate(zip(value, args)):
                if not isinstance(item, inner):
                    raise TypeError(f"Tuple[{i}]: expected {inner.__name__}, got {type(item).__name__}")
        elif origin is dict:
            ktype, vtype = args
            for k, v in value.items():
                if not isinstance(k, ktype):
                    raise TypeError(f"Dict key {k!r}: expected {ktype.__name__}")
                if not isinstance(v, vtype):
                    raise TypeError(f"Dict[{k!r}]: expected {vtype.__name__}")
    
    return value


def validate_and_cast_df(
    df: pd.DataFrame,
    columns: Schema,
    check_missing: bool = True,
    check_extra: bool = False,
    cast: bool = True,
) -> pd.DataFrame:
    """Validate and cast DataFrame against schema."""
    
    if check_missing:
        missing = set(columns.required_columns) - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
    
    if check_extra:
        extra = set(df.columns) - set(columns.column_names)
        if extra:
            raise ValueError(f"Extra columns: {extra}")
    
    dtypes = columns.dtypes
    compound_origins = (list, dict, tuple)
    scalar_cols = {}
    
    for col_name in df.columns:
        if col_name not in dtypes:
            continue
        _dtype = dtypes[col_name]
        base_type = _dtype.base_dtype  # Unwraps Optional[T] → T
        
        if cast:
            if get_origin(base_type) in compound_origins:
                df = df.copy()
                df[col_name] = df[col_name].apply(lambda x, et=base_type: parse_compound(x, et))
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

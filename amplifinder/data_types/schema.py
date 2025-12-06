"""DataFrame schema definitions with types and I/O."""

import ast
from dataclasses import dataclass, fields as dataclass_fields, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Type, TypeVar, Union, get_type_hints, get_origin, get_args

import numpy as np
import pandas as pd


# Map CSV type names to Python types
_TYPE_MAP = {
    "string": str,
    "double": float,
    "int32": int,
    "bool": bool,
    "cell": str,  # MATLAB cell arrays → strings
}

# Map generic origin to expected Python type
_ORIGIN_TO_TYPE = {
    list: list,
    dict: dict,
    tuple: tuple,
}


def _parse_compound(value: Any, expected_type: Type) -> Any:
    """Parse compound type from string and validate."""
    if pd.isna(value):
        return value

    # Parse string representation
    if isinstance(value, str):
        value = ast.literal_eval(value)

    # Check container type
    origin = get_origin(expected_type)
    expected_container = _ORIGIN_TO_TYPE.get(origin)
    if expected_container and not isinstance(value, expected_container):
        raise TypeError(f"Expected {expected_container.__name__}, got {type(value).__name__}")

    # Check inner types
    args = get_args(expected_type)
    if args and value:
        if origin is list:
            inner_type = args[0]
            for i, item in enumerate(value):
                if not isinstance(item, inner_type):
                    raise TypeError(f"List item {i}: expected {inner_type.__name__}, got {type(item).__name__}")
        elif origin is tuple:
            # Tuple[int, int] means fixed types per position
            for i, (item, inner_type) in enumerate(zip(value, args)):
                if not isinstance(item, inner_type):
                    raise TypeError(f"Tuple item {i}: expected {inner_type.__name__}, got {type(item).__name__}")
        elif origin is dict:
            key_type, val_type = args[0], args[1]
            for k, v in value.items():
                if not isinstance(k, key_type):
                    raise TypeError(f"Dict key {k!r}: expected {key_type.__name__}, got {type(k).__name__}")
                if not isinstance(v, val_type):
                    raise TypeError(f"Dict value for {k!r}: expected {val_type.__name__}, got {type(v).__name__}")

    return value

T = TypeVar("T")


@dataclass
class TypedSchema:
    """Schema with column names, dtypes, and CSV I/O.

    Can be created from CSV file (from_csv) or from a dataclass (from_dataclass).
    """

    columns: Tuple[Tuple[str, Any, bool], ...]  # (name, dtype, optional) tuples
    dataclass_type: Optional[Type] = field(default=None, repr=False)

    @classmethod
    def from_csv(cls, path: Path) -> "TypedSchema":
        """Load schema definition from CSV file.

        CSV must have columns: fields, types, optional
        """
        df = pd.read_csv(path, skipinitialspace=True)
        columns = tuple(
            (row["fields"], _TYPE_MAP[row["types"]], bool(row["optional"]))
            for _, row in df.iterrows()
        )
        return cls(columns)

    @classmethod
    def from_dataclass(cls, dc: Type[T]) -> "TypedSchema":
        """Create schema from a dataclass.

        - Field metadata {"optional": True} marks column as optional (may not exist)
        - Optional[T] type hint means value can be None (column still required)
        - Field named 'extra' (Dict[str, Any]) is skipped (catch-all for unknown columns)
        """
        hints = get_type_hints(dc)
        columns = []

        for f in dataclass_fields(dc):
            if f.name == "extra":
                continue  # skip the catch-all field

            hint = hints.get(f.name, str)
            origin = get_origin(hint)

            # Extract actual type from Optional[T] if present
            if origin is Union and type(None) in get_args(hint):
                inner_types = [t for t in get_args(hint) if t is not type(None)]
                dtype = inner_types[0] if inner_types else str
            else:
                dtype = hint

            # Column optional = field metadata, NOT type hint
            is_optional = f.metadata.get("optional", False) if f.metadata else False

            columns.append((f.name, dtype, is_optional))

        return cls(columns=tuple(columns), dataclass_type=dc)

    @property
    def names(self) -> Tuple[str, ...]:
        return tuple(c[0] for c in self.columns)

    @property
    def dtypes(self) -> Dict[str, Any]:
        return {c[0]: c[1] for c in self.columns}

    @property
    def required_names(self) -> Tuple[str, ...]:
        return tuple(c[0] for c in self.columns if not c[2])

    @property
    def optional_names(self) -> Tuple[str, ...]:
        return tuple(c[0] for c in self.columns if c[2])

    def empty(self) -> pd.DataFrame:
        """Empty DataFrame with correct columns and types."""
        df = pd.DataFrame(columns=self.names)
        # Only cast scalar types; compound types (list, dict, tuple) stay as object
        compound_origins = (list, dict, tuple)
        scalar_dtypes = {c: t for c, t in self.dtypes.items() if get_origin(t) not in compound_origins}
        if scalar_dtypes:
            df = df.astype(scalar_dtypes)
        return df

    def create(self, data: Dict[str, Any]) -> pd.DataFrame:
        """Create DataFrame from column dict with schema validation and types."""
        df = pd.DataFrame(data)
        return self._validate(df)

    def from_records(self, records: List[Dict[str, Any]]) -> pd.DataFrame:
        """Create DataFrame from list of row dicts with schema validation and types."""
        if not records:
            return self.empty()
        df = pd.DataFrame.from_records(records)
        return self._validate(df)

    def _validate(self, df: pd.DataFrame, strict: bool = False, cast: bool = True) -> pd.DataFrame:
        """Validate columns and optionally cast to schema dtypes.

        Args:
            df: DataFrame to validate
            strict: If True, reject columns not in schema
            cast: If True, cast to types; if False, assert types match
        """
        # Required columns must be present
        missing = set(self.required_names) - set(df.columns)
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

        # Extra columns: reject if strict, otherwise keep as-is
        extra = set(df.columns) - set(self.names)
        if strict and extra:
            raise ValueError(f"Extra columns: {extra}")

        schema_cols = [c for c in self.names if c in df.columns]
        if cast:
            # Separate compound types (list, dict, tuple) from scalar columns
            compound_origins = (list, dict, tuple)
            compound_cols = {c for c in schema_cols if get_origin(self.dtypes[c]) in compound_origins}
            scalar_cols = {c: self.dtypes[c] for c in schema_cols if c not in compound_cols}

            # Parse compound columns from string representation
            for col in compound_cols:
                expected_type = self.dtypes[col]
                df[col] = df[col].apply(
                    lambda x, et=expected_type: _parse_compound(x, et)
                )

            # Cast scalar columns
            if scalar_cols:
                df = df.astype(scalar_cols)
        else:
            # Assert types match
            for col in schema_cols:
                expected = self.dtypes[col]
                actual = df[col].dtype
                expected_dtype = np.dtype('O') if expected is str else np.dtype(expected)
                if actual != expected_dtype:
                    raise TypeError(f"Column {col}: expected {expected_dtype}, got {actual}")
        return df

    def read_csv(self, path: Path, headers: bool = True, strict: bool = False) -> pd.DataFrame:
        """Read CSV, validate and cast to schema.

        Args:
            path: CSV file path
            headers: If True, CSV has header row; if False, use schema column names
            strict: If True and headers=True, reject extra columns
        """
        if headers:
            df = pd.read_csv(path)
            return self._validate(df, strict)
        else:
            # For headerless CSVs, empty file means no results
            if path.stat().st_size == 0:
                return self.empty()
            df = pd.read_csv(path, header=None, names=self.names)
            return df.astype(self.dtypes)

    def to_csv(self, data: Union[pd.DataFrame, List[T]], path: Path, headers: bool = True) -> None:
        """Save DataFrame or list of dataclass instances to CSV."""
        if isinstance(data, list):
            data = self.instances_to_df(data)
        data.to_csv(path, index=False, header=headers)

    def assert_matches(self, df: pd.DataFrame) -> None:
        """Assert DataFrame matches schema columns and dtypes exactly."""
        self._validate(df, strict=True, cast=False)

    def df_to_instances(self, df: pd.DataFrame) -> List[T]:
        """Convert DataFrame to list of dataclass instances.

        Requires schema to be created from a dataclass (from_dataclass).
        Unknown columns are stored in the 'extra' field if present.
        """
        if self.dataclass_type is None:
            raise ValueError("df_to_instances requires schema created with from_dataclass()")

        if df.empty:
            return []

        known_cols = set(self.names)
        has_extra = any(f.name == "extra" for f in dataclass_fields(self.dataclass_type))
        instances = []

        for _, row in df.iterrows():
            kwargs = {}
            extra = {}

            for col, val in row.items():
                # Handle NaN → None
                if pd.isna(val):
                    val = None

                if col in known_cols:
                    kwargs[col] = val
                elif has_extra:
                    extra[col] = val
                else:
                    raise ValueError(f"Unknown column: {col}")

            if has_extra:
                kwargs["extra"] = extra

            instances.append(self.dataclass_type(**kwargs))

        return instances

    def instances_to_df(self, instances: List[T]) -> pd.DataFrame:
        """Convert list of dataclass instances to DataFrame.

        Merges 'extra' field contents into the DataFrame columns.
        """
        if not instances:
            return self.empty()

        records = []
        for inst in instances:
            d = {}
            for f in dataclass_fields(inst):
                if f.name == "extra":
                    extra = getattr(inst, "extra", None)
                    if extra:
                        d.update(extra)
                else:
                    d[f.name] = getattr(inst, f.name)
            records.append(d)

        return pd.DataFrame(records)

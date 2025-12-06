"""DataFrame schema definitions with types and I/O."""

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
        return df.astype(self.dtypes)

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
            dtypes_to_apply = {c: self.dtypes[c] for c in schema_cols}
            df = df.astype(dtypes_to_apply)
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

    def to_csv(self, df: pd.DataFrame, path: Path, headers: bool = True) -> None:
        """Save DataFrame to CSV."""
        df.to_csv(path, index=False, header=headers)

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

"""DataFrame schema definitions with types and I/O."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

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


@dataclass(frozen=True)
class TypedSchema:
    """Schema with column names, dtypes, and CSV I/O."""

    columns: Tuple[Tuple[str, Any, bool], ...]  # (name, dtype, optional) tuples

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

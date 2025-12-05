"""DataFrame schema definitions with types and I/O."""

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd


@dataclass(frozen=True)
class TypedSchema:
    """Schema with column names, dtypes, and CSV I/O."""
    
    columns: Tuple[Tuple[str, Any], ...]  # (name, dtype) pairs
    
    @property
    def names(self) -> Tuple[str, ...]:
        return tuple(c[0] for c in self.columns)
    
    @property
    def dtypes(self) -> Dict[str, Any]:
        return {c[0]: c[1] for c in self.columns}
    
    def empty(self) -> pd.DataFrame:
        """Empty DataFrame with correct columns and types."""
        df = pd.DataFrame(columns=self.names)
        return df.astype(self.dtypes)
    
    def create(self, data: Dict[str, Any]) -> pd.DataFrame:
        """Create DataFrame from column dict with schema validation and types."""
        df = pd.DataFrame(data)
        return self._validate_and_cast(df)
    
    def from_records(self, records: List[Dict[str, Any]]) -> pd.DataFrame:
        """Create DataFrame from list of row dicts with schema validation and types."""
        if not records:
            return self.empty()
        df = pd.DataFrame.from_records(records)
        return self._validate_and_cast(df)
    
    def _validate_and_cast(self, df: pd.DataFrame, strict: bool = False) -> pd.DataFrame:
        """Validate columns and cast to schema dtypes."""
        missing = set(self.names) - set(df.columns)
        if missing:
            raise ValueError(f"CSV missing columns: {missing}")
        extra = set(df.columns) - set(self.names)
        if strict and extra:
            raise ValueError(f"CSV has extra columns: {extra}")
        # Select schema columns in order and cast types
        return df[list(self.names)].astype(self.dtypes)
    
    def read_csv(self, path: Path, headers: bool = True, strict: bool = False) -> pd.DataFrame:
        """Read CSV, validate and cast to schema.
        
        Args:
            path: CSV file path
            headers: If True, CSV has header row; if False, use schema column names
            strict: If True and headers=True, reject extra columns
        """
        if headers:
            df = pd.read_csv(path)
            return self._validate_and_cast(df, strict)
        else:
            df = pd.read_csv(path, header=None, names=self.names)
            return df.astype(self.dtypes)
    
    def to_csv(self, df: pd.DataFrame, path: Path, headers: bool = True) -> None:
        """Save DataFrame to CSV."""
        df.to_csv(path, index=False, header=headers)


# BLAST output format 10 (CSV) columns
BLAST_SCHEMA = TypedSchema((
    ("query", str),
    ("subject", str),
    ("percent_identical", float),
    ("length", int),
    ("mismatch", int),
    ("gapopen", int),
    ("qstart", int),
    ("qend", int),
    ("sstart", int),
    ("send", int),
    ("evalue", float),
    ("bitscore", float),
))

# IS element locations table
IS_LOC_SCHEMA = TypedSchema((
    ("ID", int),
    ("IS_Name", str),
    ("IS_scaf", str),
    ("LocLeft", int),
    ("LocRight", int),
    ("Complement", bool),
    ("Join", bool),
))

"""Data loading utilities."""

from functools import lru_cache
from pathlib import Path
from typing import Dict

import pandas as pd

# Map MATLAB types to Python types
TYPE_MAP = {
    "string": "string",
    "double": "float",
    "int32": "int",
    "cell": "string",  # treat cell as string for parsing
}

RECORD_TYPES = ["JC", "SNP", "MOB", "DEL", "UN"]


def get_data_path() -> Path:
    """Get path to bundled data directory."""
    return Path(__file__).parent


@lru_cache(maxsize=None)
def load_field_defs(name: str) -> pd.DataFrame:
    """Load field definitions from CSV.
    
    Args:
        name: Field type name (JC, SNP, MOB, DEL, UN)
        
    Returns:
        DataFrame with columns: fields, types, optional
    """
    csv_path = get_data_path() / "fields" / f"{name}_fields.csv"
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df["types"] = df["types"].map(TYPE_MAP).fillna("string")
    df["optional"] = df["optional"].astype(bool)
    return df


def load_all_field_defs() -> Dict[str, pd.DataFrame]:
    """Load all field definitions.
    
    Returns:
        Dict mapping record type to field definitions DataFrame
    """
    return {name: load_field_defs(name) for name in RECORD_TYPES}


def get_isfinder_db_path() -> Path:
    """Get path to bundled ISfinder database (IS.fna)."""
    return get_data_path() / "ISfinderDB" / "IS.fna"

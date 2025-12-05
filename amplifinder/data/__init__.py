"""Data loading utilities."""

from functools import lru_cache
from pathlib import Path
from typing import Dict, Tuple

# Map MATLAB types to Python types
TYPE_MAP = {
    "string": "string",
    "double": "float",
    "int32": "int",
    "cell": "string",  # treat cell as string for parsing
}


def get_data_path() -> Path:
    """Get path to bundled data directory."""
    return Path(__file__).parent


@lru_cache(maxsize=None)
def load_field_defs(name: str) -> Dict[str, Tuple[str, bool]]:
    """Load field definitions from CSV.
    
    Args:
        name: Field type name (JC, SNP, MOB, DEL, UN)
        
    Returns:
        Dict mapping field name to (python_type, is_optional)
    """
    csv_path = get_data_path() / "fields" / f"{name}_fields.csv"
    
    if not csv_path.exists():
        raise FileNotFoundError(f"Field definitions not found: {csv_path}")
    
    fields = {}
    with open(csv_path) as f:
        # Skip header
        next(f)
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            if len(parts) >= 3:
                name_col, type_col, optional_col = parts[0], parts[1], parts[2]
                py_type = TYPE_MAP.get(type_col, "string")
                is_optional = optional_col == "1"
                fields[name_col] = (py_type, is_optional)
    
    return fields


# Convenience loaders for each type
def load_jc_fields() -> Dict[str, Tuple[str, bool]]:
    """Load JC (junction) field definitions."""
    return load_field_defs("JC")


def load_snp_fields() -> Dict[str, Tuple[str, bool]]:
    """Load SNP field definitions."""
    return load_field_defs("SNP")


def load_mob_fields() -> Dict[str, Tuple[str, bool]]:
    """Load MOB (mobile element) field definitions."""
    return load_field_defs("MOB")


def load_del_fields() -> Dict[str, Tuple[str, bool]]:
    """Load DEL (deletion) field definitions."""
    return load_field_defs("DEL")


def load_un_fields() -> Dict[str, Tuple[str, bool]]:
    """Load UN (unknown) field definitions."""
    return load_field_defs("UN")


def load_all_field_defs() -> Dict[str, Dict[str, Tuple[str, bool]]]:
    """Load all field definitions.
    
    Returns:
        Dict mapping record type to field definitions
    """
    return {
        "JC": load_jc_fields(),
        "SNP": load_snp_fields(),
        "MOB": load_mob_fields(),
        "DEL": load_del_fields(),
        "UN": load_un_fields(),
    }


def get_isfinder_db_path() -> Path:
    """Get path to bundled ISfinder database (IS.fna)."""
    return get_data_path() / "ISfinderDB" / "IS.fna"

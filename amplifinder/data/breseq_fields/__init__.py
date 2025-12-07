"""Breseq record type field definitions."""

from pathlib import Path
from typing import Dict

from amplifinder.data_types.records import Schema

RECORD_TYPES = ["JC", "SNP", "MOB", "DEL", "UN"]

_FIELDS_PATH = Path(__file__).parent


def load_all_field_defs() -> Dict[str, Schema]:
    """Load all field definitions as Schema objects."""
    return {name: Schema.from_csv(_FIELDS_PATH / f"{name}_fields.csv")
            for name in RECORD_TYPES}

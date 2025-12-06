"""Breseq record type field definitions."""

from pathlib import Path
from typing import Dict

from amplifinder.data_types.schema import TypedSchema

RECORD_TYPES = ["JC", "SNP", "MOB", "DEL", "UN"]

_FIELDS_PATH = Path(__file__).parent


def load_all_field_defs() -> Dict[str, TypedSchema]:
    """Load all field definitions."""
    return {name: TypedSchema.from_csv(_FIELDS_PATH / f"{name}_fields.csv")
            for name in RECORD_TYPES}

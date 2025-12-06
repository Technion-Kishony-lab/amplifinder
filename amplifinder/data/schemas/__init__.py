"""Schema definitions loaded from CSV files."""

from pathlib import Path

from amplifinder.data_types.schema import TypedSchema

_SCHEMAS_PATH = Path(__file__).parent

BLAST_SCHEMA = TypedSchema.from_csv(_SCHEMAS_PATH / "blast.csv")
IS_LOC_SCHEMA = TypedSchema.from_csv(_SCHEMAS_PATH / "is_loc.csv")

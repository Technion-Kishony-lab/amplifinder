"""Schema definitions loaded from CSV files."""

from pathlib import Path

from amplifinder.data_types.schema import TypedSchema

_SCHEMAS_PATH = Path(__file__).parent

BLAST_SCHEMA = TypedSchema.from_csv(_SCHEMAS_PATH / "blast.csv")
TN_LOC_SCHEMA = TypedSchema.from_csv(_SCHEMAS_PATH / "tn_loc.csv")

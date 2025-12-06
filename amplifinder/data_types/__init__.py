"""Data types and structures."""

from amplifinder.data_types.schema import TypedSchema
from amplifinder.data_types.tabularable import Tabularable
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.junction import Junction
from amplifinder.data.schemas import BLAST_SCHEMA, IS_LOC_SCHEMA

__all__ = [
    "TypedSchema",
    "BLAST_SCHEMA",
    "IS_LOC_SCHEMA",
    "Tabularable",
    "Genome",
    "GenomeRegistry",
    "get_genome",
    "Junction",
]

"""Data types and structures."""

from amplifinder.data_types.schema import TypedSchema
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.records import TnLoc, BlastHit, TN_LOC_SCHEMA, BLAST_SCHEMA

__all__ = [
    "TypedSchema",
    "TnLoc",
    "BlastHit",
    "TN_LOC_SCHEMA",
    "BLAST_SCHEMA",
    "Genome",
    "GenomeRegistry",
    "get_genome",
]

"""Data types and structures."""

from amplifinder.data_types.typed_df import RecordTypedDF, TypedDF
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.records import Column, Schema, Record
from amplifinder.data_types.record_types import (
    Side, Orientation,
    TnLoc, TnEndSeq, BlastHit, Junction, RefTnJunction, TnJunction, TnJc2, TnMatch,
    Coverage, JunctionCoverage, CoveredTnJc2, RawEvent, ClassifiedTnJc2,
    CandidateTnJc2, JunctionType, EventModifier, AnalyzedTnJc2,
)

__all__ = [
    "Side",
    "Orientation",
    "Column",
    "Schema",
    "RecordTypedDF",
    "TypedDF",
    "Record",
    "TnLoc",
    "TnEndSeq",
    "BlastHit",
    "Junction",
    "RefTnJunction",
    "TnJunction",
    "TnJc2",
    "TnMatch",
    "Genome",
    "GenomeRegistry",
    "get_genome",
    # Coverage types
    "Coverage",
    "JunctionCoverage",
    "CoveredTnJc2",
    "RawEvent",
    "ClassifiedTnJc2",
    "CandidateTnJc2",
    "JunctionType",
    "EventModifier",
    "AnalyzedTnJc2",
]

"""Data types and structures."""

from amplifinder.data_types.typed_df import RecordTypedDF, TypedDF
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.records import Column, Schema, Record
from amplifinder.data_types.record_types import (
    Side, Orientation,
    RefTnLoc, SeqRefTnSide, BlastHit, Junction, RefTnJunction, TnJunction, TnJc2, RefTnSide,
    Coverage, JunctionCoverage, CoveredTnJc2, RawEvent, ClassifiedTnJc2,
    CandidateTnJc2, JunctionType, EventModifier, AnalyzedTnJc2, ISJC2Export,
)

__all__ = [
    "Side",
    "Orientation",
    "Column",
    "Schema",
    "RecordTypedDF",
    "TypedDF",
    "Record",
    "RefTnLoc",
    "SeqRefTnSide",
    "BlastHit",
    "Junction",
    "RefTnJunction",
    "TnJunction",
    "TnJc2",
    "RefTnSide",
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
    "ISJC2Export",
]

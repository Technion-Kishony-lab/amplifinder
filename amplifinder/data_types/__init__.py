"""Data types and structures."""

from amplifinder.data_types.typed_df import RecordTypedDf, TypedDF
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.records import Column, Schema, Record
from amplifinder.data_types.record_types import (
    Side, Orientation,
    RefTnLoc, BlastHit, Junction, RefTnJunction, TnJunction, RawTnJc2, RefTnSide,
    Average, JunctionCoverage, CoveredTnJc2, RawEvent, ClassifiedTnJc2,
    FilteredTnJc2, JunctionType, EventModifier, AnalyzedTnJc2, ExportedTnJc2,
)

__all__ = [
    "Side",
    "Orientation",
    "Column",
    "Schema",
    "RecordTypedDf",
    "TypedDF",
    "Record",
    "RefTnLoc",
    "BlastHit",
    "Junction",
    "RefTnJunction",
    "TnJunction",
    "RawTnJc2",
    "RefTnSide",
    "Genome",
    "GenomeRegistry",
    "get_genome",
    # Coverage types
    "Average",
    "JunctionCoverage",
    "CoveredTnJc2",
    "RawEvent",
    "ClassifiedTnJc2",
    "FilteredTnJc2",
    "JunctionType",
    "EventModifier",
    "AnalyzedTnJc2",
    "ExportedTnJc2",
]

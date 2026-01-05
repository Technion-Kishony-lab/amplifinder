"""Data types and structures."""

from amplifinder.data_types.typed_df import RecordTypedDf, TypedDF
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.records import Column, Schema, Record
from amplifinder.data_types.enums import Side, Orientation, AverageMethod, JunctionType, EventModifier, JunctionCoverage
from amplifinder.data_types.scaffold import JcArm
from amplifinder.data_types.record_types import (
    RefTn, BlastHit, Junction, BreseqJunction, RefTnJunction, TnJunction, RawTnJc2, RefTnSide, OffsetRefTnSide,
    CoveredTnJc2, RawEvent, ClassifiedTnJc2,
    FilteredTnJc2, AnalyzedTnJc2, ExportedTnJc2,
)

__all__ = [
    "Side",
    "Orientation",
    "AverageMethod",
    "Column",
    "Schema",
    "RecordTypedDf",
    "TypedDF",
    "Record",
    "RefTn",
    "BlastHit",
    "Junction",
    "BreseqJunction",
    "RefTnJunction",
    "TnJunction",
    "RawTnJc2",
    "RefTnSide",
    "OffsetRefTnSide",
    "Genome",
    "GenomeRegistry",
    "get_genome",
    # Scaffold types
    "JcArm",
    # Coverage types
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

"""Data types and structures."""

from amplifinder.data_types.typed_df import RecordTypedDf, TypedDF
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.records import Column, Schema, Record
from amplifinder.data_types.enums import Terminal, Side, Orientation, AverageMethod, JunctionType, \
    EventModifier, JunctionReadCounts, BaseRawEvent
from amplifinder.data_types.scaffold import JcArm, Scaffold, SeqScaffold
from amplifinder.data_types.record_types import RefTn, BlastHit, Junction, BreseqJunction, \
    RefTnJunction, TnJunction, RawTnJc2, RefTnSide, OffsetRefTnSide, \
    CoveredTnJc2, RawEvent, SingleLocusLinkedTnJc2, SynJctsTnJc2, \
    AnalyzedTnJc2, ClassifiedTnJc2, ExportedTnJc2

__all__ = [
    "Terminal",
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
    "Scaffold",
    "SeqScaffold",

    # Coverage types
    "JunctionReadCounts",
    "CoveredTnJc2",
    "BaseRawEvent",
    "RawEvent",
    "SingleLocusLinkedTnJc2",
    "SynJctsTnJc2",
    "JunctionType",
    "EventModifier",
    "AnalyzedTnJc2",
    "ClassifiedTnJc2",
    "ExportedTnJc2",
]

"""Data types and structures."""

from amplifinder.records.typed_df import RecordTypedDf, TypedDF
from amplifinder.records.base_records import Column, Schema, Record
from amplifinder.data_types.genome import Genome, GenomeRegistry, get_genome
from amplifinder.data_types.basic_enums import Terminal, Side, Orientation, AverageMethod
from amplifinder.data_types.events import BaseEvent, Architecture, EventDescriptor
from amplifinder.data_types.read_types import ReadType, JunctionReadCounts
from amplifinder.data_types.jc_types import Element, JunctionType, JcCall
from amplifinder.data_types.scaffold import Scaffold, SeqScaffold
from amplifinder.data_types.junctions import JcArm, Junction, NumJunction, BreseqJunction
from amplifinder.data_types.ref_tn import TnId, RefTn, RefTnSide, OffsetRefTnSide, RefTnJunction, TnJunction
from amplifinder.data_types.tnjc2s import RawTnJc2, CoveredTnJc2, SingleLocusLinkedTnJc2, SynJctsTnJc2, \
    AnalyzedTnJc2, ClassifiedTnJc2, ExportedTnJc2

__all__ = [
    # Basic enums
    "Terminal",
    "Side",
    "Orientation",
    "AverageMethod",

    # Records base
    "Column",
    "Schema",
    "Record",
    "RecordTypedDf",
    "TypedDF",

    # Genome
    "Genome",
    "GenomeRegistry",
    "get_genome",

    # Scaffold types
    "JcArm",
    "Scaffold",
    "SeqScaffold",

    # Events
    "BaseEvent",
    "Architecture",
    "EventDescriptor",

    # Read types
    "ReadType",
    "JunctionReadCounts",

    # Junction types
    "Element",
    "JunctionType",
    "JcCall",

    # Reference TN
    "TnId",
    "RefTn",
    "RefTnSide",
    "OffsetRefTnSide",

    # Junctions
    "Junction",
    "NumJunction",
    "BreseqJunction",
    "RefTnJunction",
    "TnJunction",

    # TN junction pairs
    "RawTnJc2",
    "CoveredTnJc2",
    "SingleLocusLinkedTnJc2",
    "SynJctsTnJc2",
    "AnalyzedTnJc2",
    "ClassifiedTnJc2",
    "ExportedTnJc2",
]

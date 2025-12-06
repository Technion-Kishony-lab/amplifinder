"""Record dataclasses and their schemas."""
from __future__ import annotations

from dataclasses import dataclass, field, fields as dataclass_fields
from typing import Any, Dict, List, NamedTuple, TypeVar

from amplifinder.data_types.schema import TypedSchema

T = TypeVar("T", bound="Record")


@dataclass
class Record:
    """Base class for record dataclasses with from_instance support."""

    @classmethod
    def from_instance(cls: type[T], obj: Record, **kwargs) -> T:
        """Create instance from another record, copying shared fields."""
        shared = {f.name for f in dataclass_fields(cls)} & {f.name for f in dataclass_fields(obj)}
        data = {name: getattr(obj, name) for name in shared if name != "extra"}
        data.update(kwargs)
        return cls(**data)


class TnMatch(NamedTuple):
    """A single TN element match for a junction."""
    tn_id: int
    side: str      # "left" or "right"
    distance: int


@dataclass
class TnLoc(Record):
    """TN element location record."""
    ID: int
    TN_Name: str
    TN_scaf: str
    LocLeft: int
    LocRight: int
    Complement: bool
    Join: bool
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class BlastHit(Record):
    """BLAST alignment hit record."""
    query: str
    subject: str
    percent_identical: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Junction(Record):
    """Base junction record with shared positional fields."""
    num: int
    scaf1: str
    pos1: int
    dir1: int
    scaf2: str
    pos2: int
    dir2: int
    flanking_left: int
    flanking_right: int
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class RefTnJunction(Junction):
    """Synthetic junction for reference TN element."""
    refTN: int
    tn_side: str  # "left" or "right"
    extra: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Tnjc(Junction):
    """Junction matched to TN element(s)."""
    matches: List[TnMatch]  # TN matches: [(tn_id, side, distance), ...]
    switched: bool          # True if sides were swapped to normalize
    extra: Dict[str, Any] = field(default_factory=dict)


TN_LOC_SCHEMA = TypedSchema.from_dataclass(TnLoc)
BLAST_SCHEMA = TypedSchema.from_dataclass(BlastHit)
JUNCTION_SCHEMA = TypedSchema.from_dataclass(Junction)
REF_TN_JC_SCHEMA = TypedSchema.from_dataclass(RefTnJunction)
TNJC_SCHEMA = TypedSchema.from_dataclass(Tnjc)

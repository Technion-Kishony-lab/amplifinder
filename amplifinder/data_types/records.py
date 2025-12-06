"""Record dataclasses and their schemas."""

from dataclasses import dataclass, field
from typing import Any, Dict

from amplifinder.data_types.schema import TypedSchema


@dataclass
class TnLoc:
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
class BlastHit:
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


TN_LOC_SCHEMA = TypedSchema.from_dataclass(TnLoc)
BLAST_SCHEMA = TypedSchema.from_dataclass(BlastHit)

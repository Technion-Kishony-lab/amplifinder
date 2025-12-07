"""Record type definitions for AmpliFinder."""
from __future__ import annotations

from dataclasses import dataclass, replace
from typing import List, NamedTuple, TypeVar

from amplifinder.data_types.records import Record


class TnMatch(NamedTuple):
    """A single TN element match for a junction."""
    tn_id: int
    side: str      # "left" or "right"
    distance: int


@dataclass(kw_only=True)
class TnLoc(Record):
    """TN element location record."""
    ID: int
    TN_Name: str
    TN_scaf: str
    LocLeft: int
    LocRight: int
    Complement: bool
    Join: bool


@dataclass(kw_only=True)
class TnEndSeq(Record):
    """TN element end sequence for matching."""
    tn_id: int
    tn_side: str    # "left" or "right"
    seq_fwd: str    # forward sequence
    seq_rc: str     # reverse complement


@dataclass(kw_only=True)
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


JunctionT = TypeVar("JunctionT", bound="Junction")


@dataclass(kw_only=True)
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

    def switch_sides(self: JunctionT) -> JunctionT:
        """Return new junction with side 1 and side 2 swapped."""
        return replace(
            self,
            scaf1=self.scaf2, scaf2=self.scaf1,
            pos1=self.pos2, pos2=self.pos1,
            dir1=self.dir2, dir2=self.dir1,
            flanking_left=self.flanking_right, flanking_right=self.flanking_left,
        )


@dataclass(kw_only=True)
class RefTnJunction(Junction):
    """Synthetic junction for reference TN element."""
    refTN: int
    tn_side: str  # "left" or "right"


@dataclass(kw_only=True)
class TnJunction(Junction):
    """Junction matched to TN element(s)."""
    matches: List[TnMatch]  # TN matches: [(tn_id, side, distance), ...]
    switched: bool          # True if sides were swapped to normalize


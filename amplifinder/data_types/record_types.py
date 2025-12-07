"""Record type definitions for AmpliFinder."""
from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from typing import List, NamedTuple, Optional, TypeVar

from amplifinder.data_types.records import Record


class Side(str, Enum):
    """Side of a TN element (left or right)."""
    LEFT = "left"
    RIGHT = "right"

    def opposite(self) -> "Side":
        return Side.RIGHT if self == Side.LEFT else Side.LEFT


class Orientation(int, Enum):
    """Orientation relative to reference (forward, reverse, or both/mixed)."""
    FORWARD = 1
    REVERSE = -1
    BOTH = 0

    def opposite(self) -> "Orientation":
        if self == Orientation.FORWARD:
            return Orientation.REVERSE
        elif self == Orientation.REVERSE:
            return Orientation.FORWARD
        return Orientation.BOTH  # BOTH stays BOTH


class TnMatch(NamedTuple):
    """A single TN element match for a junction."""
    tn_id: int
    side: Side
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
    tn_side: Side
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
    dir1: Orientation
    scaf2: str
    pos2: int
    dir2: Orientation
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
    tn_side: Side


@dataclass(kw_only=True)
class TnJunction(Junction):
    """Junction matched to TN element(s)."""
    matches: List[TnMatch]  # TN matches: [(tn_id, side, distance), ...]
    switched: bool          # True if sides were swapped to normalize


@dataclass(kw_only=True)
class TnJunctionPair(Record):
    """Paired TN junctions (candidate amplicon).
    
    Represents two junctions that:
    - Are on the same scaffold facing opposite directions
    - Match the same TN element on different sides (left/right)
    
    Based on MATLAB combine_ISJC_pairs.m
    """
    # Junction IDs
    jc_num_L: int
    jc_num_R: int
    
    # Scaffold
    scaf_chr: str
    
    # Chromosome positions (left/right junction)
    pos_chr_L: int
    pos_chr_R: int
    
    # TN positions
    pos_tn_L: int
    pos_tn_R: int
    
    # Chromosome directions
    dir_chr_L: Orientation
    dir_chr_R: Orientation
    
    # TN directions
    dir_tn_L: Orientation
    dir_tn_R: Orientation
    
    # TN info
    tn_ids: List[int]        # matching TN element IDs
    tn_orientation: Orientation
    span_origin: bool        # True if amplicon spans circular origin
    
    # Computed fields
    amplicon_length: int
    complementary_length: int
    
    # Coverage (optional, computed later)
    amplicon_coverage: Optional[float] = None
    copy_number_mode: Optional[float] = None


"""Record type definitions for AmpliFinder."""
from __future__ import annotations

from enum import Enum
from operator import add, mul, sub, truediv
from typing import List, NamedTuple, Optional, TypeVar

from amplifinder.data_types.records import Record


TnId = int

NAMES_TO_OPERATORS = {
    "__truediv__": truediv,
    "__mul__": mul,
    "__add__": add,
    "__sub__": sub,
}


# Average types
class Average(NamedTuple):
    """Coverage statistics for a genomic region."""
    mean: float
    median: float
    mode: float


# add operator methods
for op_name, op in NAMES_TO_OPERATORS.items():
    def wrapper(self, other: Average) -> Average:
        return Average(op(self.mean, other.mean), op(self.median, other.median), op(self.mode, other.mode))
    setattr(Average, op_name, wrapper)


class JunctionCoverage(NamedTuple):
    """Read coverage at a synthetic junction."""
    spanning: int  # reads crossing junction
    left: int      # reads ending at junction
    right: int     # reads starting at junction


class ReversibleIntEnum(int, Enum):
    """Base class for int enums with opposite() method."""

    def opposite(self):
        """Return the enum member with negated value."""
        return type(self)(-self.value)


class Side(ReversibleIntEnum):
    """Side of a TN element (left or right)."""
    LEFT = -1
    RIGHT = 1


class Orientation(ReversibleIntEnum):
    """Orientation relative to reference (forward, reverse, or both/mixed)."""
    FORWARD = 1
    REVERSE = -1
    BOTH = 0


### Reference TN element sides ###

class RefTnSide(Record):
    """A reference TN element side (with optional distance for matches)."""
    tn_id: TnId
    side: Side
    distance: Optional[int] = None  # None for reference junctions


class SeqRefTnSide(RefTnSide):
    """TN element end sequence for matching."""
    offset: int         # offset of the inward-seq start from the TN boundary (>0 for inward, <0 for outward)
    seq_inward: str     # sequence inward from chromosome into TN


### Reference TN elements ###

class RefTnLoc(Record):
    """Reference TN element location in the genome."""
    tn_id: TnId
    tn_name: str
    tn_scaf: str
    loc_left: int
    loc_right: int
    orientation: Orientation
    join: bool

    @property
    def length(self) -> int:
        """TN length in bp (1-based inclusive coordinates)."""
        return self.loc_right - self.loc_left + 1

    def get_sides(self) -> tuple[RefTnSide, RefTnSide]:
        """Get left and right sides of the TN."""
        return (
            RefTnSide(tn_id=self.tn_id, side=Side.LEFT),
            RefTnSide(tn_id=self.tn_id, side=Side.RIGHT),
        )

    def get_junctions(self, out_span: int, in_span: Optional[int] = None) -> tuple[RefTnJunction, RefTnJunction]:
        """Get left and right junctions of the TN.
        
        Junction numbering: negative values for reference junctions (breseq uses positive).
        Left: -tn_id * 2 (even negative), Right: -tn_id * 2 - 1 (odd negative).
        """
        if in_span is None:
            in_span = self.length
        return (
            RefTnJunction(
                num=-self.tn_id * 2,  # Left: even negative
                scaf1=self.tn_scaf, pos1=self.loc_left, dir1=Orientation.FORWARD,
                scaf2=self.tn_scaf, pos2=self.loc_left - 1, dir2=Orientation.REVERSE,
                flanking_left=in_span, flanking_right=out_span,
                ref_tn_side=RefTnSide(tn_id=self.tn_id, side=Side.LEFT),
            ),
            RefTnJunction(
                num=-self.tn_id * 2 - 1,  # Right: odd negative
                scaf1=self.tn_scaf, pos1=self.loc_right, dir1=Orientation.REVERSE,
                scaf2=self.tn_scaf, pos2=self.loc_right + 1, dir2=Orientation.FORWARD,
                flanking_left=in_span, flanking_right=out_span,
                ref_tn_side=RefTnSide(tn_id=self.tn_id, side=Side.RIGHT),
            ),
        )


### BLAST hits ###

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


### Junctions ###

JunctionT = TypeVar("JunctionT", bound="Junction")


class Junction(Record):
    """Base junction record with shared positional fields."""
    num: int  # Junction identifier: breseq junction number (positive), or negative for reference junctions
    scaf1: str
    pos1: int
    dir1: Orientation
    scaf2: str
    pos2: int
    dir2: Orientation
    flanking_left: int   # Length of sequence flanking side 1 (used for sequence extraction)
    flanking_right: int  # Length of sequence flanking side 2 (used for sequence extraction)

    def switch_sides(self: JunctionT) -> JunctionT:
        """Return new junction with side 1 and side 2 swapped."""
        return self.model_copy(update={
            "scaf1": self.scaf2, "scaf2": self.scaf1,
            "pos1": self.pos2, "pos2": self.pos1,
            "dir1": self.dir2, "dir2": self.dir1,
            "flanking_left": self.flanking_right, "flanking_right": self.flanking_left,
        })

    def get_scaf_pos_dir_flank(self, side: int) -> tuple[str, int, Orientation, int]:
        """Get scaffold, position, direction, and flanking length for a side."""
        return (
            self.scaf1 if side == 1 else self.scaf2, 
            self.pos1 if side == 1 else self.pos2, 
            self.dir1 if side == 1 else self.dir2, 
            self.flanking_left if side == 1 else self.flanking_right
            )


class RefTnJunction(Junction):
    """Synthetic junction for reference TN element.

    For RefTnJunction, side 1 is always the TN side, side 2 is the chromosome side.
    ref_tn_side indicates which TN boundary (LEFT or RIGHT) this junction represents.
    """
    ref_tn_side: RefTnSide


class TnJunction(Junction):
    """Junction matched to TN element(s)."""
    ref_tn_sides: List[RefTnSide]  # Reference TN matches: [(tn_id, side, distance?), ...]
    switched: bool                 # True if BRESEQ sides were swapped (to normalize to TN on side 1)


### Paired TN junctions ###

class RawTnJc2(Record):
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
    scaf: str

    # Scafold positions (left/right junction)
    pos_scaf_L: int
    pos_scaf_R: int

    # TN positions
    pos_tn_L: int
    pos_tn_R: int

    # Scaffold directions
    dir_scaf_L: Orientation
    dir_scaf_R: Orientation

    # TN directions
    dir_tn_L: Orientation
    dir_tn_R: Orientation

    # TN info
    tn_ids: List[int]              # matching TN element IDs
    tn_orientations: List[Orientation]  # one per tn_id
    span_origin: bool        # True if amplicon spans circular origin

    # Computed fields
    amplicon_length: int
    complementary_length: int


class CoveredTnJc2(RawTnJc2):
    """RawTnJc2 with coverage information (Step 7 output).

    Coverage fields depend on run type:
    - anc_path=None: raw coverage only, copy_number_ratio is None
    - anc_path=set: normalized coverage, copy_number_ratio = iso/anc
    """
    iso_amplicon_coverage: Average
    iso_scaf_coverage: Average
    anc_amplicon_coverage: Optional[Average] = None
    anc_scaf_coverage: Optional[Average] = None
    copy_number: float = None  #
    copy_number_vs_anc: Optional[float] = None


class RawEvent(str, Enum):
    """Structural classification based on junction pair relationships (Step 8)."""
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"
    UNFLANKED = "unflanked"
    HEMI_FLANKED_LEFT = "hemi-flanked left"
    HEMI_FLANKED_RIGHT = "hemi-flanked right"
    FLANKED = "flanked"
    MULTIPLE_SINGLE_LOCUS = "multiple single locus"
    UNRESOLVED = "unresolved"


class ClassifiedTnJc2(CoveredTnJc2):
    """CoveredTnJc2 with structural classification (Step 8 output)."""
    raw_event: RawEvent
    shared_tn_ids: List[int]        # TN IDs shared by both junctions
    chosen_tn_id: Optional[int]     # selected TN for analysis


class FilteredTnJc2(ClassifiedTnJc2):
    """Filtered candidate for detailed analysis (Step 9 output)."""
    analysis_dir: str  # "tn_jc2_001", "tn_jc2_002", etc.


class JunctionType(int, Enum):
    """The 7 synthetic junction types.

    Amplicon structure: ~~~>>>======>>>======>>>~~~
    (1) ~~==  left reference (chromosome-cassette)
    (2) ~~>>  left IS transposition (chromosome-IS)
    (3) ==>>  left of mid IS (cassette-IS)
    (4) ====  lost IS (cassette-cassette, no IS)
    (5) >>==  right of mid IS (IS-cassette)
    (6) >>~~  right IS transposition (IS-chromosome)
    (7) ==~~  right reference (cassette-chromosome)
    """
    LEFT_REF = 1
    LEFT_IS_TRANS = 2
    LEFT_MID_IS = 3
    LOST_IS = 4
    RIGHT_MID_IS = 5
    RIGHT_IS_TRANS = 6
    RIGHT_REF = 7


class EventModifier(str, Enum):
    """Modifiers for classified events (Step 13)."""
    ANCESTRAL = "ancestral"
    DE_NOVO = "de novo"
    LOW_COVERAGE = "low coverage near junction"


class AnalyzedTnJc2(FilteredTnJc2):
    """Candidate with junction coverage analysis (Step 12 output).

    Junction coverage fields depend on run type:
    - anc_path=None: jc_cov only, anc_jc_cov is None
    - anc_path=set: both jc_cov and anc_jc_cov present
    """
    # Junction coverage (7 elements, one per JunctionType)
    jc_cov_left: List[int]      # left-side read counts per junction
    jc_cov_right: List[int]     # right-side read counts per junction
    jc_cov_spanning: List[int]  # spanning read counts per junction

    # Ancestor junction coverage (only when anc_path is set)
    anc_jc_cov_left: Optional[List[int]] = None
    anc_jc_cov_right: Optional[List[int]] = None
    anc_jc_cov_spanning: Optional[List[int]] = None

    # Architecture classification
    isolate_architecture: RawEvent
    ancestor_architecture: Optional[RawEvent] = None  # only when anc_path is set

    # Final event classification
    event: str                              # full event description
    event_modifiers: List[EventModifier]    # de novo, ancestral, etc.


class ExportedTnJc2(Record):
    """Export record for tnjc2_exported.csv (Step 14 output).

    Represents the user-facing export format with renamed/combined fields.
    All fields are optional to handle cases where input data may be missing.
    """
    isolate: Optional[str] = None
    Reference: Optional[str] = None
    Positions_in_chromosome: Optional[str] = None
    Direction_in_chromosome: Optional[str] = None
    amplicon_length: Optional[int] = None
    IS_element: Optional[str] = None
    median_copy_number: Optional[float] = None
    mode_copy_number: Optional[float] = None
    Ancestor: Optional[str] = None
    event: Optional[str] = None
    isolate_architecture: Optional[str] = None

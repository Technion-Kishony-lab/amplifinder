"""Record type definitions for AmpliFinder."""
from __future__ import annotations

from enum import Enum
from typing import List, NamedTuple, Optional, TypeVar

from amplifinder.data_types.records import Record


# Coverage types
class Coverage(NamedTuple):
    """Coverage statistics for a genomic region."""
    mean: float
    median: float
    mode: float


class JunctionCoverage(NamedTuple):
    """Read coverage at a synthetic junction."""
    spanning: int  # reads crossing junction
    left: int      # reads ending at junction
    right: int     # reads starting at junction


class Side(int, Enum):
    """Side of a TN element (left or right). Values match MATLAB convention."""
    LEFT = -1
    RIGHT = 1

    def opposite(self) -> "Side":
        return Side(-self.value)


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


class RefTnSide(Record):
    """A reference TN element side (with optional distance for matches)."""
    tn_id: int
    side: Side
    distance: Optional[int] = None


class TnLoc(Record):
    """TN element location record."""
    ID: int
    TN_Name: str
    TN_scaf: str
    LocLeft: int
    LocRight: int
    Complement: bool
    Join: bool


class SeqRefTnSide(RefTnSide):
    """TN element end sequence for matching."""
    seq_fwd: str    # forward sequence
    seq_rc: str     # reverse complement


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
        return self.model_copy(update={
            "scaf1": self.scaf2, "scaf2": self.scaf1,
            "pos1": self.pos2, "pos2": self.pos1,
            "dir1": self.dir2, "dir2": self.dir1,
            "flanking_left": self.flanking_right, "flanking_right": self.flanking_left,
        })


class RefTnJunction(Junction):
    """Synthetic junction for reference TN element."""
    ref_tn_side: RefTnSide


class TnJunction(Junction):
    """Junction matched to TN element(s)."""
    ref_tn_sides: List[RefTnSide]  # Reference TN matches: [(tn_id, side, distance?), ...]
    switched: bool          # True if sides were swapped to normalize


class TnJc2(Record):
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
    tn_ids: List[int]              # matching TN element IDs
    tn_orientations: List[Orientation]  # one per tn_id
    span_origin: bool        # True if amplicon spans circular origin

    # Computed fields
    amplicon_length: int
    complementary_length: int


class CoveredTnJc2(TnJc2):
    """TnJc2 with coverage information (Step 7 output).
    
    Coverage fields depend on run type:
    - anc_path=None: raw coverage only, copy_number_ratio is None
    - anc_path=set: normalized coverage, copy_number_ratio = iso/anc
    """
    # Context
    ref_name: str
    iso_name: str
    anc_name: Optional[str] = None  # present when anc_path was set
    
    # Coverage fields (always present)
    amplicon_coverage: float        # raw: iso_cov/genome_cov, normalized: iso_cov/anc_cov
    genome_coverage: float          # median genome coverage for this sample
    copy_number: float              # amplicon / genome (raw copy number)
    amplicon_coverage_mode: float   # mode of copy number distribution
    
    # Ancestor comparison (only when anc_path is set)
    copy_number_ratio: Optional[float] = None  # iso_copy / anc_copy


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
    """TnJc2 with structural classification (Step 8 output)."""
    raw_event: RawEvent
    shared_tn_ids: List[int]        # TN IDs shared by both junctions
    chosen_tn_id: Optional[int]     # selected TN for analysis


class CandidateTnJc2(ClassifiedTnJc2):
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


class AnalyzedTnJc2(CandidateTnJc2):
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

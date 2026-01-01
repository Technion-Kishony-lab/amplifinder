"""Record type definitions for AmpliFinder."""
from __future__ import annotations

from typing import ClassVar, List, NamedTuple, Optional, TypeVar, TYPE_CHECKING
from enum import Enum

from amplifinder.data_types.records import Record
from amplifinder.data_types.enums import Side, Orientation
from amplifinder.data_types.scaffold import SegmentScaffold, JcArm

if TYPE_CHECKING:
    from amplifinder.data_types.genome import Genome

TnId = int


class JunctionCoverage(NamedTuple):
    """Read coverage at a synthetic junction."""
    spanning: int  # reads crossing junction
    left: int      # reads ending at junction
    right: int     # reads starting at junction


# ===== Reference TN element sides =====

class RefTnSide(Record):
    """A reference TN element side."""
    NAME: ClassVar[str] = "Reference TN sides"
    tn_id: TnId
    side: Side

    def is_same_side(self, other: RefTnSide) -> bool:
        """Check if two RefTnSide objects are the same side."""
        return self.tn_id == other.tn_id and self.side == other.side


class OffsetRefTnSide(RefTnSide):
    """A reference TN element side with offset distance for matches."""
    NAME: ClassVar[str] = "Offset Reference TN sides"
    distance: int


# ===== Reference TN elements =====

class RefTnLoc(Record):
    """Reference TN element location in the genome."""
    NAME: ClassVar[str] = "TN elements"
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

    def get_junctions(
            self, out_span: int, in_span: Optional[int] = None
            ) -> tuple[RefTnJunction, RefTnJunction]:
        """Get left and right junctions of the TN.

        Junction numbering: negative values for reference junctions.
        Left: -tn_id * 2, Right: -tn_id * 2 - 1 (odd negative).
        """
        if in_span is None:
            in_span = self.length
        side_left, side_right = self.get_sides()
        return (
            RefTnJunction(
                num=-self.tn_id * 2,  # Left: even negative
                scaf1=self.tn_scaf, pos1=self.loc_left, dir1=Orientation.FORWARD,
                scaf2=self.tn_scaf, pos2=self.loc_left - 1, dir2=Orientation.REVERSE,
                flanking_left=in_span, flanking_right=out_span,
                ref_tn_side=side_left,
            ),
            RefTnJunction(
                num=-self.tn_id * 2 - 1,  # Right: odd negative
                scaf1=self.tn_scaf, pos1=self.loc_right, dir1=Orientation.REVERSE,
                scaf2=self.tn_scaf, pos2=self.loc_right + 1, dir2=Orientation.FORWARD,
                flanking_left=in_span, flanking_right=out_span,
                ref_tn_side=side_right,
            ),
        )


# ===== BLAST hits =====

class BlastHit(Record):
    """BLAST alignment hit record."""
    NAME: ClassVar[str] = "BLAST hits"
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


# ===== Junctions =====

JunctionT = TypeVar("JunctionT", bound="Junction")


class Junction(Record):
    """Base junction record with shared positional fields."""
    NAME: ClassVar[str] = "Junctions"
    num: int  # Junction identifier: breseq junction number (positive), or negative for reference junctions
    scaf1: str
    pos1: int
    dir1: Orientation
    scaf2: str
    pos2: int
    dir2: Orientation
    flanking_left: int   # Length of sequence flanking arm 1 (used for sequence extraction)
    flanking_right: int  # Length of sequence flanking arm 2 (used for sequence extraction)

    def swap_sides(self: JunctionT) -> JunctionT:
        """Return new junction with arm 1 and arm 2 swapped."""
        return self.model_copy(update={
            "scaf1": self.scaf2, "scaf2": self.scaf1,
            "pos1": self.pos2, "pos2": self.pos1,
            "dir1": self.dir2, "dir2": self.dir1,
            "flanking_left": self.flanking_right, "flanking_right": self.flanking_left,
        })

    def get_jc_arm(self, arm: int) -> JcArm:
        """Get scaffold, position, direction, and flanking length for an arm."""
        return JcArm(
            scaf=self.scaf1 if arm == 1 else self.scaf2,
            start=self.pos1 if arm == 1 else self.pos2,
            dir=self.dir1 if arm == 1 else self.dir2,
            flank=self.flanking_left if arm == 1 else self.flanking_right
        )


class BreseqJunction(Junction):
    """Breseq junction."""
    NAME: ClassVar[str] = "Breseq junctions"


class RefTnJunction(Junction):
    """Synthetic junction for reference TN element.

    For RefTnJunction, arm 1 is always the TN side, arm 2 is the chromosome side.
    ref_tn_side indicates which TN boundary (LEFT or RIGHT) this junction represents.
    """
    NAME: ClassVar[str] = "Reference TN junctions"

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Side.LEFT
    #        |------>              arm1, flanking_left (into TN)
    #     <--|                     arm2, flanking_right (out of TN)

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Side.RIGHT
    #           <------|           arm1, flanking_left (into TN)
    #                  |-->        arm2, flanking_right (out of TN)

    ref_tn_side: RefTnSide


class TnJunction(Junction):
    """Junction matched to TN element(s)."""
    NAME: ClassVar[str] = "TN-associated junctions"
    ref_tn_side: Optional[RefTnSide] = None  # None for breseq junctions
    ref_tn_sides: List[OffsetRefTnSide]      # Reference TN matches with distances
    swapped: bool                            # True if BRESEQ arms were swapped (to normalize to TN on arm 1)

    def is_ref_tn_junction(self) -> bool:
        """Return True if this is a reference TN junction."""
        return self.ref_tn_side is not None


# Paired TN junctions

class RawTnJc2(Record):
    """Paired TN junctions (candidate amplicon).

    Represents two junctions that:
    - Are on the same scaffold facing opposite directions
    - Match the same TN element on different sides (left/right)

    # nomenclature:
    # start (S): the junction from which we start the amplicon segment going on the forward strand
    # end (E):   the junction at which we end the amplicon segment when we go on the forward strand

    modified from combine_ISJC_pairs.m

    """
    NAME: ClassVar[str] = "Junction Pairs"
    # Junction IDs
    jc_num_S: int  # Junction number of the 'start' junction
    jc_num_E: int  # Junction number of the 'end' junction

    # Scaffold
    scaf: str

    # Scaffold positions (start/end of amplicon segment on the forward strand)
    start: int
    end: int

    # TN positions
    pos_tn_S: int  # TN position of the 'start' junction
    pos_tn_E: int  # TN position of the 'end' junction

    # TN directions
    dir_tn_S: Orientation  # TN direction of the 'start' junction
    dir_tn_E: Orientation  # TN direction of the 'end' junction

    # TN info
    tn_ids: List[int]              # matching TN element IDs
    tn_orientations: List[Orientation]  # one per tn_id

    # Computed fields
    amplicon_length: int = None  # None until computed

    def compute_and_store_amplicon_length(self, genome: Genome):
        """Compute amplicon length using the genome."""
        self.amplicon_length = self.get_segment_scaffold(genome).segment_length

    def get_segment_scaffold(self, genome) -> SegmentScaffold:
        """Get SegmentScaffold for this amplicon segment.

        Returns a SegmentScaffold with start/end positions that provides
        properties: left, right, span_origin, segment_length.
        """
        scaf_obj = genome.get_scaffold(self.scaf)
        return SegmentScaffold(
            is_circular=scaf_obj.is_circular,
            length=scaf_obj.length,
            start=self.start,
            end=self.end,
            orientation=Orientation.FORWARD,  # by definition, the amplicon segment is on the forward strand
        )


class CoveredTnJc2(RawTnJc2):
    """RawTnJc2 with coverage information (Step 7 output).

    Coverage fields depend on run type:
    - anc_path=None: scaffold-normalized only
    - anc_path=set: scaffold-normalized then ancestor-normalized
    """
    NAME: ClassVar[str] = "Covered Junction Pairs"
    iso_scaf_avg: float
    iso_amplicon_avg: float
    anc_scaf_avg: Optional[float] = None
    anc_amplicon_avg: Optional[float] = None
    avg_norm_cov: Optional[float] = None  # Position-by-position ancestor-normalized average

    @property
    def iso_scaf_norm_copy_number(self) -> float:
        """Isolate copy number (scaffold-normalized)."""
        return self.iso_amplicon_avg / self.iso_scaf_avg

    @property
    def anc_scaf_norm_copy_number(self) -> Optional[float]:
        """Ancestor copy number (scaffold-normalized)."""
        if self.anc_amplicon_avg is None or self.anc_scaf_avg is None:
            return None
        return self.anc_amplicon_avg / self.anc_scaf_avg

    @property
    def scaf_norm_copy_number_ratio(self) -> Optional[float]:
        """Ratio of scaffold-normalized copy numbers (iso/anc)."""
        anc_cn = self.anc_scaf_norm_copy_number
        if anc_cn is None:
            return None
        return self.iso_scaf_norm_copy_number / anc_cn

    @property
    def copy_number(self) -> float:
        """Final copy number: ancestor-normalized if available, else scaffold-normalized."""
        return self.avg_norm_cov if self.avg_norm_cov is not None else self.iso_scaf_norm_copy_number


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
    NAME: ClassVar[str] = "Classified Amplicons"
    raw_event: RawEvent
    shared_tn_ids: List[int]        # TN IDs shared by both junctions
    chosen_tn_id: Optional[int]     # selected TN for analysis


class FilteredTnJc2(ClassifiedTnJc2):
    """Filtered candidate for detailed analysis (Step 9 output)."""
    NAME: ClassVar[str] = "Filtered Amplicons"
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
    NAME: ClassVar[str] = "Analyzed Amplicons"
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
    NAME: ClassVar[str] = "Exported Amplicons"
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

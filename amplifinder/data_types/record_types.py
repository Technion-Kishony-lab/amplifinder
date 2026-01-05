"""Record type definitions for AmpliFinder."""
from __future__ import annotations

from typing import ClassVar, List, NamedTuple, Optional, TypeVar, TYPE_CHECKING
from enum import Enum

from amplifinder.data_types.records import Record
from amplifinder.data_types.enums import BaseRawEvent, RawEvent, Side, Orientation, EventModifier
from amplifinder.data_types.scaffold import SegmentScaffold, JcArm

if TYPE_CHECKING:
    from amplifinder.data_types.genome import Genome
    from amplifinder.data_types.scaffold import Scaffold, SeqScaffold

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

class RefTn(SegmentScaffold):
    """Reference TN element location in the genome.
    
    Inherits from SegmentScaffold for coordinate operations.
    """
    NAME: ClassVar[str] = "TN elements"
    tn_id: TnId
    tn_name: str
    join: bool

    def get_ref_tn_sides(self) -> tuple[RefTnSide, RefTnSide]:
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
            in_span = self.segment_length
        side_left, side_right = self.get_ref_tn_sides()
        return (
            RefTnJunction(
                num=-self.tn_id * 2,  # Left: even negative
                scaf1=self.scaf, pos1=self.start, dir1=Orientation.FORWARD,
                scaf2=self.scaf, pos2=self.start - 1, dir2=Orientation.REVERSE,
                flanking_left=in_span, flanking_right=out_span,
                ref_tn_side=side_left,
            ),
            RefTnJunction(
                num=-self.tn_id * 2 - 1,  # Right: odd negative
                scaf1=self.scaf, pos1=self.end, dir1=Orientation.REVERSE,
                scaf2=self.scaf, pos2=self.end + 1, dir2=Orientation.FORWARD,
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
    num: int = None  # Junction identifier: breseq junction number (positive), or negative for reference junctions
    scaf1: str
    pos1: int
    dir1: Orientation
    scaf2: str
    pos2: int
    dir2: Orientation
    flanking1: int   # Length of sequence flanking arm 1 (used for sequence extraction)
    flanking2: int  # Length of sequence flanking arm 2 (used for sequence extraction)

    def swap_sides(self: JunctionT) -> JunctionT:
        """Return new junction with arm 1 and arm 2 swapped."""
        return self.model_copy(update={
            "scaf1": self.scaf2, "scaf2": self.scaf1,
            "pos1": self.pos2, "pos2": self.pos1,
            "dir1": self.dir2, "dir2": self.dir1,
            "flanking1": self.flanking2, "flanking2": self.flanking1,
        })

    def get_jc_arm(self, arm: int) -> JcArm:
        """Get scaffold, position, direction, and flanking length for an arm."""
        return JcArm(
            scaf=self.scaf1 if arm == 1 else self.scaf2,
            start=self.pos1 if arm == 1 else self.pos2,
            dir=self.dir1 if arm == 1 else self.dir2,
            flank=self.flanking1 if arm == 1 else self.flanking2
        )

    def __eq__(self, other: Junction) -> bool:
        """Check if two junctions are the same."""
        return self.num == other.num

    @classmethod
    def from_jc_arms(cls, arm1: JcArm, arm2: JcArm) -> Junction:
        """Create a Junction from junction arm coordinates."""
        return cls(scaf1=arm1.scaf, pos1=arm1.start, dir1=arm1.dir, flanking1=arm1.flank, 
                   scaf2=arm2.scaf, pos2=arm2.start, dir2=arm2.dir, flanking2=arm2.flank)

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
    #        |------>              arm1, flanking1 (into TN)
    #     <--|                     arm2, flanking2 (out of TN)

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Side.RIGHT
    #           <------|           arm1, flanking1 (into TN)
    #                  |-->        arm2, flanking2 (out of TN)

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
    
    # Core fields: the two junctions
    tnjc_S: TnJunction  # Start junction (dir2 == FORWARD)
    tnjc_E: TnJunction  # End junction (dir2 == REVERSE)
    
    # Cached computed value (requires genome to compute)
    amplicon_length: Optional[int] = None
    
    # CSV export: only export derived properties, not the complex TnJunction objects
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = [
        'jc_num_S', 'jc_num_E', 'scaf', 'start', 'end',
        'pos_tn_S', 'pos_tn_E', 'dir_tn_S', 'dir_tn_E',
        'tn_ids', 'tn_orientations', 'tn_distances', 'amplicon_length'
    ]

    @staticmethod
    def find_matching_tn_sides(
        i_ref_tn_sides: List[OffsetRefTnSide],
        j_ref_tn_sides: List[OffsetRefTnSide],
    ) -> List[tuple[OffsetRefTnSide, OffsetRefTnSide]]:
        """Find TN elements that match both junctions on different sides.
        
        Returns list of OffsetRefTnSide objects from the first junction (i) that match
        the second junction (j) on different sides.
        """
        return [(i_tn_side, j_tn_side) for i_tn_side in i_ref_tn_sides for j_tn_side in j_ref_tn_sides
                if i_tn_side.tn_id == j_tn_side.tn_id    # same TN ID
                and i_tn_side.side != j_tn_side.side]    # different sides

    def _find_matching_tn_sides(self) -> List[tuple[OffsetRefTnSide, OffsetRefTnSide]]:
        """Instance method: find matching TNs using this record's junctions.
        
        Returns list of OffsetRefTnSide objects from tn_jc_S that match tn_jc_E.
        """
        return self.find_matching_tn_sides(self.tnjc_S.ref_tn_sides, self.tnjc_E.ref_tn_sides)

    @property
    def scaf(self) -> str:
        """Scaffold name."""
        assert self.tnjc_S.scaf2 == self.tnjc_E.scaf2
        return self.tnjc_S.scaf2

    @property
    def start(self) -> int:
        """Start position of amplicon segment on the forward strand."""
        return self.tnjc_S.pos2

    @property
    def end(self) -> int:
        """End position of amplicon segment on the forward strand."""
        return self.tnjc_E.pos2

    @property
    def tn_ids(self) -> List[int]:
        """Matching TN element IDs."""
        return [tn_side_S.tn_id for tn_side_S, tn_side_E in self._find_matching_tn_sides()]

    @property
    def tn_distances(self) -> List[tuple[int, int]]:
        """TN distances, one per tn_id."""
        return [(tn_side_S.distance, tn_side_E.distance) for tn_side_S, tn_side_E in self._find_matching_tn_sides()]

    def compute_and_store_amplicon_length(self, genome: Genome):
        """Compute amplicon length using the genome."""
        self.amplicon_length = self.get_segment_scaffold(genome).segment_length

    def get_segment_scaffold(self, genome) -> SegmentScaffold:
        """Get SegmentScaffold for this amplicon segment.

        Returns a SegmentScaffold with start/end positions that provides
        properties: left, right, span_origin, segment_length.
        """
        scaf_obj = genome.get_scaffold(self.scaf)
        return SegmentScaffold.from_other(
            scaf_obj,
            start=self.start,
            end=self.end,
            orientation=Orientation.FORWARD,  # by definition, the amplicon segment is on the forward strand
        )

    def __str__(self) -> str:
        """String representation of RawTnJc2."""
        return f"RawTnJc2({self.start}-{self.end}, scaf={self.scaf}, len={self.amplicon_length}, tn_ids={self.tn_ids})"


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


class ClassifiedTnJc2(CoveredTnJc2):
    """CoveredTnJc2 with structural classification (Step 8 output)."""
    NAME: ClassVar[str] = "Classified Amplicons"
    tnjc2_matching_S: Optional[CoveredTnJc2] = None
    tnjc2_matching_E: Optional[CoveredTnJc2] = None
    base_raw_event: BaseRawEvent

    @property
    def raw_event(self) -> RawEvent:
        """Raw event classification."""
        if self.base_raw_event == BaseRawEvent.REFERENCE:
            return RawEvent.REFERENCE
        elif self.base_raw_event == BaseRawEvent.TRANSPOSITION:
            return RawEvent.TRANSPOSITION
        assert self.base_raw_event == BaseRawEvent.LOCUS_JOINING
        has_match_S = self.tnjc2_matching_S is not None 
        has_match_E = self.tnjc2_matching_E is not None
        if has_match_S and has_match_E:
            return RawEvent.FLANKED
        elif has_match_S:
            return RawEvent.HEMI_FLANKED_LEFT
        elif has_match_E:
            return RawEvent.HEMI_FLANKED_RIGHT
        else:
            return RawEvent.UNFLANKED
    
    @property
    def chosen_tn_id(self) -> TnId:
        """Selected TN for analysis."""
        # Start with TN IDs from this tnjc2
        tn_id_set = set(self.tn_ids)
        
        # Intersect with matching S/E tnjc2 if exists
        if self.tnjc2_matching_S is not None:
            tn_id_set &= set(self.tnjc2_matching_S.tn_ids)
        if self.tnjc2_matching_E is not None:
            tn_id_set &= set(self.tnjc2_matching_E.tn_ids)
        
        # Return first if available, otherwise None
        return list(tn_id_set)[0] if tn_id_set else None

    def get_sides_of_chosen_tn(self) -> tuple[Optional[OffsetRefTnSide], Optional[OffsetRefTnSide]]:
        """Get sides of chosen TN."""
        chosen_id = self.chosen_tn_id
        if chosen_id is None:
            return None, None
        matching_tns = self._find_matching_tn_sides()
        for tn_side_S, tn_side_E in matching_tns:
            if tn_side_S.tn_id == chosen_id:
                return tn_side_S, tn_side_E
        assert False

    @property
    def chosen_tn_side_S(self) -> Optional[Side]:
        """Side of chosen TN that the Start junction connects to."""
        chosen_id = self.chosen_tn_id
        if chosen_id is None:
            return None
        matching_tns = self._find_matching_tn_sides()
        for tn_side_S, tn_side_E in matching_tns:
            assert tn_side_S.side == -tn_side_E.side
            assert tn_side_S.tn_id == tn_side_E.tn_id
            if tn_side_S.tn_id == chosen_id:
                return tn_side_S.side
        return None


class FilteredTnJc2(ClassifiedTnJc2):
    """Filtered candidate for detailed analysis (Step 9 output)."""
    NAME: ClassVar[str] = "Filtered Amplicons"
    analysis_dir: str  # "tn_jc2_001", "tn_jc2_002", etc.


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

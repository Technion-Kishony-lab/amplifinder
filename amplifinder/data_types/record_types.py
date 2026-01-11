"""Record type definitions for AmpliFinder."""
from __future__ import annotations
from typing import Any, ClassVar, Dict, List, Optional, TypeVar
from pathlib import Path
from pydantic import ConfigDict, field_validator
import numpy as np

from amplifinder.data_types.records import Record
from amplifinder.data_types.enums import BaseRawEvent, RawEvent, Side, Orientation, EventModifier, JunctionType, \
    JunctionReadCounts
from amplifinder.data_types.scaffold import SegmentScaffold, JcArm, SeqScaffold, SeqSegmentScaffold

TnId = int


# ===== Reference TN element sides =====

class RefTnSide(Record):
    """A reference TN element side."""
    NAME: ClassVar[str] = "Reference TN sides"
    tn_id: TnId
    side: Side

    def is_same_side(self, other: RefTnSide) -> bool:
        """Check if two RefTnSide objects are the same side."""
        return self.tn_id == other.tn_id and self.side == other.side

    def is_opposite_side(self, other: RefTnSide) -> bool:
        """Check if two RefTnSide objects are on opposite sides."""
        return self.tn_id == other.tn_id and self.side != other.side


class OffsetRefTnSide(RefTnSide):
    """A reference TN element side with offset for matches."""
    NAME: ClassVar[str] = "Offset Reference TN sides"
    offset: int  # >0 into TN (inward), <0 away from TN into the chromosome (outward)


# ===== Reference TN elements =====

class RefTn(SegmentScaffold):
    """Reference TN element location in the genome.
    """
    NAME: ClassVar[str] = "Reference TN elements"
    tn_id: TnId
    tn_name: str
    join: bool

    def get_ref_tn_sides(self) -> tuple[RefTnSide, RefTnSide]:
        """Get start and end sides of the TN."""
        return (
            RefTnSide(tn_id=self.tn_id, side=Side.START),
            RefTnSide(tn_id=self.tn_id, side=Side.END),
        )

    def get_junctions(self, out_flanks: int | tuple[int, int],
                      in_flanks: int | tuple[int, int] | None = None,
                      junction_class=None):
        """Get start and end reference TN junctions.

        Note:
            - Junction numbering: Start = -tn_id*2 (even), End = -tn_id*2-1 (odd)
            - Arm 1 is always the TN side, Arm 2 is always the chromosome side
        """
        if in_flanks is None:
            in_flanks = self.segment_length

        # Get base junctions from parent
        jc_start, jc_end = super().get_junctions(out_flanks=out_flanks, in_flanks=in_flanks)

        # Add TN-specific fields
        tn_side_start, tn_side_end = self.get_ref_tn_sides()

        return (
            RefTnJunction.from_other(
                jc_start,
                num=-self.tn_id * 2,  # Start: even negative
                ref_tn_side=tn_side_start,
            ),
            RefTnJunction.from_other(
                jc_end,
                num=-self.tn_id * 2 - 1,  # End: odd negative
                ref_tn_side=tn_side_end,
            ),
        )

    def get_inward_arm_by_ref_tn_side(self, ref_tn_side: RefTnSide, flank: int) -> JcArm:
        """Get inward arm by reference TN side with optional offset."""
        ref_arms = self.get_inward_arms(flanks=flank)
        arm = ref_arms[0] if ref_tn_side.side == Side.START else ref_arms[1]

        if isinstance(ref_tn_side, OffsetRefTnSide):
            return arm.shift_by_offset(ref_tn_side.offset)

        return arm


# ===== BLAST hits =====

class BlastHit(Record):
    """BLAST alignment hit record."""
    NAME: ClassVar[str] = "BLAST hits"
    CSV_FIELD_FORMATS: ClassVar[Dict[str, str]] = {
        'evalue': '.6e',  # Scientific notation for e-values
    }
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
    """Base junction record with coordinate fields only."""
    NAME: ClassVar[str] = "Junctions"

    # Arm 1 fields
    scaf1: str
    pos1: int
    dir1: Orientation
    flanking1: int   # Length of sequence flanking arm 1 (used for sequence extraction)

    # Arm 2 fields
    scaf2: str
    pos2: int
    dir2: Orientation
    flanking2: int  # Length of sequence flanking arm 2 (used for sequence extraction)

    @classmethod
    def _get_extra_fields(cls) -> List[str]:
        """Get extra fields not part of the arm coordinates."""
        arm_fields = Junction.model_fields.keys()
        return [f for f in cls.model_fields if f not in arm_fields]

    def _get_extra_kwargs(self) -> Dict[str, Any]:
        """Get extra kwargs not part of the arm coordinates."""
        return {field: getattr(self, field) for field in self._get_extra_fields()}

    def swap_sides(self: JunctionT, **kwargs) -> JunctionT:
        """Return new junction with arm 1 and arm 2 swapped."""
        jc_arms = self.get_jc_arms()
        extra_kwargs = self._get_extra_kwargs()
        extra_kwargs.update(kwargs)
        return self.from_jc_arms(jc_arms[1], jc_arms[0], **extra_kwargs)

    def get_jc_arm(self, arm: int) -> JcArm:
        """Get scaffold, position, direction, and flanking length for an arm."""
        return self.get_jc_arms()[arm - 1]

    def get_jc_arms(self) -> tuple[JcArm, JcArm]:
        """Get junction arms."""
        return (
            JcArm(scaf=self.scaf1, start=self.pos1, dir=self.dir1, flank=self.flanking1),
            JcArm(scaf=self.scaf2, start=self.pos2, dir=self.dir2, flank=self.flanking2),
        )

    @classmethod
    def from_jc_arms(cls, arm1: JcArm, arm2: JcArm, **kwargs) -> Junction:
        """Create a Junction from junction arm coordinates."""
        return cls(scaf1=arm1.scaf, pos1=arm1.start, dir1=arm1.dir, flanking1=arm1.flank,
                   scaf2=arm2.scaf, pos2=arm2.start, dir2=arm2.dir, flanking2=arm2.flank,
                   **kwargs)


class NumJunction(Junction):
    """Junction with identifier."""
    NAME: ClassVar[str] = "Numbered Junctions"
    # Junction identifier: breseq junction number (positive), or negative for reference junctions
    num: Optional[int] = None


class BreseqJunction(NumJunction):
    """Breseq junction."""
    NAME: ClassVar[str] = "Breseq junctions"
    model_config = ConfigDict(extra='allow')
    ALLOW_EXTRA: ClassVar[bool] = True


class RefTnJunction(NumJunction):
    """Synthetic junction for reference TN element.

    For RefTnJunction, arm 1 is always the TN side, arm 2 is the chromosome side.
    ref_tn_side indicates which TN boundary (START or END) this junction represents.
    """
    NAME: ClassVar[str] = "Reference TN junctions"

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Side.START
    #        |------>              arm1, flanking1 (into TN)
    #     <--|                     arm2, flanking2 (out of TN)

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Side.END
    #           <------|           arm1, flanking1 (into TN)
    #                  |-->        arm2, flanking2 (out of TN)

    ref_tn_side: RefTnSide


class TnJunction(NumJunction):
    """Junction matched to TN element(s)."""
    NAME: ClassVar[str] = "TN-associated junctions"
    ref_tn_side: Optional[RefTnSide] = None  # None for breseq junctions
    ref_tn_sides: List[OffsetRefTnSide]      # Reference TN matches with offsets
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
    # left: the junction from which we start the amplicon segment going on the forward strand (dir2 == FORWARD)
    # right: the junction at which we end the amplicon segment when we go on the forward strand (dir2 == REVERSE)

    modified from combine_ISJC_pairs.m

    """
    NAME: ClassVar[str] = "Junction Pairs"

    # Unique identifier
    pair_id: Optional[int] = None

    # Core fields: the two junctions
    tnjc_left: TnJunction  # Left junction (dir2 == FORWARD)
    tnjc_right: TnJunction  # Right junction (dir2 == REVERSE)

    # Scaffold object (not exported to CSV)
    scaffold: SeqScaffold

    # CSV export: only export derived properties, not the complex TnJunction/Scaffold objects
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = [
        'pair_id', 'jc_num_left', 'jc_num_right', 'scaf', 'left', 'right',
        'tn_ids', 'tn_offsets', 'amplicon_length'
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
                if i_tn_side.is_opposite_side(j_tn_side)]

    def _find_matching_tn_sides(self) -> List[tuple[OffsetRefTnSide, OffsetRefTnSide]]:
        """Instance method: find matching TNs using this record's junctions.

        Returns list of (left_tn_side, right_tn_side) tuples where both sides match the same TN ID.
        """
        return self.find_matching_tn_sides(self.tnjc_left.ref_tn_sides, self.tnjc_right.ref_tn_sides)

    @property
    def jc_num_left(self) -> int:
        """Junction number for left junction."""
        return self.tnjc_left.num

    @property
    def jc_num_right(self) -> int:
        """Junction number for right junction."""
        return self.tnjc_right.num

    @property
    def scaf(self) -> str:
        """Scaffold name."""
        assert self.tnjc_left.scaf2 == self.tnjc_right.scaf2 == self.scaffold.scaf
        return self.scaffold.scaf

    @property
    def left(self) -> int:
        """Left position of amplicon segment on the forward strand."""
        return self.tnjc_left.pos2

    @property
    def right(self) -> int:
        """Right position of amplicon segment on the forward strand."""
        return self.tnjc_right.pos2

    @property
    def tn_ids(self) -> List[int]:
        """Matching TN element IDs."""
        return [tn_side_S.tn_id for tn_side_S, tn_side_E in self._find_matching_tn_sides()]

    @property
    def tn_offsets(self) -> List[tuple[int, int]]:
        """TN offsets, one per tn_id."""
        return [(tn_side_S.offset, tn_side_E.offset) for tn_side_S, tn_side_E in self._find_matching_tn_sides()]

    @property
    def amplicon_length(self) -> int:
        """Amplicon length computed from segment scaffold."""
        return self.get_segment_scaffold().segment_length

    @property
    def span_origin(self) -> bool:
        """True if amplicon segment spans the origin of the scaffold."""
        span_origin = self.left > self.right
        assert span_origin == self.get_segment_scaffold().span_origin
        return span_origin

    def get_segment_scaffold(self) -> SeqSegmentScaffold:
        """Get SegmentScaffold for this amplicon segment.

        Returns a SeqSegmentScaffold with start=left, end=right, orientation=FORWARD
        that provides properties: segment_length, span_origin.
        """
        return SeqSegmentScaffold.from_other(
            self.scaffold,
            start=self.left,
            end=self.right,
            orientation=Orientation.FORWARD,  # by definition, the amplicon segment is on the forward strand
        )

    def get_outward_arms(self, flank: int) -> tuple[JcArm, JcArm]:
        """Get chromosome arms flanking the amplicon (outside boundaries)."""
        return self.get_segment_scaffold().get_outward_arms(flank)

    def get_inward_arms(self, flank: int) -> tuple[JcArm, JcArm]:
        """Get chromosome arms at the amplicon boundaries (inside boundaries)."""
        return self.get_segment_scaffold().get_inward_arms(flank)

    def __str__(self) -> str:
        """String representation of RawTnJc2."""
        cls_name = self.__class__.__name__
        return (f"{cls_name}({self.left}-{self.right}, scaf={self.scaf}, "
                f"len={self.amplicon_length}, tn_ids={self.tn_ids})")


class CoveredTnJc2(RawTnJc2):
    """RawTnJc2 with coverage.

    Coverage values:
    None = not applicable (anc fields when we don't have an ancestor)
    np.nan = value not calculated (amplicon too short/long)
    """
    NAME: ClassVar[str] = "Covered Junction Pairs"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = RawTnJc2.CSV_EXPORT_FIELDS + [
        'iso_scaf_avg', 'iso_amplicon_avg', 'anc_scaf_avg', 'anc_amplicon_avg', 'avg_norm_cov', 'copy_number'
    ]
    CSV_FIELD_FORMATS: ClassVar[Dict[str, str]] = {
        'iso_scaf_avg': '.1f',
        'iso_amplicon_avg': '.1f',
        'anc_scaf_avg': '.1f',
        'anc_amplicon_avg': '.1f',
        'avg_norm_cov': '.1f',
        'copy_number': '.1f',
    }
    iso_scaf_avg: float
    iso_amplicon_avg: float
    anc_scaf_avg: Optional[float] = None  # None if no ancestor
    anc_amplicon_avg: Optional[float] = None  # None if no ancestor
    avg_norm_cov: Optional[float] = None  # Position-by-position ancestor-normalized average

    @field_validator('iso_scaf_avg', 'iso_amplicon_avg', mode='before')
    @classmethod
    def _none_to_nan(cls, v):
        """Convert None to np.nan for required float fields."""
        return np.nan if v is None else v

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


class SingleLocusLinkedTnJc2(CoveredTnJc2):
    """CoveredTnJc2 with structural classification (Step 8 output)."""
    NAME: ClassVar[str] = "Classified Amplicons"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = CoveredTnJc2.CSV_EXPORT_FIELDS + [
        'single_locus_left_pair_id', 'single_locus_right_pair_id', 'raw_event', 'chosen_tn_id'
    ]
    single_locus_tnjc2_left_matchings: List[tuple[SingleLocusLinkedTnJc2, Side]]
    single_locus_tnjc2_right_matchings: List[tuple[SingleLocusLinkedTnJc2, Side]]
    base_raw_event: BaseRawEvent

    @property
    def single_locus_tnjc2_matching_left(self) -> Optional[SingleLocusLinkedTnJc2]:
        """First single-locus TN junction matching the left junction."""
        return self.single_locus_tnjc2_left_matchings[0][0] if self.single_locus_tnjc2_left_matchings else None

    @property
    def single_locus_tnjc2_matching_right(self) -> Optional[SingleLocusLinkedTnJc2]:
        """First single-locus TN junction matching the right junction."""
        return self.single_locus_tnjc2_right_matchings[0][0] if self.single_locus_tnjc2_right_matchings else None

    @property
    def is_multiple_single_locus_tnjc2_left(self) -> bool:
        """Whether there are multiple single-locus matchings for the left junction."""
        # TODO: This is not used yet
        return len(self.single_locus_tnjc2_left_matchings) > 1

    @property
    def is_multiple_single_locus_tnjc2_right(self) -> bool:
        """Whether there are multiple single-locus matchings for the right junction."""
        # TODO: This is not used yet
        return len(self.single_locus_tnjc2_right_matchings) > 1

    @property
    def single_locus_left_pair_id(self) -> Optional[int]:
        """Pair ID of the single-locus TN junction matching the left junction."""
        return self.single_locus_tnjc2_matching_left.pair_id if self.single_locus_tnjc2_matching_left else None

    @property
    def single_locus_right_pair_id(self) -> Optional[int]:
        """Pair ID of the single-locus TN junction matching the right junction."""
        return self.single_locus_tnjc2_matching_right.pair_id if self.single_locus_tnjc2_matching_right else None

    @property
    def raw_event(self) -> RawEvent:
        """Raw event classification."""
        if self.base_raw_event == BaseRawEvent.REFERENCE:
            return RawEvent.REFERENCE
        elif self.base_raw_event == BaseRawEvent.TRANSPOSITION:
            return RawEvent.TRANSPOSITION
        assert self.base_raw_event == BaseRawEvent.LOCUS_JOINING
        has_match_left = self.single_locus_tnjc2_matching_left is not None
        has_match_right = self.single_locus_tnjc2_matching_right is not None
        if has_match_left and has_match_right:
            return RawEvent.FLANKED
        elif has_match_left:
            return RawEvent.HEMI_FLANKED_LEFT
        elif has_match_right:
            return RawEvent.HEMI_FLANKED_RIGHT
        else:
            return RawEvent.UNFLANKED

    @property
    def chosen_tn_id(self) -> Optional[TnId]:
        """Selected TN for analysis."""
        # Start with TN IDs from this tnjc2
        tn_id_set = set(self.tn_ids)

        # Intersect with matching left/right tnjc2 if exists
        if self.single_locus_tnjc2_matching_left is not None:
            tn_id_set &= set(self.single_locus_tnjc2_matching_left.tn_ids)
        if self.single_locus_tnjc2_matching_right is not None:
            tn_id_set &= set(self.single_locus_tnjc2_matching_right.tn_ids)

        # Return first if available, otherwise None
        return list(tn_id_set)[0] if tn_id_set else None

    def get_sides_of_chosen_tn(self) -> tuple[Optional[OffsetRefTnSide], Optional[OffsetRefTnSide]]:
        """Get sides of chosen TN (left and right amplicon sides)."""
        chosen_id = self.chosen_tn_id
        if chosen_id is None:
            return None, None
        matching_tns = self._find_matching_tn_sides()
        for tn_side_left_amplicon, tn_side_right_amplicon in matching_tns:
            if tn_side_left_amplicon.tn_id == chosen_id:
                return tn_side_left_amplicon, tn_side_right_amplicon
        assert False


class SynJctsTnJc2(SingleLocusLinkedTnJc2):
    """Candidate with synthetic junction folder names."""
    NAME: ClassVar[str] = "Synthetic Junction Amplicons"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = SingleLocusLinkedTnJc2.CSV_EXPORT_FIELDS + [
        'analysis_dir', 'analysis_dir_anc'
    ]
    analysis_dir: str
    analysis_dir_anc: Optional[str] = None

    def analysis_dir_name(self, *, is_ancestor: bool = False) -> str:
        """Return the analysis directory name (ancestor-aware)."""
        return self.analysis_dir_anc if is_ancestor else self.analysis_dir

    def analysis_dir_path(self, base_dir: Path, *, is_ancestor: bool = False) -> Path:
        """Directory path for this candidate under base_dir/junctions."""
        return base_dir / "junctions" / self.analysis_dir_name(is_ancestor=is_ancestor)

    def fasta_path(self, base_dir: Path, *, is_ancestor: bool = False) -> Path:
        """Path to junctions FASTA for this candidate."""
        return self.analysis_dir_path(base_dir, is_ancestor=is_ancestor) / "junctions.fasta"

    def bam_path(self, base_dir: Path, *, is_ancestor: bool = False) -> Path:
        """Path to BAM for this candidate."""
        return self.analysis_dir_path(base_dir, is_ancestor=is_ancestor) / "sorted.bam"

    def alignment_plot_path(self, base_dir: Path, *, is_ancestor: bool = False) -> Path:
        """Path to alignment distribution plot PNG for this candidate."""
        return self.analysis_dir_path(base_dir, is_ancestor=is_ancestor) / f"alignment_distribution.png"


class AnalyzedTnJc2(SynJctsTnJc2):
    """Candidate with junction coverage analysis (Step 12 output).

    Junction coverage fields depend on run type:
    - anc_path=None: jc_cov only, jc_cov_anc is None
    - anc_path=set: both jc_cov and jc_cov_anc present
    """
    NAME: ClassVar[str] = "Analyzed Amplicons"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = SingleLocusLinkedTnJc2.CSV_EXPORT_FIELDS + [
        'jc_cov_vector', 'jc_cov_anc_vector'
    ]
    # Junction coverage: JunctionType -> JunctionReadCounts
    jc_cov: Dict[JunctionType, JunctionReadCounts]
    jc_cov_anc: Optional[Dict[JunctionType, JunctionReadCounts]] = None

    # For each junction type, whether it is covered (True), ambiguous (None), or not covered (False)
    jc_calls: Dict[JunctionType, Optional[bool]] = None
    jc_calls_anc: Optional[Dict[JunctionType, Optional[bool]]] = None

    @property
    def jc_cov_vector(self) -> List[tuple[int, int, int]]:
        """Junction coverage vector."""
        return [(self.jc_cov[jt].left, self.jc_cov[jt].spanning, self.jc_cov[jt].right)
                for jt in JunctionType.sorted()]

    @property
    def jc_cov_anc_vector(self) -> Optional[List[tuple[int, int, int]]]:
        """Ancestor junction coverage vector."""
        if self.jc_cov_anc is None:
            return None
        return [(self.jc_cov_anc[jt].left, self.jc_cov_anc[jt].spanning, self.jc_cov_anc[jt].right)
                for jt in JunctionType.sorted()]


class ClassifiedTnJc2(AnalyzedTnJc2):
    """AnalyzedTnJc2 with architecture/event classification (Step 13 output)."""
    NAME: ClassVar[str] = "Classified Amplicons"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = AnalyzedTnJc2.CSV_EXPORT_FIELDS + [
        'event'
    ]
    # Architecture classification
    isolate_architecture: RawEvent
    ancestor_architecture: Optional[RawEvent] = None  # only when anc_path is set

    # Final event classification
    event_modifiers: List[EventModifier]    # de novo left/right, ancestral, etc.

    @property
    def event(self) -> str:
        """Full event description derived from architecture and modifiers."""
        modifiers = [m.value for m in self.event_modifiers if m]
        return f"{self.isolate_architecture.value} ({', '.join(modifiers)})"


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

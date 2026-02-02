"""TN junction pairs and amplicon records for AmpliFinder."""
from __future__ import annotations
import numpy as np

from pathlib import Path
from typing import ClassVar, Dict, List, NamedTuple, Optional
from pydantic import field_validator

from amplifinder.records.base_records import Record
from amplifinder.data_types.basic_enums import Side, Orientation
from amplifinder.data_types.events import BaseEvent, Architecture, EventDescriptor
from amplifinder.data_types.jc_types import JcCall, JunctionType
from amplifinder.data_types.read_types import JunctionReadCounts
from amplifinder.data_types.ref_tn import TnId, OffsetRefTnSide, TnJunction
from amplifinder.data_types.scaffold import SeqScaffold, SeqSegmentScaffold
from amplifinder.data_types.junctions import JcArm


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

    # Base event classification
    base_event: BaseEvent

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


class TnJc2AndSide(NamedTuple):
    tnjc2: SingleLocusLinkedTnJc2
    side: Side


class SingleLocusLinkedTnJc2(RawTnJc2):
    """RawTnJc2 with structural classification (Step 8 output)."""
    NAME: ClassVar[str] = "Classified Junction Pairs"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = RawTnJc2.CSV_EXPORT_FIELDS + [
        'single_locus_left_pair_id', 'single_locus_right_pair_id', 'raw_event', 'chosen_tn_id'
    ]
    single_locus_tnjc2_left_matchings: List[TnJc2AndSide]
    single_locus_tnjc2_right_matchings: List[TnJc2AndSide]

    @property
    def single_locus_tnjc2_and_side_matching_left(self) -> Optional[TnJc2AndSide]:
        """First single-locus TN junction matching the left junction."""
        return self.single_locus_tnjc2_left_matchings[0] if self.single_locus_tnjc2_left_matchings else None

    @property
    def single_locus_tnjc2_and_side_matching_right(self) -> Optional[TnJc2AndSide]:
        """First single-locus TN junction matching the right junction."""
        return self.single_locus_tnjc2_right_matchings[0] if self.single_locus_tnjc2_right_matchings else None

    def get_matching_single_locus_tnjc2_and_side(self, side: Side) -> Optional[TnJc2AndSide]:
        """Get the single-locus TN junction matching the given side."""
        return self.single_locus_tnjc2_and_side_matching_left if side == Side.LEFT \
            else self.single_locus_tnjc2_and_side_matching_right

    @property
    def single_locus_left_pair_id(self) -> Optional[int]:
        """Pair ID of the single-locus TN junction matching the left junction."""
        match = self.single_locus_tnjc2_and_side_matching_left
        return match.tnjc2.pair_id if match is not None else None

    @property
    def single_locus_right_pair_id(self) -> Optional[int]:
        """Pair ID of the single-locus TN junction matching the right junction."""
        match = self.single_locus_tnjc2_and_side_matching_right
        return match.tnjc2.pair_id if match is not None else None

    @property
    def raw_event(self) -> Architecture:
        """Raw event classification."""
        if self.base_event == BaseEvent.REFERENCE_TN:
            return Architecture.REFERENCE_TN
        elif self.base_event == BaseEvent.TRANSPOSITION:
            return Architecture.TRANSPOSITION
        assert self.base_event == BaseEvent.LOCUS_JOINING
        has_match_left = self.single_locus_tnjc2_and_side_matching_left is not None
        has_match_right = self.single_locus_tnjc2_and_side_matching_right is not None
        if has_match_left and has_match_right:
            return Architecture.FLANKED
        elif has_match_left:
            return Architecture.HEMI_FLANKED_LEFT
        elif has_match_right:
            return Architecture.HEMI_FLANKED_RIGHT
        else:
            return Architecture.UNFLANKED

    @property
    def chosen_tn_id(self) -> Optional[TnId]:
        """Selected TN for analysis."""

        # TODO: Is the logic below correct?

        # For each junction, get the reference TN side ID if it it is a reference
        # TN junction (None for breseq junctions)
        ref_tn_id_left = self.tnjc_left.ref_tn_side.tn_id if self.tnjc_left.is_ref_tn_junction() else None
        ref_tn_id_right = self.tnjc_right.ref_tn_side.tn_id if self.tnjc_right.is_ref_tn_junction() else None

        if ref_tn_id_left is None and ref_tn_id_right is not None:
            # The right junction is a reference TN junction. Left is de novo.
            # We assume that the chosen TN is the one on the right side (i.e. the reference TN)
            return ref_tn_id_right

        if ref_tn_id_left is not None and ref_tn_id_right is None:
            # The left junction is a reference TN junction. Right is de novo.
            # We assume that the chosen TN is the one on the left side (i.e. the reference TN)
            return ref_tn_id_left

        if ref_tn_id_left is not None and ref_tn_id_right is not None:
            # Both junctions are sides of a reference TN
            # We are looking at an amplification between two reference TNs
            # We assume they are very similar, and arbitrarily choose one of them (left)
            return ref_tn_id_left

        # Intersect with matching left/right single-locus tnjc2 if exists
        tn_id_set = set(self.tn_ids)  # Start with TN IDs from this tnjc2
        if self.single_locus_tnjc2_and_side_matching_left is not None:
            # Intersect with matching left single-locus tnjc2
            tn_id_set &= set(self.single_locus_tnjc2_and_side_matching_left.tnjc2.tn_ids)
        if self.single_locus_tnjc2_and_side_matching_right is not None:
            # Intersect with matching right single-locus tnjc2
            tn_id_set &= set(self.single_locus_tnjc2_and_side_matching_right.tnjc2.tn_ids)

        # If we have matches, return one of them arbitrarily
        if tn_id_set:
            return tn_id_set.pop()  # Arbitrarily choose one of the matches

        # Otherwise return None (no matches)
        return None

    def get_sides_of_chosen_tn(self) -> tuple[Optional[OffsetRefTnSide], Optional[OffsetRefTnSide]]:
        """Get sides of chosen TN (left and right amplicon sides)."""
        chosen_id = self.chosen_tn_id
        if chosen_id is None:
            raise ValueError("This tnjc2 does not have a chosen tn")
        matching_tns = self._find_matching_tn_sides()
        for tn_side_left_amplicon, tn_side_right_amplicon in matching_tns:
            if tn_side_left_amplicon.tn_id == chosen_id:
                return tn_side_left_amplicon, tn_side_right_amplicon
        assert False


class CoveredTnJc2(SingleLocusLinkedTnJc2):
    """SingleLocusLinkedTnJc2 with coverage.

    Coverage values:
    None = not applicable (anc fields when we don't have an ancestor)
    np.nan = value not calculated (amplicon too short/long)
    """
    NAME: ClassVar[str] = "Covered Junction Pairs"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = SingleLocusLinkedTnJc2.CSV_EXPORT_FIELDS + [
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


class SynJctsTnJc2(CoveredTnJc2):
    """Candidate with synthetic junction folder names."""
    NAME: ClassVar[str] = "Junction Pairs with Synthetic Junctions"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = CoveredTnJc2.CSV_EXPORT_FIELDS + [
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


class AnalyzedTnJc2(SynJctsTnJc2):
    """Candidate with junction coverage analysis (Step 12 output).

    Junction coverage fields depend on run type:
    - anc_fastq_path=None: jc_cov only, jc_cov_anc is None
    - anc_fastq_path=set: both jc_cov and jc_cov_anc present
    """
    NAME: ClassVar[str] = "Junction Pairs with Synthetic Junctions Coverage Analysis"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = CoveredTnJc2.CSV_EXPORT_FIELDS + [
        'jc_cov_vector', 'jc_cov_anc_vector'
    ]
    # Junction coverage: JunctionType -> JunctionReadCounts
    jc_covs: Dict[JunctionType, JunctionReadCounts]
    jc_covs_anc: Optional[Dict[JunctionType, JunctionReadCounts]] = None

    # For each junction type, whether it is covered (True), ambiguous (None), or not covered (False)
    jc_calls: Optional[Dict[JunctionType, JcCall]] = None
    jc_calls_anc: Optional[Dict[JunctionType, JcCall]] = None

    @property
    def jc_cov_vector(self) -> List[tuple[int, int, int]]:
        """Junction coverage vector."""
        return [(self.jc_covs[jt].left, self.jc_covs[jt].spanning, self.jc_covs[jt].right)
                for jt in JunctionType]

    @property
    def jc_cov_anc_vector(self) -> Optional[List[tuple[int, int, int]]]:
        """Ancestor junction coverage vector."""
        if self.jc_covs_anc is None:
            return None
        return [(self.jc_covs_anc[jt].left, self.jc_covs_anc[jt].spanning, self.jc_covs_anc[jt].right)
                for jt in JunctionType]


class ClassifiedTnJc2(AnalyzedTnJc2):
    """AnalyzedTnJc2 with architecture/event classification (Step 13 output)."""
    NAME: ClassVar[str] = "Classified Architectures"
    CSV_EXPORT_FIELDS: ClassVar[List[str]] = AnalyzedTnJc2.CSV_EXPORT_FIELDS + [
        'event_str', 'iso_architecture', 'anc_architecture', 'event_descriptors'
    ]
    # Architecture classification
    iso_architecture: Architecture
    anc_architecture: Optional[Architecture] = None  # only when anc_fastq_path is set

    # Event descriptors
    event_descriptors: List[EventDescriptor]

    @property
    def event_str(self) -> str:
        """Full event description derived from architecture and descriptors."""
        iso_architecture_str = self.iso_architecture.description
        if self.event_descriptors:
            descriptors_str = ', '.join(d.value for d in self.event_descriptors)
            return f"{iso_architecture_str} ({descriptors_str})"
        else:
            return iso_architecture_str

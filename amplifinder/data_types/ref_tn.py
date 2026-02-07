"""Reference TN elements and related types for AmpliFinder."""
from __future__ import annotations

from typing import ClassVar, List, Optional

from amplifinder.records.base_records import Record
from amplifinder.data_types.basic_enums import Terminal
from amplifinder.data_types.scaffold import SegmentScaffold
from amplifinder.data_types.junctions import JcArm, NumJunction

TnId = int


class RefTnSide(Record):
    """A reference TN element side."""
    NAME: ClassVar[str] = "Reference TN sides"
    ref_tn: RefTn
    side: Terminal

    @property
    def tn_id(self) -> TnId:
        """Get TN ID from the reference TN."""
        return self.ref_tn.tn_id

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

    def __hash__(self) -> int:
        """Make RefTn hashable based on tn_id for use in sets."""
        return hash(self.tn_id)

    def __eq__(self, other) -> bool:
        """Compare RefTn objects by tn_id."""
        if not isinstance(other, RefTn):
            return False
        return self.tn_id == other.tn_id

    def get_ref_tn_sides(self) -> tuple[RefTnSide, RefTnSide]:
        """Get start and end sides of the TN."""
        return (
            RefTnSide(ref_tn=self, side=Terminal.START),
            RefTnSide(ref_tn=self, side=Terminal.END),
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
        arm = ref_arms[0] if ref_tn_side.side == Terminal.START else ref_arms[1]

        if isinstance(ref_tn_side, OffsetRefTnSide):
            return arm.shift_by_offset(ref_tn_side.offset)

        return arm


class RefTnJunction(NumJunction):
    """Synthetic junction for reference TN element.

    For RefTnJunction, arm 1 is always the TN side, arm 2 is the chromosome side.
    ref_tn_side indicates which TN boundary (START or END) this junction represents.
    """
    NAME: ClassVar[str] = "Reference TN junctions"

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Terminal.START
    #        |------>              arm1, flanking1 (into TN)
    #     <--|                     arm2, flanking2 (out of TN)

    #   chr      TN       chr
    # ~~~~~~~|>>>>>>>>>|~~~~~~~    ref_tn_side.side == Terminal.END
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

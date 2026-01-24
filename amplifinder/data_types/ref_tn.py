
TnId = int


class RefTnSide(Record):
    """A reference TN element side."""
    NAME: ClassVar[str] = "Reference TN sides"
    tn_id: TnId
    side: Terminal

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
            RefTnSide(tn_id=self.tn_id, side=Terminal.START),
            RefTnSide(tn_id=self.tn_id, side=Terminal.END),
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


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


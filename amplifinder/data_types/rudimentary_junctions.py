"""Rudimentary junction values for creating synthetic junctions."""
from typing import NamedTuple

from amplifinder.data_types.basic_enums import Orientation
from amplifinder.data_types.jc_types import JunctionType
from amplifinder.data_types.junctions import JcArm, Junction
from amplifinder.data_types.ref_tn import Terminal


class RudimentaryJunctionValues(NamedTuple):
    """
    Primitive parameters for creating synthetic junctions.
    Encodes the minimal parameter set for creating synthetic junctions.
    Guarentees a 1-to-1 mapping between the fasta file name and the junction sequences.
    """
    amp_left_pos: int
    amp_right_pos: int
    amp_scaf: str
    tn_start_pos: int
    tn_end_pos: int
    tn_scaf: str
    tn_side_left_amp_side: Terminal
    flank: int

    # chromosome arm start positions are indicated as offsets from the amplicon
    # 0 is precisly flanking the amplicon
    # ~~~~~~~~~~~~~~~~||===================...
    #            <----0
    chr_left_pos_offset: int
    chr_right_pos_offset: int
    chr_left_pos_for_tn_offset: int
    chr_right_pos_for_tn_offset: int

    def create_syn_junctions(self) -> dict[JunctionType, Junction]:
        """Create 7 synthetic junctions from primitive coordinates and strands."""
        # Amplicon inward arms (left/right)
        amp_left_arm = JcArm(scaf=self.amp_scaf, start=self.amp_left_pos, dir=Orientation.FORWARD, flank=self.flank)
        amp_right_arm = JcArm(scaf=self.amp_scaf, start=self.amp_right_pos, dir=Orientation.REVERSE, flank=self.flank)

        # Chromosome outward arms; mirror of amplicon if not flanked by ancestral TN
        chr_left_arm = amp_left_arm.mirror().shift_by_offset(self.chr_left_pos_offset)
        chr_right_arm = amp_right_arm.mirror().shift_by_offset(self.chr_right_pos_offset)

        # Chromosome outward arms for TN
        chr_left_arm_for_tn = amp_left_arm.mirror().shift_by_offset(self.chr_left_pos_for_tn_offset)
        chr_right_arm_for_tn = amp_right_arm.mirror().shift_by_offset(self.chr_right_pos_for_tn_offset)

        # TN inward arms
        tn_orientation = Orientation.FORWARD if self.tn_start_pos < self.tn_end_pos else Orientation.REVERSE
        tn_start_arm = JcArm(scaf=self.tn_scaf, start=self.tn_start_pos, dir=tn_orientation, flank=self.flank)
        tn_end_arm = JcArm(scaf=self.tn_scaf, start=self.tn_end_pos, dir=tn_orientation.opposite(), flank=self.flank)

        if self.tn_side_left_amp_side == Terminal.START:
            tn_left_arm, tn_right_arm = tn_end_arm, tn_start_arm
        else:
            tn_left_arm, tn_right_arm = tn_start_arm, tn_end_arm

        return build_7_junctions_from_8_arms(
            chr_left_arm_for_tn, chr_right_arm_for_tn, chr_left_arm, chr_right_arm,
            amp_left_arm, amp_right_arm, tn_left_arm, tn_right_arm)

    def get_name(self) -> str:
        side_str = "S" if self.tn_side_left_amp_side == Terminal.START else "E"
        chr_left_pos_str = f"L{self.chr_left_pos_offset:+}L{self.chr_left_pos_for_tn_offset:+}"
        chr_right_pos_str = f"R{self.chr_right_pos_offset:+}R{self.chr_right_pos_for_tn_offset:+}"
        return (f"jc_{self.amp_scaf}_{self.amp_left_pos}-{self.amp_right_pos}_"
                f"{self.tn_scaf}_{self.tn_start_pos}-{self.tn_end_pos}_"
                f"{chr_left_pos_str}_{chr_right_pos_str}"
                f"{side_str}_{self.flank}bp")


def build_7_junctions_from_8_arms(
    chr_left_for_tn: JcArm, chr_right_for_tn: JcArm,
    chr_left: JcArm, chr_right: JcArm,
    amp_left: JcArm, amp_right: JcArm,
    tn_left: JcArm, tn_right: JcArm,
) -> dict[JunctionType, Junction]:
    """Create Junction objects for each junction type."""
    return {
        JunctionType.CHR_AMP: Junction.from_jc_arms(chr_left, amp_left),
        JunctionType.CHR_TN: Junction.from_jc_arms(chr_left_for_tn, tn_left),
        JunctionType.AMP_TN: Junction.from_jc_arms(amp_right, tn_left),
        JunctionType.AMP_AMP: Junction.from_jc_arms(amp_right, amp_left),
        JunctionType.TN_AMP: Junction.from_jc_arms(tn_right, amp_left),
        JunctionType.TN_CHR: Junction.from_jc_arms(tn_right, chr_right_for_tn),
        JunctionType.AMP_CHR: Junction.from_jc_arms(amp_right, chr_right),
    }

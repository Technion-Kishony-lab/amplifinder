"""
Create synthetic junction sequences for alignment.

Amplicon structure:

~~~~~~~~~>>>======>>>======>>>~~~~~~~~~

legend:
~~~    chromosome
>>>    IS
====== cassette (amplicon)
(1) ~~-==  left reference (chromosome-cassette)
(2) ~~->>  left IS transposition (chromosome-IS)
(3) ==->>  left of mid IS (cassette-IS)
(4) ==-==  lost IS (cassette-cassette, no IS)
(5) >>-==  right of mid IS (IS-cassette)
(6) >>-~~  right IS transposition (IS-chromosome)
(7) ==-~~  right reference (cassette-chromosome)
"""

from pathlib import Path
from typing import Optional, NamedTuple

from amplifinder.data_types import (
    BaseEvent, RecordTypedDf, CoveredTnJc2, Side, SynJctsTnJc2, Genome,
    JunctionType, RefTn, Junction, Terminal, JcArm, Orientation
)

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.utils.file_utils import ensure_dir, remove_file_or_dir
from amplifinder.utils.fasta import write_fasta
from amplifinder.utils.file_lock import locked_resource


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

        return _build_7_junctions_from_8_arms(
            chr_left_arm_for_tn, chr_right_arm_for_tn, chr_left_arm, chr_right_arm, amp_left_arm, amp_right_arm, tn_left_arm, tn_right_arm)

    def get_name(self) -> str:
        side_str = "S" if self.tn_side_left_amp_side == Terminal.START else "E"
        chr_left_pos_str = f"L{self.chr_left_pos_offset:+}L{self.chr_left_pos_for_tn_offset:+}"
        chr_right_pos_str = f"R{self.chr_right_pos_offset:+}R{self.chr_right_pos_for_tn_offset:+}"
        return (f"jc_{self.amp_scaf}_{self.amp_left_pos}-{self.amp_right_pos}_"
                f"{self.tn_scaf}_{self.tn_start_pos}-{self.tn_end_pos}_"
                f"{chr_left_pos_str}_{chr_right_pos_str}"
                f"{side_str}_{self.flank}bp")


def _build_7_junctions_from_8_arms(
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


def _handle_flanked_side(
    tnjc2: CoveredTnJc2, side: Side, jc_arm_len: int, chr_arm: JcArm
) -> tuple[JcArm, JcArm]:
    """
    Handle tnjc2 flanked by ancestral or transposed TN on one side.
    Returns: (chr_arm_for_tn, chr_arm, is_ref_tn)
    """
    paired = tnjc2.get_matching_single_locus_tnjc2_and_side(side)
    if paired is None:
        # Not paired with a single-locus TN junction
        chr_arm_for_tn = chr_arm
    else:
        assert paired.side == side
        arm_index = 1 if side == Side.LEFT else 0
        chr_arm_for_tn = paired.tnjc2.get_inward_arms(flank=jc_arm_len)[arm_index]
        if paired.tnjc2.base_event == BaseEvent.REFERENCE_TN:
            chr_arm = chr_arm_for_tn
    return chr_arm_for_tn, chr_arm


def create_synthetic_junctions_and_name(
    tnjc2: CoveredTnJc2, jc_arm_len: int, ref_tns: RecordTypedDf[RefTn]
) -> tuple[dict[JunctionType, Junction], str]:

    # Get TN sides that S and E junctions connect to (with offsets)
    tn_side_left_amp, tn_side_right_amp = tnjc2.get_sides_of_chosen_tn()

    # Get the chosen RefTn from chosen_tn_id (O(1) lookup using indexed ref_tns)
    chosen_tn_id = tnjc2.chosen_tn_id
    ref_tn = ref_tns[chosen_tn_id]
    assert ref_tn.tn_id == chosen_tn_id

    # Get inward arms for Amplicon
    amp_left_arm, amp_right_arm = tnjc2.get_inward_arms(flank=jc_arm_len)

    # Get chromosome arms (outward amplicon arms)
    chr_left_arm, chr_right_arm = tnjc2.get_outward_arms(flank=jc_arm_len)

    # Handle tnjc2 flanked by ancestral or transposed TN
    chr_left_arm_for_tn, chr_left_arm = _handle_flanked_side(tnjc2, Side.LEFT, jc_arm_len, chr_left_arm)
    chr_right_arm_for_tn, chr_right_arm = _handle_flanked_side(tnjc2, Side.RIGHT, jc_arm_len, chr_right_arm)

    # Get inward arms for RefTn with offset adjustments via ref_tn_side
    # The right-side of the TN is one that connects to the left-side of the amplicon
    tn_right_arm = ref_tn.get_inward_arm_by_ref_tn_side(tn_side_left_amp, jc_arm_len)
    tn_left_arm = ref_tn.get_inward_arm_by_ref_tn_side(tn_side_right_amp, jc_arm_len)

    # Create Junction objects for each junction type
    direct_jc = _build_7_junctions_from_8_arms(
        chr_left_arm_for_tn, chr_right_arm_for_tn, chr_left_arm, chr_right_arm, amp_left_arm, amp_right_arm, tn_left_arm, tn_right_arm)

    # Create rudimentary junction values for naming
    tn_side_left_amp_side = tn_side_left_amp.side
    rudimentary = RudimentaryJunctionValues(
        amp_left_pos=amp_left_arm.start, amp_right_pos=amp_right_arm.start, amp_scaf=amp_left_arm.scaf,
        tn_start_pos=tn_left_arm.start if tn_side_left_amp_side == Terminal.END else tn_right_arm.start,
        tn_end_pos=tn_right_arm.start if tn_side_left_amp_side == Terminal.END else tn_left_arm.start,
        tn_scaf=tn_left_arm.scaf,
        tn_side_left_amp_side=tn_side_left_amp_side,
        chr_left_pos_offset=amp_left_arm.mirror().get_distance_to(chr_left_arm.start),
        chr_right_pos_offset=amp_right_arm.mirror().get_distance_to(chr_right_arm.start),
        chr_left_pos_for_tn_offset=amp_left_arm.mirror().get_distance_to(chr_left_arm_for_tn.start),
        chr_right_pos_for_tn_offset=amp_right_arm.mirror().get_distance_to(chr_right_arm_for_tn.start),
        flank=jc_arm_len,
    )

    assert direct_jc == rudimentary.create_syn_junctions()

    return direct_jc, rudimentary.get_name()


class CreateSyntheticJunctionsStep(RecordTypedDfStep[SynJctsTnJc2]):
    """
    Creates 7 junction sequences per candidate for read alignment analysis.
    """
    NAME = "Create synthetic junctions"
    dir_field: str = 'analysis_dir'
    is_ancestor: bool = False

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[CoveredTnJc2],
        genome: Genome,
        ref_tns: RecordTypedDf[RefTn],
        output_dir: Path,
        jc_arm_len: int = 150,
        force: Optional[bool] = None,
        csv_output_dir: Optional[Path] = None,
    ):
        self.genome = genome
        self.ref_tns = ref_tns
        self.jc_arm_len = jc_arm_len
        self._artifact_dir = Path(output_dir)
        self._tnjc2s_and_junctions: list[tuple[SynJctsTnJc2, dict[JunctionType, Junction]]] = []
        for tnjc2 in filtered_tnjc2s:
            junctions, analysis_dir_name = create_synthetic_junctions_and_name(
                tnjc2=tnjc2, jc_arm_len=self.jc_arm_len, ref_tns=self.ref_tns
            )
            tnjc2 = SynJctsTnJc2.from_other(tnjc2, **{self.dir_field: analysis_dir_name})
            self._tnjc2s_and_junctions.append((tnjc2, junctions))

        artifact_files = [
            tnjc2.fasta_path(output_dir, is_ancestor=self.is_ancestor)
            for tnjc2, _ in self._tnjc2s_and_junctions]
        super().__init__(
            input_files=[genome.fasta_path],
            artifact_files=artifact_files,
            output_dir=csv_output_dir or output_dir,
            force=force,
        )

    def _generate_artifacts(self) -> None:
        """Create junction FASTA files for each candidate."""
        for tnjc2, junctions in self._tnjc2s_and_junctions:
            fasta_path = tnjc2.fasta_path(self._artifact_dir, is_ancestor=self.is_ancestor)

            # Lock per-junction for ancestor (None for isolate = no lock)
            lock_path = fasta_path.parent if self.is_ancestor else None

            with locked_resource(lock_path, "junction_fasta", timeout=300):
                if fasta_path.exists() and not self.force:
                    continue
                if self.force:
                    remove_file_or_dir(fasta_path)

                ensure_dir(fasta_path.parent)
                sequences = [
                    (jtype.name, self.genome.get_junction_sequence_arm1_to_arm2(jc))
                    for jtype, jc in junctions.items()
                ]
                write_fasta(sequences, fasta_path)

    def _calculate_output(self) -> RecordTypedDf[SynJctsTnJc2]:
        """Return records with analysis dirs (artifacts are FASTA files)."""
        tnjc2s = [tnjc2 for tnjc2, _ in self._tnjc2s_and_junctions]
        return RecordTypedDf.from_records(tnjc2s, SynJctsTnJc2)

    def _artifact_labels(self) -> list[str]:
        """Summarize outputs as count."""
        n = len(self._tnjc2s_and_junctions)
        return [f"{n} fasta files of synthetic junctions"]


class AncCreateSyntheticJunctionsStep(CreateSyntheticJunctionsStep):
    """Creates synthetic junctions for ancestor with its own junction length."""
    NAME = "Create synthetic junctions (ancestor)"
    dir_field: str = 'analysis_dir_anc'
    is_ancestor: bool = True

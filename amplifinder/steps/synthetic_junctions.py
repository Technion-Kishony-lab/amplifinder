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

from amplifinder.data_types import BaseRawEvent, RecordTypedDf, SingleLocusLinkedTnJc2, SynJctsTnJc2, Genome, JunctionType, RefTn, \
    Junction, Side, JcArm, Orientation

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.utils.file_utils import ensure_dir, remove_file_or_dir
from amplifinder.utils.fasta import write_fasta


class RudimentaryJunctionValues(NamedTuple):
    amp_left_pos: int
    amp_right_pos: int
    amp_scaf: str
    tn_start_pos: int
    tn_end_pos: int
    tn_scaf: str
    tn_side_left_amp_side: Side
    flank: int
    # Optional: chr positions in case of flanked TNJC2. Otherwise, we use the mirror of the amplicon arms
    chr_left_pos: Optional[int] = None
    chr_right_pos: Optional[int] = None

    def create_syn_junctions(self) -> dict[JunctionType, Junction]:
        """Create 7 synthetic junctions from primitive coordinates and strands."""
        # Amplicon inward arms (left/right)
        amp_left_arm = JcArm(scaf=self.amp_scaf, start=self.amp_left_pos, dir=Orientation.FORWARD, flank=self.flank)
        amp_right_arm = JcArm(scaf=self.amp_scaf, start=self.amp_right_pos, dir=Orientation.REVERSE, flank=self.flank)

        # Chromosome outward arms; mirror of amplicon if not flanked by ancestral TN
        if self.chr_left_pos is None:
            chr_left_arm = amp_left_arm.mirror()
        else:
            chr_left_arm = JcArm(scaf=self.amp_scaf, start=self.chr_left_pos, dir=Orientation.REVERSE, flank=self.flank)

        if self.chr_right_pos is None:
            chr_right_arm = amp_right_arm.mirror()
        else:
            chr_right_arm = JcArm(scaf=self.amp_scaf, start=self.chr_right_pos, dir=Orientation.FORWARD, flank=self.flank)

        # TN inward arms
        tn_orientation = Orientation.FORWARD if self.tn_start_pos < self.tn_end_pos else Orientation.REVERSE
        tn_start_arm = JcArm(scaf=self.tn_scaf, start=self.tn_start_pos, dir=tn_orientation, flank=self.flank)
        tn_end_arm = JcArm(scaf=self.tn_scaf, start=self.tn_end_pos, dir=tn_orientation.opposite(), flank=self.flank)

        if self.tn_side_left_amp_side == Side.START:
            tn_left_arm, tn_right_arm = tn_end_arm, tn_start_arm
        else:
            tn_left_arm, tn_right_arm = tn_start_arm, tn_end_arm

        return _build_7_junctions_from_6_arms(
            chr_left_arm, chr_right_arm, amp_left_arm, amp_right_arm, tn_left_arm, tn_right_arm)

    def get_name(self) -> str:
        side_str = "S" if self.tn_side_left_amp_side == Side.START else "E"
        chr_left_pos_str = f"chrL={self.chr_left_pos}_" if self.chr_left_pos is not None else ""
        chr_right_pos_str = f"chrR={self.chr_right_pos}_" if self.chr_right_pos is not None else ""
        return (f"jc_{self.amp_scaf}_{self.amp_left_pos}-{self.amp_right_pos}_"
                f"{self.tn_scaf}_{self.tn_start_pos}-{self.tn_end_pos}_"
                f"{chr_left_pos_str}{chr_right_pos_str}"
                f"{side_str}_{self.flank}bp")


def _build_7_junctions_from_6_arms(
    chr_left: JcArm, chr_right: JcArm,
    amp_left: JcArm, amp_right: JcArm,
    tn_left: JcArm, tn_right: JcArm,
) -> dict[JunctionType, Junction]:
    """Create Junction objects for each junction type."""
    return {
        JunctionType.CHR_TO_AMP_LEFT: Junction.from_jc_arms(chr_left, amp_left),
        JunctionType.CHR_TO_TN_LEFT: Junction.from_jc_arms(chr_left, tn_left),
        JunctionType.AMP_RIGHT_TO_TN_LEFT: Junction.from_jc_arms(amp_right, tn_left),
        JunctionType.AMP_RIGHT_TO_AMP_LEFT: Junction.from_jc_arms(amp_right, amp_left),
        JunctionType.TN_RIGHT_TO_AMP_LEFT: Junction.from_jc_arms(tn_right, amp_left),
        JunctionType.TN_RIGHT_TO_CHR: Junction.from_jc_arms(tn_right, chr_right),
        JunctionType.AMP_RIGHT_TO_CHR: Junction.from_jc_arms(amp_right, chr_right),
    }


def create_synthetic_junctions_and_name(
    tnjc2: SingleLocusLinkedTnJc2, jc_width: int, ref_tns: RecordTypedDf[RefTn]
) -> tuple[dict[JunctionType, Junction], str]:

    # Get TN sides that S and E junctions connect to (with offsets)
    tn_side_left_amp, tn_side_right_amp = tnjc2.get_sides_of_chosen_tn()

    # Get the chosen RefTn from chosen_tn_id (O(1) lookup using indexed ref_tns)
    chosen_tn_id = tnjc2.chosen_tn_id
    ref_tn = ref_tns[chosen_tn_id]
    assert ref_tn.tn_id == chosen_tn_id

    # Get inward arms for Amplicon
    amp_left_arm, amp_right_arm = tnjc2.get_inward_arms(flank=jc_width)

    # Get chromosome arms (outward amplicon arms)
    chr_left_arm, chr_right_arm = tnjc2.get_outward_arms(flank=jc_width)

    # Handle tnjc2 flanked by ancestral TN
    paired_left = tnjc2.single_locus_tnjc2_matching_left
    is_left_ref_tn = paired_left is not None and paired_left.raw_event == BaseRawEvent.REFERENCE    
    if is_left_ref_tn:
        chr_left_arm = paired_left.get_outward_arms(flank=jc_width)[0]

    paired_right = tnjc2.single_locus_tnjc2_matching_right
    is_right_ref_tn = paired_right is not None and paired_right.raw_event == BaseRawEvent.REFERENCE
    if is_right_ref_tn:
        chr_right_arm = paired_right.get_outward_arms(flank=jc_width)[1]

    # Get inward arms for RefTn with offset adjustments via ref_tn_side
    # The right-side of the TN is one that connects to the left-side of the amplicon
    tn_right_arm = ref_tn.get_inward_arm_by_ref_tn_side(tn_side_left_amp, jc_width)
    tn_left_arm = ref_tn.get_inward_arm_by_ref_tn_side(tn_side_right_amp, jc_width)

    # Create Junction objects for each junction type
    direct_jc = _build_7_junctions_from_6_arms(
        chr_left_arm, chr_right_arm, amp_left_arm, amp_right_arm, tn_left_arm, tn_right_arm)

    # Create rudimentary junction values for naming
    tn_side_left_amp_side = tn_side_left_amp.side
    rudimentary = RudimentaryJunctionValues(
        amp_left_pos=amp_left_arm.start, amp_right_pos=amp_right_arm.start, amp_scaf=amp_left_arm.scaf,
        tn_start_pos=tn_left_arm.start if tn_side_left_amp_side == Side.END else tn_right_arm.start,
        tn_end_pos=tn_right_arm.start if tn_side_left_amp_side == Side.END else tn_left_arm.start,
        tn_scaf=tn_left_arm.scaf,
        tn_side_left_amp_side=tn_side_left_amp_side,
        chr_left_pos=chr_left_arm.start if is_left_ref_tn else None,
        chr_right_pos=chr_right_arm.start if is_right_ref_tn else None,
        flank=jc_width,
    )

    assert direct_jc == rudimentary.create_syn_junctions()

    return direct_jc, rudimentary.get_name()


class CreateSyntheticJunctionsStep(RecordTypedDfStep[SynJctsTnJc2]):
    """
    Creates 7 junction sequences per candidate for read alignment analysis.
    """
    dir_field: str = 'analysis_dir'
    is_ancestor: bool = False

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[SingleLocusLinkedTnJc2],
        genome: Genome,
        ref_tns: RecordTypedDf[RefTn],
        output_dir: Path,
        junction_length: int = 150,
        force: Optional[bool] = None,
    ):
        self.genome = genome
        self.ref_tns = ref_tns
        self.junction_length = junction_length
        self._tnjc2s_and_junctions: list[tuple[SynJctsTnJc2, dict[JunctionType, Junction]]] = []
        for tnjc2 in filtered_tnjc2s:
            junctions, analysis_dir_name = create_synthetic_junctions_and_name(
                tnjc2=tnjc2, jc_width=self.junction_length * 2, ref_tns=self.ref_tns
            )
            tnjc2 = SynJctsTnJc2.from_other(tnjc2, **{self.dir_field: analysis_dir_name})
            self._tnjc2s_and_junctions.append((tnjc2, junctions))

        artifact_files = [
            tnjc2.fasta_path(output_dir, is_ancestor=self.is_ancestor)
            for tnjc2, _ in self._tnjc2s_and_junctions]
        super().__init__(
            input_files=[genome.fasta_path],
            artifact_files=artifact_files,
            output_dir=output_dir,
            force=force,
        )

    def _generate_artifacts(self) -> None:
        """Create junction FASTA files for each candidate."""
        for tnjc2, junctions in self._tnjc2s_and_junctions:
            fasta_path = tnjc2.fasta_path(self.output_dir, is_ancestor=self.is_ancestor)
            if self.force:
                remove_file_or_dir(fasta_path)
            if fasta_path.exists():
                continue

            ensure_dir(fasta_path.parent)
            sequences = {
                jtype.name: self.genome.get_junction_sequence_arm1_to_arm2(jc)
                for jtype, jc in junctions.items()
            }
            write_fasta(sequences, fasta_path, sort_keys=True)

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
    dir_field: str = 'analysis_dir_anc'
    is_ancestor: bool = True

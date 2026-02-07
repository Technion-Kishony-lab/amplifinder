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
from typing import Optional

from amplifinder.data_types import (
    BaseEvent, RecordTypedDf, CoveredTnJc2, Side, SynJctsTnJc2, Genome,
    JunctionType, Junction, Terminal, JcArm
)
from amplifinder.data_types.rudimentary_junctions import RudimentaryJunctionValues, build_7_junctions_from_8_arms

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.utils.file_utils import ensure_dir, remove_file_or_dir
from amplifinder.utils.fasta import write_fasta
from amplifinder.utils.file_lock import locked_resource


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
    tnjc2: CoveredTnJc2, jc_arm_len: int
) -> RudimentaryJunctionValues:

    # Get TN sides that S and E junctions connect to (with offsets)
    tn_side_left_amp, tn_side_right_amp = tnjc2.get_sides_of_chosen_tn()

    # Get the chosen RefTn directly
    ref_tn = tnjc2.chosen_tn
    assert ref_tn is not None

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
    direct_jc = build_7_junctions_from_8_arms(
        chr_left_arm_for_tn, chr_right_arm_for_tn, chr_left_arm, chr_right_arm,
        amp_left_arm, amp_right_arm, tn_left_arm, tn_right_arm)

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

    return rudimentary


class CreateSyntheticJunctionsStep(RecordTypedDfStep[SynJctsTnJc2]):
    """
    Creates 7 junction sequences per candidate for read alignment analysis.
    """
    NAME = "Create synthetic junctions"
    dir_field: str = 'rudimentary_junction_values'
    is_ancestor: bool = False

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[CoveredTnJc2],
        genome: Genome,
        output_dir: Path,
        jc_arm_len: int = 150,
        force: Optional[bool] = None,
        csv_output_dir: Optional[Path] = None,
    ):
        self.genome = genome
        self.jc_arm_len = jc_arm_len
        self._artifact_dir = Path(output_dir)
        self._tnjc2s_and_junctions: list[tuple[SynJctsTnJc2, dict[JunctionType, Junction]]] = []
        for tnjc2 in filtered_tnjc2s:
            rudimentary = create_synthetic_junctions_and_name(
                tnjc2=tnjc2, jc_arm_len=self.jc_arm_len
            )
            junctions = rudimentary.create_syn_junctions()
            tnjc2 = SynJctsTnJc2.from_other(
                tnjc2,
                **{self.dir_field: rudimentary}
            )
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
    dir_field: str = 'anc_rudimentary_junction_values'
    is_ancestor: bool = True

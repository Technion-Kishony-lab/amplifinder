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

from amplifinder.data_types import RecordTypedDf, ClassifiedTnJc2, SynJctsTnJc2, Genome, JunctionType, RefTn, \
    Junction, Side, JcArm, Orientation

from amplifinder.steps.base import Step
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

    def create_syn_junctions(self) -> dict[JunctionType, Junction]:
        """Create 7 synthetic junctions from primitive coordinates and strands."""
        # Amplicon inward arms (left/right)
        amp_left_arm = JcArm(scaf=self.amp_scaf, start=self.amp_left_pos, dir=Orientation.FORWARD, flank=self.flank)
        amp_right_arm = JcArm(scaf=self.amp_scaf, start=self.amp_right_pos, dir=Orientation.REVERSE, flank=self.flank)

        # Chromosome outward arms (mirror of amplicon)
        chr_left_arm = amp_left_arm.mirror()
        chr_right_arm = amp_right_arm.mirror()

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
        return f"jc_{self.amp_scaf}_{self.amp_left_pos}-{self.amp_right_pos}_{self.tn_scaf}_{self.tn_start_pos}-{self.tn_end_pos}_{side_str}_{self.flank}bp"


def _build_7_junctions_from_6_arms(
    chr_left: JcArm, chr_right: JcArm,
    amp_left: JcArm, amp_right: JcArm,
    tn_left: JcArm, tn_right: JcArm,
) -> dict[JunctionType, Junction]:
    """Create Junction objects for each junction type."""
    return {
        JunctionType.CHR_TO_AMP_LEFT:       Junction.from_jc_arms(chr_left,  amp_left),
        JunctionType.CHR_TO_TN_LEFT:        Junction.from_jc_arms(chr_left,  tn_left ),
        JunctionType.AMP_RIGHT_TO_TN_LEFT:  Junction.from_jc_arms(amp_right, tn_left ),
        JunctionType.AMP_RIGHT_TO_AMP_LEFT: Junction.from_jc_arms(amp_right, amp_left),
        JunctionType.TN_RIGHT_TO_AMP_LEFT:  Junction.from_jc_arms(tn_right,  amp_left),
        JunctionType.TN_RIGHT_TO_CHR:       Junction.from_jc_arms(tn_right,  chr_right),
        JunctionType.AMP_RIGHT_TO_CHR:      Junction.from_jc_arms(amp_right, chr_right),
    }


def create_synthetic_junctions(tnjc2: ClassifiedTnJc2, jc_width: int, ref_tns: RecordTypedDf[RefTn]) -> tuple[dict[JunctionType, Junction], str]:

    # Get TN sides that S and E junctions connect to (with offsets)
    tn_side_left_amp, tn_side_right_amp = tnjc2.get_sides_of_chosen_tn()
    
    # Get RefTn from chosen_tn_id (O(1) lookup using indexed ref_tns)
    chosen_tn_id = tnjc2.chosen_tn_id
    ref_tn = ref_tns[chosen_tn_id]
    assert ref_tn.tn_id == chosen_tn_id
    
    # Get inward arms for Amplicon
    amp_left_arm, amp_right_arm = tnjc2.get_inward_arms(flank=jc_width)

    # Get chromosome arms (outward amplicon arms)
    chr_left_arm, chr_right_arm = tnjc2.get_outward_arms(flank=jc_width)
   
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
        flank=jc_width,
    )

    return direct_jc, rudimentary.get_name()


class CreateSyntheticJunctionsStep(Step[RecordTypedDf[SynJctsTnJc2]]):
    """
    Creates 7 junction sequences per candidate for read alignment analysis.
    """
    dir_field: str = 'analysis_dir'

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        genome: Genome,
        ref_tns: RecordTypedDf[RefTn],
        output_dir: Path,
        junction_length: int = 150,
        force: Optional[bool] = None,
    ):
        self.filtered_tnjc2s = filtered_tnjc2s
        self.genome = genome
        self.ref_tns = ref_tns
        self.output_dir = Path(output_dir)
        self.junction_length = junction_length

        super().__init__(input_files=[genome.fasta_path], output_files=[], force=force)

    def _get_analysis_dir(self, analysis_dir: str) -> Path:
        """Get analysis directory path for a candidate."""
        return self.output_dir / "junctions" / analysis_dir

    def _calculate_output(self) -> RecordTypedDf[SynJctsTnJc2]:
        """Create synthetic junctions for each candidate."""
        syn_records: list[SynJctsTnJc2] = []
        for tnjc2 in self.filtered_tnjc2s:
            # Skip candidates with no chosen TN
            if tnjc2.chosen_tn_id is None:
                continue
                
            junctions, analysis_dir = create_synthetic_junctions(tnjc2=tnjc2, jc_width=self.junction_length * 2, ref_tns=self.ref_tns)
            analysis_path = self._get_analysis_dir(analysis_dir)

            # Extract sequences and write FASTA
            sequences = {
                jtype.name: self.genome.get_junction_sequence_arm1_to_arm2(jc)
                for jtype, jc in junctions.items()
            }
            write_fasta(sequences, analysis_path / "junctions.fasta", sort_keys=True)

            syn_records.append(
                SynJctsTnJc2.from_other(tnjc2, **{self.dir_field: analysis_dir})
            )
        return RecordTypedDf.from_records(syn_records, SynJctsTnJc2)

    def _save_output(self, output: RecordTypedDf[SynJctsTnJc2]) -> None:
        """Output already saved in _calculate_output."""
        pass

    def _output_labels(self) -> list[str]:
        """Summarize outputs as count."""
        if self.output_files:
            n = len(self.filtered_tnjc2s)
            return [f"{n} fasta files of synthetic junctions"]
        return []

    def load_outputs(self) -> RecordTypedDf[SynJctsTnJc2]:
        """Return records with analysis dirs (can't reload, must recompute)."""
        return self._calculate_output()


class AncCreateSyntheticJunctionsStep(CreateSyntheticJunctionsStep):
    """Creates synthetic junctions for ancestor with its own junction length."""
    dir_field: str = 'analysis_dir_anc'

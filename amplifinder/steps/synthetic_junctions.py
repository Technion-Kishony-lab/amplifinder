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
    RecordTypedDf, FilteredTnJc2, Genome, JunctionType, RefTn, 
    Junction, Orientation, Side, JcArm, TnJunction,
)
from amplifinder.steps.base import Step
from amplifinder.utils.fasta import write_fasta


def create_synthetic_junctions(tnjc2: FilteredTnJc2, jc_width: int, ref_tns: RecordTypedDf[RefTn]) -> dict[JunctionType, Junction]:

    # Get TN sides that S and E junctions connect to (with offsets)
    tn_side_left_amp, tn_side_right_amp = tnjc2.get_sides_of_chosen_tn()
    
    # Get RefTn from chosen_tn_id (O(1) lookup using indexed ref_tns)
    chosen_tn_id = tnjc2.chosen_tn_id
    ref_tn = ref_tns[chosen_tn_id]
    assert ref_tn.tn_id == chosen_tn_id
    
    # Get inward arms for Amplicon
    amp_left, amp_right = tnjc2.get_inward_arms(flank=jc_width)

    # Get chromosome arms (ouward amplicon arms)
    chr_left, chr_right = tnjc2.get_outward_arms(flank=jc_width)
   
    # Get inward arms for RefTn with offset adjustments via ref_tn_side
    # The right-side of the TN is one that connects to the left-side of the amplicon
    tn_right = ref_tn.get_inward_arm_by_ref_tn_side(tn_side_left_amp, jc_width)
    tn_left = ref_tn.get_inward_arm_by_ref_tn_side(tn_side_right_amp, jc_width)
    
    # Create Junction objects for each junction type
    return {
        JunctionType.CHR_TO_AMP_LEFT:       Junction.from_jc_arms(chr_left,  amp_left),
        JunctionType.CHR_TO_TN_LEFT:        Junction.from_jc_arms(chr_left,  tn_left ),
        JunctionType.AMP_RIGHT_TO_TN_LEFT:  Junction.from_jc_arms(amp_right, tn_left ),
        JunctionType.AMP_RIGHT_TO_AMP_LEFT: Junction.from_jc_arms(amp_right, amp_left),
        JunctionType.TN_RIGHT_TO_AMP_LEFT:  Junction.from_jc_arms(tn_right,  amp_left),
        JunctionType.TN_RIGHT_TO_CHR:       Junction.from_jc_arms(tn_right,  chr_right),
        JunctionType.AMP_RIGHT_TO_CHR:      Junction.from_jc_arms(amp_right, chr_right),
    }


class CreateSyntheticJunctionsStep(Step[None]):
    """
    Creates 7 junction sequences per candidate for read alignment analysis.
    """

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        genome: Genome,
        ref_tns: RecordTypedDf[RefTn],
        output_dir: Path,
        read_length: int = 150,
        force: Optional[bool] = None,
    ):
        self.filtered_tnjc2s = filtered_tnjc2s
        self.genome = genome
        self.ref_tns = ref_tns
        self.output_dir = Path(output_dir)
        self.read_length = read_length

        # Output files are per-candidate analysis directories
        self.analysis_dirs = []
        output_files = []
        for filtered_tnjc2 in filtered_tnjc2s:
            analysis_dir = self._get_analysis_dir(filtered_tnjc2)
            self.analysis_dirs.append(analysis_dir)
            output_files.append(analysis_dir / "junctions.fasta")

        super().__init__(
            input_files=[genome.fasta_path],
            output_files=output_files,
            force=force,
        )

    def _get_analysis_dir(self, filtered_tnjc2: FilteredTnJc2) -> Path:
        """Get analysis directory for a candidate."""
        return self.output_dir / "junctions" / filtered_tnjc2.analysis_dir

    def _calculate_output(self) -> None:
        """Create synthetic junctions for each candidate."""
        for tnjc2 in self.filtered_tnjc2s:
            analysis_dir = self._get_analysis_dir(tnjc2)

            junctions = create_synthetic_junctions(tnjc2=tnjc2, jc_width=self.read_length * 2, ref_tns=self.ref_tns)

            # Extract sequences and write FASTA
            sequences = {
                jtype.name: self.genome.get_junction_sequence_arm1_to_arm2(jc)
                for jtype, jc in junctions.items()
            }
            write_fasta(sequences, analysis_dir / "junctions.fasta", sort_keys=True)
        return None

    def _save_output(self, output: None) -> None:
        """Output already saved in _calculate_output."""
        pass

    def _output_labels(self) -> list[str]:
        """Summarize outputs as count."""
        if self.output_files:
            n = len(self.filtered_tnjc2s)
            return [f"{n} fasta files of synthetic junctions"]
        return []

    def load_outputs(self) -> None:
        """Junction files are side effects, nothing to return."""
        return None


class AncCreateSyntheticJunctionsStep(CreateSyntheticJunctionsStep):
    """
    Creates synthetic junctions for ancestor with its own read_length.
    Avoids copying junctions from isolate.
    """
    pass

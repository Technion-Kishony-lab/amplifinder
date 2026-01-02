"""
Create synthetic junction sequences for alignment.
Junction types (amplicon: ~~~>>>======>>>======>>>~~~):
1. ~~==  left reference (chromosome-cassette)
2. ~~>>  left IS transposition (chromosome-IS)
3. ==>>  left of mid IS (cassette-IS)
4. ====  lost IS (cassette-cassette, no IS)
5. >>==  right of mid IS (IS-cassette)
6. >>~~  right IS transposition (IS-chromosome)
7. ==~~  right reference (cassette-chromosome)
"""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, Genome, JunctionType, RefTn, 
    Junction, Orientation, Side, JcArm,
)
from amplifinder.steps.base import Step
from amplifinder.utils.file_utils import ensure_parent_dir


def create_synthetic_junctions(tnjc2: FilteredTnJc2, jc_width: int, genome: Genome) -> dict[JunctionType, Junction]:

    # Get amplicon boundaries (1-based)
    amp_scaf = tnjc2.scaf
    amp_S = tnjc2.start
    amp_E = tnjc2.end

    # Outside positions (flanking the amplicon, 1-based)
    chr_S = amp_S - 1
    chr_E = amp_E + 1

    # Get TN sides that S and E junctions connect to
    tn_side_S, tn_side_E = tnjc2.get_sides_of_chosen_tn()
    
    # Compute TN boundaries from junction positions and distances
    tn_scaf = amp_scaf  # TN is on same scaffold as amplicon
    if tn_side_S.side == Side.RIGHT:
        tn_left = amp_S - tn_side_S.distance
        tn_right = amp_E - tn_side_E.distance
    else:  # Side.LEFT
        tn_right = amp_S - tn_side_S.distance
        tn_left = amp_E - tn_side_E.distance

    # Define junction arms (chromosome/amplicon)
    chr_S_left = JcArm(scaf=amp_scaf, start=chr_S, dir=Orientation.REVERSE, flank=jc_width)
    amp_S_right = JcArm(scaf=amp_scaf, start=amp_S, dir=Orientation.FORWARD, flank=jc_width)
    amp_E_left = JcArm(scaf=amp_scaf, start=amp_E, dir=Orientation.REVERSE, flank=jc_width)
    chr_E_right = JcArm(scaf=amp_scaf, start=chr_E, dir=Orientation.FORWARD, flank=jc_width)

    # Define TN arms based on which sides connect to S/E junctions
    # Side.RIGHT (=1) → tn_left with FORWARD, Side.LEFT (=-1) → tn_right with REVERSE
    tn_S_arm = JcArm(
        scaf=tn_scaf, 
        start=tn_left if tn_side_S.side == Side.RIGHT else tn_right,
        dir=Orientation(tn_side_S.side.value),  # RIGHT=1=FORWARD, LEFT=-1=REVERSE
        flank=jc_width
    )
    tn_E_arm = JcArm(
        scaf=tn_scaf,
        start=tn_left if tn_side_E.side == Side.RIGHT else tn_right,
        dir=Orientation(tn_side_E.side.value),  # RIGHT=1=FORWARD, LEFT=-1=REVERSE
        flank=jc_width
    )

    # Create Junction objects for each junction type
    return {
        JunctionType.LEFT_REF:       Junction.from_jc_arms(chr_S_left, amp_S_right),
        JunctionType.LEFT_IS_TRANS:  Junction.from_jc_arms(chr_S_left, tn_S_arm),
        JunctionType.LEFT_MID_IS:    Junction.from_jc_arms(amp_E_left, tn_S_arm),
        JunctionType.LOST_IS:        Junction.from_jc_arms(amp_E_left, amp_S_right),
        JunctionType.RIGHT_MID_IS:   Junction.from_jc_arms(tn_E_arm,   amp_S_right),
        JunctionType.RIGHT_IS_TRANS: Junction.from_jc_arms(tn_E_arm,   chr_E_right),
        JunctionType.RIGHT_REF:      Junction.from_jc_arms(amp_E_left, chr_E_right),
    }


def write_junctions_fasta(
    junctions: dict[JunctionType, str],
    output_path: Path,
) -> None:
    """Write junction sequences to FASTA file.

    Args:
        junctions: Dict mapping JunctionType to sequence
        output_path: Output FASTA file path
    """
    ensure_parent_dir(output_path)

    with open(output_path, 'w') as f:
        for jtype, seq in sorted(junctions.items(), key=lambda x: x[0].value):
            f.write(f">{jtype.value}\n{seq}\n")


class CreateSyntheticJunctionsStep(Step[RecordTypedDf[FilteredTnJc2]]):
    """
    Creates 7 junction sequences per candidate for read alignment analysis.
    """

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        genome: Genome,
        output_dir: Path,
        read_length: int = 150,
        force: Optional[bool] = None,
    ):
        self.filtered_tnjc2s = filtered_tnjc2s
        self.genome = genome
        self.output_dir = Path(output_dir)
        self.read_length = read_length

        # Output files are per-candidate analysis directories
        # Only track outputs for candidates with valid chosen_tn_id
        self.analysis_dirs = []
        output_files = []
        for filtered_tnjc2 in filtered_tnjc2s:
            if filtered_tnjc2.chosen_tn_id is not None:
                analysis_dir = output_dir / "junctions" / filtered_tnjc2.analysis_dir
                self.analysis_dirs.append(analysis_dir)
                output_files.append(analysis_dir / "junctions.fasta")

        super().__init__(
            input_files=[genome.fasta_path],
            output_files=output_files,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[FilteredTnJc2]:
        """Create synthetic junctions for each candidate."""
        for filtered_tnjc2 in self.filtered_tnjc2s:
            # Skip candidates without valid chosen TN
            if filtered_tnjc2.chosen_tn_id is None:
                continue

            analysis_dir = self.output_dir / "junctions" / filtered_tnjc2.analysis_dir

            # Create junctions using new Genome/Scaffold methods
            junctions = create_synthetic_junctions(
                tnjc2=filtered_tnjc2,
                jc_width=self.read_length * 2,
                genome=self.genome,
            )

            # Write FASTA
            write_junctions_fasta(junctions, analysis_dir / "junctions.fasta")

        return self.filtered_tnjc2s

    def _save_output(self, output: RecordTypedDf[FilteredTnJc2]) -> None:
        """Output already saved in _calculate_output."""
        pass

    def _output_labels(self) -> list[str]:
        """Summarize outputs as count."""
        if self.output_files:
            n = len(self.filtered_tnjc2s)
            return [f"{n} junctions"]
        return []

    def load_outputs(self) -> RecordTypedDf[FilteredTnJc2]:
        """Return candidates (junction files are side effects)."""
        return self.filtered_tnjc2s

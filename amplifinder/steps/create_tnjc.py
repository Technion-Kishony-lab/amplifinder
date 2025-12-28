"""Step: Match junctions to TN elements (TnJc)."""

from pathlib import Path
from typing import Optional, List

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDf, Junction, SeqRefTnSide, RefTnSide, TnJunction, Orientation
from amplifinder.data_types.genome import Genome


class CreateTnJcStep(RecordTypedDfStep[TnJunction]):
    """Match junctions to TN elements by sequence comparison.

    For each junction, extracts flanking sequences and compares them to
    precomputed TN end sequences. Junctions matching a TN end are marked
    as TN-associated (TnJc).

    Based on assign_potential_ISs.m
    """

    def __init__(
        self,
        junctions: RecordTypedDf[Junction],
        seq_ref_tn_sides: RecordTypedDf[SeqRefTnSide],
        genome: Genome,
        output_dir: Path,
        max_dist_to_tn: int,
        trim_jc_flanking: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            jc_df: All junctions (breseq + reference TN junctions combined)
            ref_tn_end_seqs: Precomputed TN end sequences (from CreateRefTnEndSeqsStep)
            genome: Reference genome (for junction sequence extraction)
            output_dir: Directory to write output
            max_dist_to_tn: Maximum distance from junction to TN boundary
            trim_jc_flanking: Trim junction edges to avoid misalignment
            force: Force re-run
        """
        self.junctions = junctions
        self.seq_ref_tn_sides = seq_ref_tn_sides
        self.genome = genome
        self.max_dist_to_tn = max_dist_to_tn
        self.trim_jc_flanking = trim_jc_flanking

        super().__init__(
            output_dir=output_dir,
            input_files=[p for p in [genome.genbank_path, genome.fasta_path] if p],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[TnJunction]:
        """Match junctions to TN elements."""

        self.ref_seqs = self.genome.scaffold_sequences
        tnjc_records = []

        for jc in self.junctions:
            # Match against TN sequences
            matches1 = self._find_ref_tn_sides_matches(jc, 1)
            matches2 = self._find_ref_tn_sides_matches(jc, 2)

            is_arm1_tn = len(matches1) > 0
            is_arm2_tn = len(matches2) > 0

            # Skip if neither arm matches a TN
            if not is_arm1_tn and not is_arm2_tn:
                continue

            # Normalize: arm 1 should be the TN side
            switched = not is_arm1_tn and is_arm2_tn
            if switched:
                matches = matches2
                jc = jc.switch_sides()
            else:
                matches = matches1

            tnjc_records.append(TnJunction.from_other(jc, ref_tn_sides=matches, switched=switched))

        tnjcs = RecordTypedDf.from_records(tnjc_records, TnJunction)
        tnjcs = tnjcs.pipe(lambda df: df.sort_values(["scaf2", "pos2"]))

        info(f"Found {len(tnjcs)} TN-associated junctions (TnJc)")
        return tnjcs

    def _get_junction_arm_seq(self, jc: Junction, arm: int) -> str:
        """Extract sequence at junction arm."""
        seq = self.genome.get_junction_arm_sequence(jc, arm)
        if self.trim_jc_flanking == 0:
            return seq
        return seq[:-self.trim_jc_flanking]

    def _find_ref_tn_sides_matches(self, jc: Junction, arm: int) -> List[RefTnSide]:
        """Find TN elements matching a junction arm sequence."""
        jc_arm_seq = self._get_junction_arm_seq(jc, arm)
        matches = []

        for tn_side in self.seq_ref_tn_sides:
            # Check inward sequence (towards TN)
            pos = tn_side.seq_inward.find(jc_arm_seq)
            if pos >= 0:
                distance = pos + tn_side.offset
                if distance <= self.max_dist_to_tn:
                    matches.append(RefTnSide.from_other(tn_side, distance=distance))
                    continue

            # TODO: Do we need to also check outward sequence (away from TN)?

        return matches

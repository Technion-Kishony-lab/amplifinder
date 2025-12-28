"""Step: Match junctions to TN elements (TnJc)."""

from pathlib import Path
from typing import Optional, List

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDf, Junction, RefTnSide, TnJunction, Orientation, RefTnJunction
from Bio.Seq import reverse_complement
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
        ref_tnjcs: RecordTypedDf[RefTnJunction],
        genome: Genome,
        output_dir: Path,
        max_dist_to_tn: int,
        trim_jc_flanking: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            junctions: All junctions (breseq + reference TN junctions combined)
            ref_tnjcs: Reference TN junctions (used to get sequences on-the-fly)
            genome: Reference genome (for junction sequence extraction)
            output_dir: Directory to write output
            max_dist_to_tn: Maximum distance from junction to TN boundary
            trim_jc_flanking: Trim junction edges to avoid misalignment
            force: Force re-run
        """
        self.junctions = junctions
        self.ref_tnjcs = ref_tnjcs
        self.genome = genome
        self.max_dist_to_tn = max_dist_to_tn
        self.trim_jc_flanking = trim_jc_flanking
        self._ref_tn_seqs: Optional[dict[RefTnJunction, tuple[str, int]]] = None

        super().__init__(
            output_dir=output_dir,
            input_files=[p for p in [genome.genbank_path, genome.fasta_path] if p],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[TnJunction]:
        """Match junctions to TN elements."""

        self._precompute_ref_tnjcs_sequences()
        
        tnjcs = []

        for jc in self.junctions:
            # Match against TN sequences
            matches1 = self._find_ref_tn_sides_matches(jc, arm=1)
            matches2 = self._find_ref_tn_sides_matches(jc, arm=2)

            is_arm1_tn = len(matches1) > 0
            is_arm2_tn = len(matches2) > 0

            # Skip if neither arm matches a TN
            if not is_arm1_tn and not is_arm2_tn:
                continue

            # Normalize: arm 1 should be the TN side
            swapped = not is_arm1_tn and is_arm2_tn
            if swapped:
                matches = matches2
                jc = jc.swap_sides()
            else:
                matches = matches1

            tnjcs.append(TnJunction.from_other(jc, ref_tn_sides=matches, swapped=swapped))

        tnjcs = RecordTypedDf.from_records(tnjcs, TnJunction)
        tnjcs = tnjcs.pipe(lambda df: df.sort_values(["scaf2", "pos2"]))

        info(f"Found {len(tnjcs)} TN-associated junctions (TnJc)")
        return tnjcs

    def _precompute_ref_tnjcs_sequences(self) -> None:
        """Pre-compute sequences for all reference TN junctions (cache to avoid recomputing)."""
        self._ref_tn_seqs = {}
        for ref_jc in self.ref_tnjcs:
            seq_inward = self.genome.get_junction_sequence_arm2_to_arm1(ref_jc)
            offset = -ref_jc.flanking_right
            self._ref_tn_seqs[ref_jc] = (seq_inward, offset)

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

        for ref_jc in self.ref_tnjcs:
            # Use pre-computed sequence (cached in _calculate_output)
            seq_inward, offset = self._ref_tn_seqs[ref_jc]

            # Check inward sequence (towards TN)
            pos = seq_inward.find(jc_arm_seq)
            if pos >= 0:
                distance = pos + offset
                if distance <= self.max_dist_to_tn:
                    matches.append(RefTnSide.from_other(ref_jc.ref_tn_side, distance=distance))
                    continue
        return matches

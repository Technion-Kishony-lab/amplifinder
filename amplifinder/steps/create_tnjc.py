"""Step: Match junctions to TN elements (TnJc)."""

from pathlib import Path
from typing import Optional, List, Union

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.data_types import RecordTypedDf, Junction, RefTnSide, OffsetRefTnSide, TnJunction, RefTnJunction, BreseqJunction
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
        junctions: List[Union[BreseqJunction, RefTnJunction]],
        ref_tnjcs: RecordTypedDf[RefTnJunction],
        genome: Genome,
        output_dir: Path,
        max_dist_to_tn: int,
        trim_jc_flanking: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            junctions: All junctions (breseq + reference TN junctions combined) as a list to preserve types
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
        self._ref_tn_seqs: Optional[list[tuple[RefTnJunction, str, int]]] = None

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

            # If this is a reference TN junction, assert we get a distance==0 match on arm 1 (TN side)
            if isinstance(jc, RefTnJunction):
                ref_tn_side_matches = [m for m in matches1 if m.distance == 0 and m.is_same_side(jc.ref_tn_side)]
                assert len(ref_tn_side_matches) == 1, f"Reference TN junction {jc} should match itself with distance==0"

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
        # Sort but drop the old positional index so CSV won't write an Unnamed column
        tnjcs = tnjcs.pipe(lambda df: df.sort_values(["scaf2", "pos2"]).reset_index(drop=True))

        return tnjcs

    def _precompute_ref_tnjcs_sequences(self) -> None:
        """Pre-compute sequences for all reference TN junctions (cache to avoid recomputing)."""
        self._ref_tn_seqs = []
        for ref_tnjc in self.ref_tnjcs:
            seq_inward = self.genome.get_junction_sequence_arm2_to_arm1(ref_tnjc)
            self._ref_tn_seqs.append((ref_tnjc, seq_inward))

    def _get_junction_arm_seq(self, jc: Junction, arm: int) -> str:
        """Extract sequence at junction arm."""
        seq = self.genome.get_junction_arm_sequence(jc, arm)
        if not self.trim_jc_flanking:
            return seq
        return seq[:-self.trim_jc_flanking]

    def _find_ref_tn_sides_matches(self, jc: Junction, arm: int) -> List[OffsetRefTnSide]:
        """Find TN elements matching a junction arm sequence."""
        jc_arm_seq = self._get_junction_arm_seq(jc, arm)
        ref_tn_sides_matches = []

        # Use pre-computed ref-tn-side inward sequences (cached in _precompute_ref_tnjcs_sequences)
        for ref_tnjc, seq_inward in self._ref_tn_seqs:
            pos = seq_inward.find(jc_arm_seq)
            if pos >= 0:
                distance = pos - ref_tnjc.flanking_right - 1
                if abs(distance) <= self.max_dist_to_tn:
                    ref_tn_sides_matches.append(OffsetRefTnSide.from_other(ref_tnjc.ref_tn_side, distance=distance))
        return ref_tn_sides_matches

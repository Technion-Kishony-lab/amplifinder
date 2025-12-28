"""Steps: Create reference TN data (junctions and end sequences)."""

from pathlib import Path
from typing import Optional

from Bio.Seq import reverse_complement

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDf, RefTnLoc, SeqRefTnSide, RefTnJunction, RefTnSide, Side, Orientation
from amplifinder.data_types.genome import Genome


class CreateRefTnJcStep(RecordTypedDfStep[RefTnJunction]):
    """Create synthetic junction records (ref_tn_jc) from reference TN locations.

    For each TN element, creates two junction records representing the
    left and right TN-chromosome boundaries. These are used to detect
    reference TN elements that may not show up as novel junctions in breseq.

    Based on create_JC_of_reference_IS.m
    """

    def __init__(
        self,
        ref_tn_locs: RecordTypedDf[RefTnLoc],
        genome: Genome,
        output_dir: Path,
        source: str,
        reference_tn_out_span: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            ref_tn_locs: RecordDF with reference TN locations
            genome: Reference genome
            output_dir: Directory to write output
            source: Source name (genbank/isfinder) for filename prefix
            reference_tn_out_span: bp outside TN for unique chromosome seq
            force: Force re-run even if output exists
        """
        self.ref_tn_locs = ref_tn_locs
        self.genome = genome
        self.reference_tn_out_span = reference_tn_out_span
        output_file = output_dir / f"{source}_ref_tn_jc.csv"
        super().__init__(output_file=output_file, force=force)

    def _calculate_output(self) -> RecordTypedDf[RefTnJunction]:
        """Create synthetic junction records from TN locations."""

        jc_records = []

        for tn in self.ref_tn_locs:
            # pos1 and pos2 are stored as 1-based (genomic coordinates)
            # Left junction: TN left boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.tn_scaf, pos1=tn.loc_left, dir1=Orientation.FORWARD,
                scaf2=tn.tn_scaf, pos2=tn.loc_left - 1, dir2=Orientation.REVERSE,
                flanking_left=tn.length, flanking_right=self.reference_tn_out_span,
                ref_tn_side=RefTnSide(tn_id=tn.tn_id, side=Side.LEFT),
            ))

            # Right junction: TN right boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.tn_scaf, pos1=tn.loc_right, dir1=Orientation.REVERSE,
                scaf2=tn.tn_scaf, pos2=tn.loc_right + 1, dir2=Orientation.FORWARD,
                flanking_left=tn.length, flanking_right=self.reference_tn_out_span,
                ref_tn_side=RefTnSide(tn_id=tn.tn_id, side=Side.RIGHT),
            ))

        result = RecordTypedDf.from_records(jc_records, RefTnJunction)
        info(f"Created {len(jc_records)} reference TN junctions")
        return result


class CreateRefTnEndSeqsStep(RecordTypedDfStep[SeqRefTnSide]):
    """Extract full TN + flanks per junction side.

    For each `RefTnJunction` (side), returns the entire TN sequence with
    `max_dist_to_tn` margin on both sides (IS_seqs_with_margins in MATLAB),
    storing fwd and reverse-complement.

    Based on create_IS_end_seqs.m - creates IS_seqs_with_margins
    """

    def __init__(
        self,
        ref_tn_jcs: RecordTypedDf[RefTnJunction],
        genome: Genome,
        output_dir: Path,
        source: str,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            ref_tn_jcs: RecordDF with reference TN junctions (used to get tn_id and side)
            genome: Reference genome
            output_dir: Directory to write output
            source: Source name (genbank/isfinder) for filename prefix
            force: Force re-run
        """
        self.ref_tn_jcs = ref_tn_jcs
        self.genome = genome
        output_file = output_dir / f"{source}_tn_end_seqs.csv"
        super().__init__(output_file=output_file, force=force)

    def _calculate_output(self) -> RecordTypedDf[SeqRefTnSide]:
        """Extract TN sequences with margins (matching MATLAB IS_seqs_with_margins)."""
        seq_ref_tn_sides = [SeqRefTnSide.from_other(
            jc.ref_tn_side,
            offset=-jc.flanking_right,
            seq_inward=reverse_complement(self.genome.get_junction_sequence(jc)),
            seq_inward_rc=self.genome.get_junction_sequence(jc),
        ) for jc in self.ref_tn_jcs]

        seq_ref_tn_sides = RecordTypedDf.from_records(seq_ref_tn_sides, SeqRefTnSide)
        info(f"Created {len(seq_ref_tn_sides)} TN sequences with margins")
        return seq_ref_tn_sides

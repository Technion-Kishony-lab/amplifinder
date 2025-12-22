"""Steps: Create reference TN data (junctions and end sequences)."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDF, TnLoc, SeqRefTnSide, RefTnJunction, RefTnSide, Side, Orientation
from amplifinder.data_types.genome import Genome
from amplifinder.utils.fasta import reverse_complement


class CreateRefTnJcStep(Step[RecordTypedDF[RefTnJunction]]):
    """Create synthetic junction records (ref_tn_jc) from reference TN locations.

    For each TN element, creates two junction records representing the
    left and right TN-chromosome boundaries. These are used to detect
    reference TN elements that may not show up as novel junctions in breseq.

    Based on create_JC_of_reference_IS.m
    """

    def __init__(
        self,
        tn_loc: RecordTypedDF[TnLoc],
        genome: Genome,
        output_dir: Path,
        source: str,
        reference_tn_out_span: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tn_loc: RecordDF with TN locations
            genome: Reference genome
            output_dir: Directory to write output
            source: Source name (genbank/isfinder) for filename prefix
            reference_tn_out_span: bp outside TN for unique chromosome seq
            force: Force re-run even if output exists
        """
        self.tn_loc = tn_loc
        self.genome = genome
        self.output_dir = Path(output_dir)
        self.reference_tn_out_span = reference_tn_out_span

        self.output_file = self.output_dir / f"{source}_ref_tn_jc.csv"

        super().__init__(
            input_files=[],  # tn_loc passed in memory
            output_files=[self.output_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDF[RefTnJunction]:
        """Create synthetic junction records from TN locations."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        jc_records = []

        for tn in self.tn_loc:
            tn_length = tn.LocRight - tn.LocLeft + 1

            # Left junction: TN left boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.TN_scaf, pos1=tn.LocLeft, dir1=Orientation.FORWARD,
                scaf2=tn.TN_scaf, pos2=tn.LocLeft - 1, dir2=Orientation.REVERSE,
                flanking_left=tn_length, flanking_right=self.reference_tn_out_span,
                ref_tn_side=RefTnSide(tn_id=tn.ID, side=Side.LEFT),
            ))

            # Right junction: TN right boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.TN_scaf, pos1=tn.LocRight, dir1=Orientation.REVERSE,
                scaf2=tn.TN_scaf, pos2=tn.LocRight + 1, dir2=Orientation.FORWARD,
                flanking_left=self.reference_tn_out_span, flanking_right=tn_length,
                ref_tn_side=RefTnSide(tn_id=tn.ID, side=Side.RIGHT),
            ))

        result = RecordTypedDF.from_records(jc_records, RefTnJunction)
        info(f"Created {len(jc_records)} reference TN junctions")
        return result

    def _save_output(self, output: RecordTypedDF[RefTnJunction]) -> None:
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[RefTnJunction]:
        """Load reference TN junctions from output file."""
        return RecordTypedDF.from_csv(self.output_file, RefTnJunction)


class CreateRefTnEndSeqsStep(Step[RecordTypedDF[SeqRefTnSide]]):
    """Extract TN sequences with margins for junction matching.

    For each TN element, creates full TN sequence with margins (IS_seqs_with_margins),
    matching MATLAB behavior. This is the full TN element plus margins on both sides,
    not just the end sequences.

    Based on create_IS_end_seqs.m - creates IS_seqs_with_margins
    """

    def __init__(
        self,
        ref_tn_jc: RecordTypedDF[RefTnJunction],
        tn_loc: RecordTypedDF[TnLoc],
        genome: Genome,
        output_dir: Path,
        source: str,
        max_dist_to_tn: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            ref_tn_jc: RecordDF with reference TN junctions (used to get tn_id and side)
            tn_loc: RecordDF with TN locations (used to get full TN boundaries)
            genome: Reference genome
            output_dir: Directory to write output
            source: Source name (genbank/isfinder) for filename prefix
            max_dist_to_tn: Margin around TN (out_span in MATLAB)
            force: Force re-run
        """
        self.ref_tn_jc = ref_tn_jc
        self.tn_loc = tn_loc
        self.genome = genome
        self.output_dir = Path(output_dir)
        self.max_dist_to_tn = max_dist_to_tn

        self.output_file = self.output_dir / f"{source}_tn_end_seqs.csv"

        super().__init__(
            input_files=[],
            output_files=[self.output_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDF[SeqRefTnSide]:
        """Extract TN sequences with margins (matching MATLAB IS_seqs_with_margins)."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        ref_seqs = self.genome.sequences
        out_span = self.max_dist_to_tn  # MATLAB: out_span = max_dist_to_IS
        
        # Create mapping from tn_id to TN location
        tn_loc_map = {tn.ID: tn for tn in self.tn_loc}
        
        records = []

        for jc in self.ref_tn_jc:
            if jc.scaf1 not in ref_seqs:
                raise ValueError(f"TN scaffold {jc.scaf1} not in genome")
            
            # Get TN location for this junction
            tn_id = jc.ref_tn_side.tn_id
            if tn_id not in tn_loc_map:
                continue  # Skip if TN not found
            
            tn = tn_loc_map[tn_id]
            seq = ref_seqs[jc.scaf1]
            
            # MATLAB: IS_seqs_with_margins{i,1} = seq(l-out_span:r+out_span)
            # Where l=LocLeft, r=LocRight (1-based inclusive in MATLAB)
            # Python slicing: seq[l-1:r] where l-1 is 0-based start, r is 0-based exclusive end
            # But LocRight is 1-based inclusive, so for Python we need LocRight (exclusive end)
            l = tn.LocLeft  # 1-based inclusive
            r = tn.LocRight  # 1-based inclusive
            
            # Create full TN sequence with margins (like MATLAB IS_seqs_with_margins)
            # MATLAB: seq(l-out_span:r+out_span) - 1-based inclusive
            # Python: seq[(l-out_span-1):(r+out_span)] - 0-based, exclusive end
            start = max(0, (l - out_span - 1))  # Convert to 0-based
            end = min(len(seq), r + out_span)  # r is 1-based inclusive, so r+out_span is exclusive end
            seq_with_margins = seq[start:end]
            
            # MATLAB: IS_seqs_with_margins{i,2} = seqrcomplement(seq(l-out_span:r+out_span))
            seq_with_margins_rc = reverse_complement(seq_with_margins)

            records.append(SeqRefTnSide.from_other(
                jc.ref_tn_side,
                seq_fwd=seq_with_margins,
                seq_rc=seq_with_margins_rc,
            ))

        result = RecordTypedDF.from_records(records, SeqRefTnSide)
        info(f"Created {len(records)} TN sequences with margins")
        return result

    def _save_output(self, output: RecordTypedDF[SeqRefTnSide]) -> None:
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[SeqRefTnSide]:
        """Load TN end sequences from output file."""
        return RecordTypedDF.from_csv(self.output_file, SeqRefTnSide)

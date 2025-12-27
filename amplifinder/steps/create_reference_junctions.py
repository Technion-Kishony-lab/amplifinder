"""Steps: Create reference TN data (junctions and end sequences)."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDf, RefTnLoc, SeqRefTnSide, RefTnJunction, RefTnSide, Side, Orientation
from amplifinder.data_types.genome import Genome
from amplifinder.utils.fasta import reverse_complement


class CreateRefTnJcStep(RecordTypedDfStep[RefTnJunction]):
    """Create synthetic junction records (ref_tn_jc) from reference TN locations.

    For each TN element, creates two junction records representing the
    left and right TN-chromosome boundaries. These are used to detect
    reference TN elements that may not show up as novel junctions in breseq.

    Based on create_JC_of_reference_IS.m
    """

    def __init__(
        self,
        tn_loc: RecordTypedDf[RefTnLoc],
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
        self.reference_tn_out_span = reference_tn_out_span
        
        output_file = output_dir / f"{source}_ref_tn_jc.csv"
        
        super().__init__(
            output_file=output_file,
            input_files=[],  # tn_loc passed in memory
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[RefTnJunction]:
        """Create synthetic junction records from TN locations."""

        jc_records = []

        for tn in self.tn_loc:
            tn_length = tn.loc_right - tn.loc_left + 1

            # Left junction: TN left boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.tn_scaf, pos1=tn.loc_left, dir1=Orientation.FORWARD,
                scaf2=tn.tn_scaf, pos2=tn.loc_left - 1, dir2=Orientation.REVERSE,
                flanking_left=tn_length, flanking_right=self.reference_tn_out_span,
                ref_tn_side=RefTnSide(tn_id=tn.tn_id, side=Side.LEFT),
            ))

            # Right junction: TN right boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.tn_scaf, pos1=tn.loc_right, dir1=Orientation.REVERSE,
                scaf2=tn.tn_scaf, pos2=tn.loc_right + 1, dir2=Orientation.FORWARD,
                flanking_left=self.reference_tn_out_span, flanking_right=tn_length,
                ref_tn_side=RefTnSide(tn_id=tn.tn_id, side=Side.RIGHT),
            ))

        result = RecordTypedDf.from_records(jc_records, RefTnJunction)
        info(f"Created {len(jc_records)} reference TN junctions")
        return result


class CreateRefTnEndSeqsStep(RecordTypedDfStep[SeqRefTnSide]):
    """Extract TN sequences with margins for junction matching.

    For each TN element, creates full TN sequence with margins (IS_seqs_with_margins),
    matching MATLAB behavior. This is the full TN element plus margins on both sides,
    not just the end sequences.

    Based on create_IS_end_seqs.m - creates IS_seqs_with_margins
    """

    def __init__(
        self,
        ref_tn_jc: RecordTypedDf[RefTnJunction],
        tn_loc: RecordTypedDf[RefTnLoc],
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
        self.max_dist_to_tn = max_dist_to_tn
        
        output_file = output_dir / f"{source}_tn_end_seqs.csv"
        
        super().__init__(
            output_file=output_file,
            input_files=[],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[SeqRefTnSide]:
        """Extract TN sequences with margins (matching MATLAB IS_seqs_with_margins)."""

        ref_seqs = self.genome.sequences
        out_span = self.max_dist_to_tn  # MATLAB: out_span = max_dist_to_IS
        
        # Create mapping from tn_id to TN location
        tn_loc_map = {tn.tn_id: tn for tn in self.tn_loc}
        
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
            l = tn.loc_left  # 1-based inclusive
            r = tn.loc_right  # 1-based inclusive
            
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

        result = RecordTypedDf.from_records(records, SeqRefTnSide)
        info(f"Created {len(records)} TN sequences with margins")
        return result

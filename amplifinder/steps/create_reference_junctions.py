"""Steps: Create reference TN data (junctions and end sequences)."""

from pathlib import Path
from typing import Optional

from Bio.Seq import Seq

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDF, TnLoc, TnEndSeq, RefTnJunction, Side, Orientation
from amplifinder.data_types.genome import Genome


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
                refTN=tn.ID, tn_side=Side.LEFT,
            ))

            # Right junction: TN right boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=tn.TN_scaf, pos1=tn.LocRight, dir1=Orientation.REVERSE,
                scaf2=tn.TN_scaf, pos2=tn.LocRight + 1, dir2=Orientation.FORWARD,
                flanking_left=self.reference_tn_out_span, flanking_right=tn_length,
                refTN=tn.ID, tn_side=Side.RIGHT,
            ))

        result = RecordTypedDF.from_records(jc_records, RefTnJunction)
        info(f"Created {len(jc_records)} reference TN junctions")
        return result

    def _save_output(self, output: RecordTypedDF[RefTnJunction]) -> None:
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[RefTnJunction]:
        """Load reference TN junctions from output file."""
        return RecordTypedDF.from_csv(self.output_file, RefTnJunction)


class CreateRefTnEndSeqsStep(Step[RecordTypedDF[TnEndSeq]]):
    """Extract TN end sequences for junction matching.

    For each reference TN junction, extracts the sequence around the TN boundary
    plus its reverse complement. Used by matching step.

    Based on create_IS_end_seqs.m
    """

    def __init__(
        self,
        ref_tn_jc: RecordTypedDF[RefTnJunction],
        genome: Genome,
        output_dir: Path,
        source: str,
        max_dist_to_tn: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            ref_tn_jc: RecordDF with reference TN junctions
            genome: Reference genome
            output_dir: Directory to write output
            source: Source name (genbank/isfinder) for filename prefix
            max_dist_to_tn: Margin around TN boundary
            force: Force re-run
        """
        self.ref_tn_jc = ref_tn_jc
        self.genome = genome
        self.output_dir = Path(output_dir)
        self.max_dist_to_tn = max_dist_to_tn

        self.output_file = self.output_dir / f"{source}_tn_end_seqs.csv"

        super().__init__(
            input_files=[],
            output_files=[self.output_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDF[TnEndSeq]:
        """Extract TN end sequences."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        ref_seqs = self.genome.sequences
        margin = self.max_dist_to_tn
        records = []

        for jc in self.ref_tn_jc:
            if jc.scaf1 not in ref_seqs:
                raise ValueError(f"TN scaffold {jc.scaf1} not in genome")

            seq = ref_seqs[jc.scaf1]
            pos = jc.pos1 - 1  # 0-based
            start = max(0, pos - margin)
            end = min(len(seq), pos + margin + 1)
            seq_fwd = seq[start:end]

            records.append(TnEndSeq(
                tn_id=jc.refTN,
                tn_side=jc.tn_side,
                seq_fwd=seq_fwd,
                seq_rc=str(Seq(seq_fwd).reverse_complement()),
            ))

        result = RecordTypedDF.from_records(records, TnEndSeq)
        info(f"Created {len(records)} TN end sequences")
        return result

    def _save_output(self, output: RecordTypedDF[TnEndSeq]) -> None:
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[TnEndSeq]:
        """Load TN end sequences from output file."""
        return RecordTypedDF.from_csv(self.output_file, TnEndSeq)

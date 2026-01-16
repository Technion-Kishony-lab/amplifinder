"""Steps: Create reference TN data (junctions)."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.data_types import RecordTypedDf, RefTn, RefTnJunction
from amplifinder.data_types.genome import Genome


class CreateRefTnJcStep(RecordTypedDfStep[RefTnJunction]):
    """Create synthetic junction records (ref_tnjc) from reference TN locations.

    For each TN element, creates two junction records representing the
    left and right TN-chromosome boundaries. These are used to detect
    reference TN elements that may not show up as novel junctions in breseq.

    Based on create_JC_of_reference_IS.m
    """

    def __init__(
        self,
        ref_tn_locs: RecordTypedDf[RefTn],
        genome: Genome,
        output_dir: Path,
        source: str,
        reference_IS_out_span: int,
        reference_IS_in_span: Optional[int] = None,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            ref_tn_locs: RecordDF with reference TN locations
            genome: Reference genome
            output_dir: Directory to write output
            source: Source name (genbank/isfinder) for filename prefix
            reference_IS_out_span: bp outside IS for unique chromosome seq
            reference_IS_in_span: bp into IS element for junction flanking sequence
                                  None means use the entire IS element (default)
            force: Force re-run even if output exists
        """
        self.ref_tn_locs = ref_tn_locs
        self.genome = genome
        self.reference_IS_out_span = reference_IS_out_span
        self.reference_IS_in_span = reference_IS_in_span
        output_file = output_dir / f"{source}_ref_tnjc.csv"
        super().__init__(output_file=output_file, force=force)

    def _calculate_output(self) -> RecordTypedDf[RefTnJunction]:
        """Create synthetic junction records from TN locations."""
        ref_tnjcs = []
        for tn in self.ref_tn_locs:
            left_jc, right_jc = tn.get_junctions(
                out_flanks=self.reference_IS_out_span,
                in_flanks=self.reference_IS_in_span,
            )
            ref_tnjcs.append(left_jc)
            ref_tnjcs.append(right_jc)

        ref_tnjcs = RecordTypedDf.from_records(ref_tnjcs, RefTnJunction)
        return ref_tnjcs

"""Step: Create synthetic junctions from reference TN locations."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types.records import RefTnJunction, REF_TN_JC_SCHEMA


class CreateReferenceJunctionsStep(Step[pd.DataFrame]):
    """Create synthetic junction records (ref_tn_jc) from reference TN locations.

    For each TN element, creates two junction records representing the
    left and right TN-chromosome boundaries. These are used to detect
    reference TN elements that may not show up as novel junctions in breseq.

    Based on create_JC_of_reference_IS.m
    """

    def __init__(
        self,
        tn_loc: pd.DataFrame,
        ref_name: str,
        output_dir: Path,
        reference_tn_out_span: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tn_loc: DataFrame with TN locations (ID, TN_Name, TN_scaf, etc.)
            ref_name: Reference genome name
            output_dir: Directory to write output
            reference_tn_out_span: bp outside TN for unique chromosome seq
            force: Force re-run even if output exists
        """
        self.tn_loc = tn_loc
        self.ref_name = ref_name
        self.output_dir = Path(output_dir)
        self.reference_tn_out_span = reference_tn_out_span

        self.output_file = self.output_dir / "ref_tn_jc.csv"

        super().__init__(
            inputs=[],  # tn_loc passed in memory
            outputs=[self.output_file],
            force=force,
        )

    def _run(self) -> None:
        """Create synthetic junction records from TN locations."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        jc_records = []

        for idx, tn_row in self.tn_loc.iterrows():
            tn_length = tn_row["LocRight"] - tn_row["LocLeft"] + 1
            scaf = tn_row["TN_scaf"]

            # Left junction: TN left boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=scaf, pos1=int(tn_row["LocLeft"]), dir1=1,
                scaf2=scaf, pos2=int(tn_row["LocLeft"]) - 1, dir2=-1,
                flanking_left=tn_length, flanking_right=self.reference_tn_out_span,
                refTN=int(idx), tn_side="left",
            ))

            # Right junction: TN right boundary -> chromosome
            jc_records.append(RefTnJunction(
                num=0,
                scaf1=scaf, pos1=int(tn_row["LocRight"]), dir1=-1,
                scaf2=scaf, pos2=int(tn_row["LocRight"]) + 1, dir2=1,
                flanking_left=self.reference_tn_out_span, flanking_right=tn_length,
                refTN=int(idx), tn_side="right",
            ))

        REF_TN_JC_SCHEMA.to_csv(jc_records, self.output_file)
        info(f"Created {len(jc_records)} reference TN junctions")

    def read_outputs(self) -> pd.DataFrame:
        """Load reference TN junctions from output file."""
        return REF_TN_JC_SCHEMA.read_csv(self.output_file)

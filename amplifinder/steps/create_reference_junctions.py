"""Step: Create synthetic junctions from reference TN locations."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.logger import info


class CreateReferenceJunctionsStep(Step):
    """Create synthetic junction records (ref_tn_jc) from reference TN locations.

    For each TN element, creates two junction records representing the
    left and right TN-chromosome boundaries. These are used to detect
    reference TN elements that may not show up as novel junctions in breseq.

    Based on create_JC_of_reference_IS.m
    """

    REFERENCE_TN_OUT_SPAN = 100  # bp outside TN for unique chromosome seq

    def __init__(
        self,
        tn_loc: pd.DataFrame,
        ref_name: str,
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tn_loc: DataFrame with TN locations (ID, TN_Name, TN_scaf, etc.)
            ref_name: Reference genome name
            output_dir: Directory to write output
            force: Force re-run even if output exists
        """
        self.tn_loc = tn_loc
        self.ref_name = ref_name
        self.output_dir = Path(output_dir)

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

            # Left junction: TN left boundary -> chromosome
            jc_records.append(self._make_jc_record(
                tn_row=tn_row,
                ref_tn_idx=idx,
                tn_side="left",
                pos1=tn_row["LocLeft"],
                dir1=1,
                pos2=tn_row["LocLeft"] - 1,
                dir2=-1,
                flanking_tn=tn_length,
            ))

            # Right junction: TN right boundary -> chromosome
            jc_records.append(self._make_jc_record(
                tn_row=tn_row,
                ref_tn_idx=idx,
                tn_side="right",
                pos1=tn_row["LocRight"],
                dir1=-1,
                pos2=tn_row["LocRight"] + 1,
                dir2=1,
                flanking_tn=tn_length,
            ))

        jc_reftn = pd.DataFrame(jc_records)
        jc_reftn.to_csv(self.output_file, index=False)
        info(f"Created {len(jc_reftn)} reference TN junctions")

    def _make_jc_record(
        self,
        tn_row: pd.Series,
        ref_tn_idx: int,
        tn_side: str,
        pos1: int,
        dir1: int,
        pos2: int,
        dir2: int,
        flanking_tn: int,
    ) -> dict:
        """Create a single JC record for a TN boundary.

        Following the structure from create_JC_of_reference_IS.m:
        - Side 1 (scaf1, pos1, dir1) = TN side
        - Side 2 (scaf2, pos2, dir2) = chromosome side
        """
        return {
            "num": 0,  # 0 indicates reference junction
            "scaf1": tn_row["TN_scaf"],
            "pos1": int(pos1),
            "dir1": int(dir1),
            "scaf2": tn_row["TN_scaf"],
            "pos2": int(pos2),
            "dir2": int(dir2),
            "refTN": int(ref_tn_idx),  # Link back to TN_loc index
            "tn_side": tn_side,
            "flanking_left": flanking_tn if tn_side == "left" else self.REFERENCE_TN_OUT_SPAN,
            "flanking_right": flanking_tn if tn_side == "right" else self.REFERENCE_TN_OUT_SPAN,
            "reject": None,
            "new_junction_read_count": 0,
            "new_junction_coverage": 0.0,
        }

    def read_outputs(self) -> pd.DataFrame:
        """Load reference TN junctions from output file."""
        return pd.read_csv(self.output_file)

    def run_and_read_outputs(self) -> pd.DataFrame:
        """Run step and return reference TN junctions."""
        return super().run_and_read_outputs()

"""Step: Find TN elements from GenBank annotations."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data.schemas import TN_LOC_SCHEMA
from amplifinder.utils.genbank import find_TN_elements


class LocateTNsUsingGenbank(Step):
    """Extract TN elements from GenBank file annotations using BioPython.

    Parses GenBank file for 'insertion sequence' features (based on findISinRef.m).
    """

    def __init__(
        self,
        genbank_path: Path,
        ref_name: str,
        ref_path: Path,
        force: Optional[bool] = None,
    ):
        self.genbank_path = Path(genbank_path)
        self.ref_name = ref_name
        self.ref_path = Path(ref_path)

        # Output paths
        self.genbank_dir = self.ref_path / "genbank"
        self.TN_loc_output = self.genbank_dir / f"{ref_name}_TN_loc.csv"

        super().__init__(
            inputs=[self.genbank_path],
            outputs=[self.TN_loc_output],
            force=force,
        )

    def _run(self) -> None:
        """Parse GenBank file and extract TN locations."""
        self.genbank_dir.mkdir(parents=True, exist_ok=True)

        # Use BioPython to parse features
        records = find_TN_elements(self.genbank_path, self.ref_name)

        TN_loc = TN_LOC_SCHEMA.from_records(records)
        TN_LOC_SCHEMA.to_csv(TN_loc, self.TN_loc_output)
        info(f"Found {len(TN_loc)} TN elements in GenBank annotations")

    def read_outputs(self) -> pd.DataFrame:
        """Load TN locations from output file."""
        return TN_LOC_SCHEMA.read_csv(self.TN_loc_output)

    def run_and_read_outputs(self) -> pd.DataFrame:
        """Run step and return TN locations."""
        return super().run_and_read_outputs()

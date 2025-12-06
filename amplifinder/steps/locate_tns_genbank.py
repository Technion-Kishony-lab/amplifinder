"""Step: Find TN elements from GenBank annotations."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types import TN_LOC_SCHEMA
from amplifinder.utils.genbank import find_tn_elements


class LocateTNsUsingGenbank(Step[Optional[pd.DataFrame]]):
    """Extract TN elements from GenBank file annotations using BioPython.

    Parses GenBank file for 'insertion sequence' features (based on findISinRef.m).
    Returns None if no GenBank file is provided.
    """

    def __init__(
        self,
        genbank_path: Optional[Path],
        ref_name: str,
        ref_path: Path,
        force: Optional[bool] = None,
    ):
        self.genbank_path = Path(genbank_path) if genbank_path else None
        self.ref_name = ref_name
        self.ref_path = Path(ref_path)

        # Output paths
        self.genbank_dir = self.ref_path / "genbank"
        self.tn_loc_output = self.genbank_dir / f"{ref_name}_tn_loc.csv"

        super().__init__(
            inputs=[self.genbank_path] if self.genbank_path else [],
            outputs=[self.tn_loc_output],
            force=force,
        )

    def _run(self) -> None:
        """Parse GenBank file and extract TN locations."""
        if self.genbank_path is None:
            info("No GenBank file provided - skipping GenBank TN annotation")
            return

        self.genbank_dir.mkdir(parents=True, exist_ok=True)
        records = find_tn_elements(self.genbank_path, self.ref_name)
        tn_loc = TN_LOC_SCHEMA.from_records(records)
        TN_LOC_SCHEMA.to_csv(tn_loc, self.tn_loc_output)
        info(f"Found {len(tn_loc)} TN elements in GenBank annotations")

    def read_outputs(self) -> Optional[pd.DataFrame]:
        """Load TN locations from output file, or None if no GenBank."""
        if self.genbank_path is None:
            return None
        return TN_LOC_SCHEMA.read_csv(self.tn_loc_output)

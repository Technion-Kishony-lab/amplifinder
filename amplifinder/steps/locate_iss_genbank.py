"""Step: Find IS elements from GenBank annotations."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types.schema import IS_LOC_SCHEMA
from amplifinder.utils.genbank import find_IS_elements


class LocateISsUsingGenbank(Step):
    """Extract IS elements from GenBank file annotations using BioPython.
    
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
        self.IS_loc_output = self.genbank_dir / f"{ref_name}_IS_loc.csv"
        
        super().__init__(
            inputs=[self.genbank_path],
            outputs=[self.IS_loc_output],
            force=force,
        )
    
    def _run(self) -> None:
        """Parse GenBank file and extract IS locations."""
        self.genbank_dir.mkdir(parents=True, exist_ok=True)
        
        # Use BioPython to parse features
        records = find_IS_elements(self.genbank_path, self.ref_name)
        
        IS_loc = IS_LOC_SCHEMA.from_records(records)
        IS_LOC_SCHEMA.to_csv(IS_loc, self.IS_loc_output)
        info(f"Found {len(IS_loc)} IS elements in GenBank annotations")
    
    def read_outputs(self) -> pd.DataFrame:
        """Load IS locations from output file."""
        return IS_LOC_SCHEMA.read_csv(self.IS_loc_output)
    
    def run_and_read_outputs(self) -> pd.DataFrame:
        """Run step and return IS locations."""
        return super().run_and_read_outputs()

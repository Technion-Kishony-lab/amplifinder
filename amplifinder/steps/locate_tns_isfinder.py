"""Step: Find TN elements by BLAST against ISfinder database."""

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.tools.blast import run_blastn, parse_blast_csv, make_blast_db
from amplifinder.utils.fasta import read_fasta_lengths
from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data.schemas import TN_LOC_SCHEMA


class LocateTNsUsingISfinder(Step):
    """BLAST reference genome against ISfinder database to find TN elements."""

    def __init__(
        self,
        ref_fasta: Path,
        ref_name: str,
        ref_path: Path,
        isdb_path: Path,
        evalue: float = 1e-4,
        critical_coverage: float = 0.9,
        force: Optional[bool] = None,
    ):
        self.ref_fasta = Path(ref_fasta)
        self.ref_name = ref_name
        self.ref_path = Path(ref_path)
        self.isdb_path = Path(isdb_path)
        self.evalue = evalue
        self.critical_coverage = critical_coverage

        # Output paths
        self.isfinder_dir = self.ref_path / "ISfinder"
        self.blast_output = self.isfinder_dir / f"{ref_name}_blast.txt"
        self.TN_loc_output = self.isfinder_dir / f"{ref_name}_TN_loc.csv"

        super().__init__(
            inputs=[self.ref_fasta, self.isdb_path],
            outputs=[self.TN_loc_output],
            force=force,
        )

    def _run(self) -> None:
        """Run BLAST and parse results."""
        self.isfinder_dir.mkdir(parents=True, exist_ok=True)

        # Create BLAST DB if needed
        if not self.isdb_path.with_suffix(".nhr").exists():
            make_blast_db(self.isdb_path, self.isdb_path)

        # Run BLAST
        run_blastn(
            query=self.ref_fasta,
            db=self.isdb_path,
            output=self.blast_output,
            evalue=self.evalue,
        )

        # Parse BLAST results and save
        TN_loc = self._parse_blast()
        TN_LOC_SCHEMA.to_csv(TN_loc, self.TN_loc_output)
        info(f"Found {len(TN_loc)} TN elements via ISfinder")

    def _parse_blast(self) -> pd.DataFrame:
        """Parse BLAST output and filter to TN locations."""
        df = parse_blast_csv(self.blast_output)
        if df.empty:
            return TN_LOC_SCHEMA.empty()

        # Load ISfinder DB to get TN lengths
        TN_lengths = read_fasta_lengths(self.isdb_path)

        # Convert to TN locations
        return self._blast_to_TN_loc(df, TN_lengths)

    def _blast_to_TN_loc(self, blast_df: pd.DataFrame, tn_lengths: dict) -> pd.DataFrame:
        """Convert BLAST results to TN location table (matlab blast2ISloc logic)."""
        df = blast_df.copy()

        # Remove duplicate alignments (keep first by qstart, qend)
        df = df.drop_duplicates(subset=["qstart", "qend"], keep="first")

        # Add subject length
        df["subject_length"] = df["subject"].map(tn_lengths).fillna(0).astype(int)

        # Sort by query, qstart, bitscore
        df = df.sort_values(["query", "qstart", "bitscore"])

        # Filter by coverage (alignment covers >90% of TN element)
        df["coverage"] = abs(df["sstart"] - df["send"]) / df["subject_length"]
        df = df[df["coverage"] > self.critical_coverage]

        # Build TN_loc table
        # Detect strand: sstart > send means reverse complement
        is_complement = (df["sstart"] > df["send"]).values

        # Ensure LocLeft < LocRight
        loc_left = df[["qstart", "qend"]].min(axis=1).values
        loc_right = df[["qstart", "qend"]].max(axis=1).values

        return TN_LOC_SCHEMA.create({
            "ID": range(1, len(df) + 1),
            "TN_Name": df["subject"].values,
            "TN_scaf": df["query"].values,
            "LocLeft": loc_left,
            "LocRight": loc_right,
            "Complement": is_complement,
            "Join": False,
        })

    def read_outputs(self) -> pd.DataFrame:
        """Load TN locations from output file."""
        return TN_LOC_SCHEMA.read_csv(self.TN_loc_output)

    def run_and_read_outputs(self) -> pd.DataFrame:
        """Run step and return TN locations."""
        return super().run_and_read_outputs()

"""Step: Find TN elements by BLAST against ISfinder database."""

import pandas as pd

from pathlib import Path
from typing import Optional

from abc import abstractmethod

from amplifinder.tools.blast import run_blastn, parse_blast_csv, make_blast_db
from amplifinder.utils.fasta import read_fasta_lengths
from amplifinder.utils.genbank import find_tn_elements
from amplifinder.utils.file_lock import locked_resource
from amplifinder.utils.tools import ensure_dir
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDf, RefTnLoc, Genome
from amplifinder.steps.base import Step


# Base class for TN location steps

class LocateTNsStep(Step[Optional[RecordTypedDf[RefTnLoc]]]):
    """Base class for steps that locate TN elements."""

    def __init__(
        self,
        genome: Genome,
        output_dir: Path,
        source: str,
        input_files: list[Path],
        force: Optional[bool] = None,
    ):
        self.genome = genome
        self.output_dir = Path(output_dir)

        self.output_file = self.output_dir / f"{source}_tn_loc.csv"

        super().__init__(
            input_files=input_files,
            output_files=[self.output_file],
            force=force,
        )

    def _save_output(self, output: Optional[RecordTypedDf[RefTnLoc]]) -> None:
        if output is not None:
            output.to_csv(self.output_file)

    def load_outputs(self) -> Optional[RecordTypedDf[RefTnLoc]]:
        """Load TN locations from output file."""
        if not self.output_file.exists():
            return None
        return RecordTypedDf.from_csv(self.output_file, RefTnLoc)

    @abstractmethod
    def _calculate_output(self) -> Optional[RecordTypedDf[RefTnLoc]]:
        """Run the TN location logic."""
        pass


# Locate TNs using GenBank annotations

class LocateTNsUsingGenbankStep(LocateTNsStep):
    """Extract TN elements from GenBank file annotations using BioPython.

    Parses GenBank file for 'insertion sequence' features (based on findISinRef.m).
    Returns None if no GenBank file is provided.
    """

    def __init__(
        self,
        genome: Genome,
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        super().__init__(
            genome=genome,
            output_dir=output_dir,
            source="genbank",
            input_files=[genome.genbank_path] if genome.genbank_path else [],
            force=force,
        )

    def _calculate_output(self) -> Optional[RecordTypedDf[RefTnLoc]]:
        """Parse GenBank file and extract TN locations."""
        if self.genome.genbank_path is None:
            info("No GenBank file provided - skipping GenBank TN annotation")
            return None

        ensure_dir(self.output_dir)
        records = find_tn_elements(self.genome.genbank_path, self.genome.name)
        tn_loc = RecordTypedDf.from_records(records, RefTnLoc)
        info(f"Found {len(tn_loc)} TN elements in GenBank annotations")
        return tn_loc


# Locate TNs using ISfinder database

class LocateTNsUsingISfinderStep(LocateTNsStep):
    """BLAST reference genome against ISfinder database to find TN elements."""

    def __init__(
        self,
        genome: Genome,
        output_dir: Path,
        isdb_path: Path,
        evalue: float,
        critical_coverage: float,
        force: Optional[bool] = None,
    ):
        self.isdb_path = Path(isdb_path)
        self.evalue = evalue
        self.critical_coverage = critical_coverage

        super().__init__(
            genome=genome,
            output_dir=output_dir,
            source="isfinder",
            input_files=[genome.fasta_path, self.isdb_path],
            force=force,
        )

        # Additional output
        self.blast_output = self.output_dir / "isfinder_blast.txt"
        self.output_files.append(self.blast_output)

    def _calculate_output(self) -> RecordTypedDf[RefTnLoc]:
        """Run BLAST and parse results."""
        ensure_dir(self.output_dir)

        # Create BLAST DB if needed (with lock to prevent parallel creation)
        with locked_resource(self.isdb_path, "blast_db", timeout=300):
            if not self.isdb_path.with_suffix(".nhr").exists():
                info("Creating ISfinder BLAST database...")
                make_blast_db(self.isdb_path, self.isdb_path)

        # Run BLAST
        run_blastn(
            query=self.genome.fasta_path,
            db=self.isdb_path,
            output=self.blast_output,
            evalue=self.evalue,
        )

        # Parse BLAST results
        tn_loc = self._parse_blast()
        info(f"Found {len(tn_loc)} TN elements via ISfinder")
        return tn_loc

    def _parse_blast(self) -> RecordTypedDf[RefTnLoc]:
        """Parse BLAST output and convert to TN locations."""
        blast_df = parse_blast_csv(self.blast_output)

        tn_lengths = read_fasta_lengths(self.isdb_path)
        df = blast_df.df.copy()

        # Remove duplicate alignments (keep first by qstart, qend)
        df = df.drop_duplicates(subset=["qstart", "qend"], keep="first")

        if df.empty:
            return RecordTypedDf.empty(RefTnLoc)

        # Add subject length
        df["subject_length"] = df["subject"].map(tn_lengths).fillna(0).astype(int)

        # Sort by query, qstart, bitscore
        df = df.sort_values(["query", "qstart", "bitscore"])

        # Filter by coverage (alignment covers >90% of TN element)
        df["coverage"] = abs(df["sstart"] - df["send"]) / df["subject_length"]
        df = df[df["coverage"] > self.critical_coverage]

        # Detect strand: sstart > send means reverse complement
        is_complement = (df["sstart"] > df["send"]).values

        # Ensure LocLeft < LocRight
        loc_left = df[["qstart", "qend"]].min(axis=1).values
        loc_right = df[["qstart", "qend"]].max(axis=1).values

        # Build tn_loc table
        return RecordTypedDf(pd.DataFrame({
            "tn_id": range(1, len(df) + 1),
            "tn_name": df["subject"].values,
            "tn_scaf": df["query"].values,
            "loc_left": loc_left,
            "loc_right": loc_right,
            "complement": is_complement,
            "join": False,
        }), RefTnLoc)

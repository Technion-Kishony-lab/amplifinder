"""Step: Find TN elements by BLAST against ISfinder database."""

from pathlib import Path
from typing import Optional

from abc import abstractmethod

from amplifinder.data_types.record_types import TnId
from amplifinder.tools.blast import run_blastn, parse_blast_csv, make_blast_db
from amplifinder.utils.fasta import read_fasta_lengths
from amplifinder.utils.file_lock import locked_resource
from amplifinder.utils.file_utils import ensure_dir
from amplifinder.data_types import Orientation, RecordTypedDf, RefTn, Genome
from amplifinder.steps.base import Step
from amplifinder.steps.locate_tns.find_tn_in_genbank import find_tn_elements

"""GenBank TN element parsing."""
from __future__ import annotations

import re
from typing import List, Optional

from Bio.SeqFeature import SeqFeature

from amplifinder.data_types.enums import Orientation
from amplifinder.data_types.record_types import RefTn, TnId
from amplifinder.data_types.genome import Genome

 
# Base class for TN location steps

class LocateTNsStep(Step[Optional[RecordTypedDf[RefTn]]]):
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

    def _save_output(self, output: Optional[RecordTypedDf[RefTn]]) -> None:
        if output is not None:
            ensure_dir(self.output_file.parent)
            output.to_csv(self.output_file)

    def load_outputs(self) -> Optional[RecordTypedDf[RefTn]]:
        """Load TN locations from output file (indexed by tn_id)."""
        if not self.output_file.exists():
            return None
        # Load with index_col=0 to preserve tn_id index
        return RecordTypedDf.from_csv(self.output_file, RefTn, index_col=0)

    @abstractmethod
    def _calculate_output(self) -> Optional[RecordTypedDf[RefTn]]:
        """Run the TN location logic."""
        pass

    def _build_ref_tns_dict(
        self,
        items: list[tuple[str, int, int, Orientation, str]],
        start_id: int = 1,
    ) -> RecordTypedDf[RefTn]:
        """Build ref_tns dict from list of (scaf, left, right, orientation, tn_name) tuples.
        
        Args:
            items: List of tuples (scaffold_name, left, right, orientation, tn_name)
            start_id: Starting tn_id (default 1)
        
        Returns:
            RecordTypedDf[RefTn] indexed by tn_id
        """
        ref_tns: dict[TnId, RefTn] = {}
        tn_id = start_id - 1
        
        for scaf, left, right, orientation, tn_name in items:
            scaffold = self.genome.get_scaffold(scaf)
            tn_id += 1
            ref_tns[tn_id] = RefTn.from_scaffold_left_right_orientation(
                scaffold, left=left, right=right, orientation=orientation,
                tn_id=tn_id, tn_name=tn_name, join=False
            )
        
        return RecordTypedDf.from_dict(ref_tns, RefTn)


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

    @staticmethod
    def _extract_tn_name_from_feature(feature: SeqFeature) -> Optional[str]:
        """Extract TN name from a GenBank feature if it's an insertion sequence.

        Returns TN name if feature is a TN element, None otherwise.
        """
        # Get text to search based on feature type
        qualifiers = feature.qualifiers
        if feature.type == "mobile_element":
            text = qualifiers.get("mobile_element_type", [""])[0]
        elif feature.type in ("misc_feature", "repeat_region"):
            text = qualifiers.get("note", [""])[0]
        else:
            return None

        # Check for "insertion sequence" and extract name
        if "insertion sequence" not in text.lower():
            return None

        # Try "insertion sequence:IS1" format
        match = re.search(r'insertion sequence[:\s]*(\S+)', text, re.IGNORECASE)
        if match:
            return match.group(1)

        # Try ISxxx pattern anywhere
        match = re.search(r'(IS\d+\w*)', text, re.IGNORECASE)
        if match:
            return match.group(1)

        return "unknown"

    def _calculate_output(self) -> Optional[RecordTypedDf[RefTn]]:
        """Parse GenBank file and extract TN locations."""
        items = []
        
        for record in self.genome.gb_records:
            scaf = record.name
            scaffold = self.genome.get_scaffold(scaf)
            for feature in record.features:
                tn_name = self._extract_tn_name_from_feature(feature)
                if tn_name is None:
                    continue

                location = feature.location
                assert 1 <= location.start <= location.end <= len(scaffold)  # GenBank start/end are left/right
                orientation = Orientation.REVERSE if location.strand == -1 else Orientation.FORWARD
                items.append((scaf, location.start, location.end, orientation, tn_name))

        return self._build_ref_tns_dict(items)

    def report_output_message(self, output: Optional[RecordTypedDf[RefTn]], *, from_cache: bool) -> Optional[str]:
        if output is None:
            return "No GenBank file provided - skipping GenBank TN annotation."
        return f"GenBank: found {len(output)} TN elements"


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

    def _calculate_output(self) -> RecordTypedDf[RefTn]:
        """Run BLAST and parse results."""
        ensure_dir(self.blast_output.parent)

        # Create BLAST DB if needed (with lock to prevent parallel creation)
        with locked_resource(self.isdb_path, "blast_db", timeout=300):
            if not self.isdb_path.with_suffix(".nhr").exists():
                self.log("Creating ISfinder BLAST database...")
                make_blast_db(self.isdb_path, self.isdb_path)

        # Run BLAST
        run_blastn(
            query=self.genome.fasta_path,
            db=self.isdb_path,
            output=self.blast_output,
            evalue=self.evalue,
        )

        # Parse BLAST results
        return self._parse_blast()

    def report_output_message(self, output: RecordTypedDf[RefTn], *, from_cache: bool) -> Optional[str]:
        return f"ISfinder: found {len(output)} TN elements"

    def _parse_blast(self) -> RecordTypedDf[RefTn]:
        """Parse BLAST output and convert to TN locations."""
        # Parse BLAST output
        blast_hits = parse_blast_csv(self.blast_output)

        tn_lengths = read_fasta_lengths(self.isdb_path)
        df = blast_hits.df.copy()

        # qstart and qend are 1-based from BLAST output, qstart < qend always
        assert (df["qstart"] < df["qend"]).all(), "BLAST query coordinates must have qstart < qend"

        # Remove duplicate alignments (keep first by qstart, qend)
        df = df.drop_duplicates(subset=["qstart", "qend"], keep="first")

        # Add subject length
        df["subject_length"] = df["subject"].map(tn_lengths).fillna(0).astype(int)

        # Sort by query, qstart, bitscore
        df = df.sort_values(["query", "qstart", "bitscore"])

        # Filter by coverage (alignment covers >90% of TN element)
        df["coverage"] = abs(df["sstart"] - df["send"]) / df["subject_length"]
        df = df[df["coverage"] > self.critical_coverage]

        # Build ref_tns table indexed by tn_id
        items = []
        for _, row in df.iterrows():
            # Detect strand: sstart > send means reverse complement
            is_complement = row["sstart"] > row["send"]
            orientation = Orientation.REVERSE if is_complement else Orientation.FORWARD
            items.append((row["query"], row["qstart"], row["qend"], orientation, row["subject"]))
        
        return self._build_ref_tns_dict(items)

"""Step: Find TN elements by BLAST against ISfinder database."""
from __future__ import annotations

from pathlib import Path
from typing import Optional
import re

from Bio.SeqFeature import SeqFeature

from amplifinder.data_types import TnId
from amplifinder.tools.blast import run_blastn, parse_blast_csv, make_blast_db
from amplifinder.utils.fasta import read_fasta_lengths
from amplifinder.utils.file_lock import locked_resource
from amplifinder.utils.file_utils import ensure_dir
from amplifinder.data_types import Orientation, RecordTypedDf, RefTn, Genome
from amplifinder.steps.base import OutputStep
from amplifinder.logger import logger


# Base class for TN location steps

class LocateTNsStep(OutputStep[Optional[RecordTypedDf[RefTn]]]):
    """Base class for steps that locate TN elements."""
    NAME = "Locate TNs"

    def __init__(
        self,
        genome: Genome,
        output_dir: Path,
        source: str,
        input_files: list[Path],
        artifact_files: Optional[list[Path]] = None,
        force: Optional[bool] = None,
    ):
        self.genome = genome
        self.output_dir = Path(output_dir)

        self.output_file = self.output_dir / f"{source}_tn_loc.csv"

        super().__init__(
            input_files=input_files,
            artifact_files=artifact_files,
            output_files=[self.output_file],
            force=force,
        )

    def _save_output(self, output: Optional[RecordTypedDf[RefTn]]) -> None:
        if output is not None:
            ensure_dir(self.output_file.parent)
            output.to_csv(self.output_file)

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
        tn_id = start_id

        for scaf, left, right, orientation, tn_name in items:
            scaffold = self.genome.get_scaffold(scaf)
            ref_tns[tn_id] = RefTn.from_scaffold_left_right_orientation(
                scaffold, left=left, right=right, orientation=orientation,
                tn_id=tn_id, tn_name=tn_name, join=False
            )
            tn_id += 1

        return RecordTypedDf.from_dict(ref_tns, RefTn)


# Locate TNs using GenBank annotations

class LocateTNsUsingGenbankStep(LocateTNsStep):
    """Extract TN elements from GenBank file annotations using BioPython.

    Parses GenBank file for 'insertion sequence' features (based on findISinRef.m).
    Returns None if no GenBank file is provided.
    """
    NAME = "locate TNs (genbank)"

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
        if self.genome.gb_records is None:
            return None

        items = []

        for record in self.genome.gb_records:
            scaf = record.name
            scaffold = self.genome.get_scaffold(scaf)
            for feature in record.features:
                tn_name = self._extract_tn_name_from_feature(feature)
                if tn_name is None:
                    continue

                location = feature.location
                # BioPython uses 0-based start, 1-based exclusive end; convert to 1-based inclusive
                left = location.start + 1  # Convert 0-based to 1-based
                right = location.end  # Already 1-based exclusive == 1-based inclusive
                assert 1 <= left <= right <= len(scaffold)
                orientation = Orientation.REVERSE if location.strand == -1 else Orientation.FORWARD
                items.append((scaf, left, right, orientation, tn_name))

        return self._build_ref_tns_dict(items)

    def report_output_message(self, output: Optional[RecordTypedDf[RefTn]]) -> Optional[str]:
        if output is None:
            return "No GenBank file provided - skipping GenBank TN annotation."
        return f"GenBank: found {len(output)} TN elements"


# Locate TNs using ISfinder database

class LocateTNsUsingISfinderStep(LocateTNsStep):
    """BLAST reference genome against ISfinder database to find TN elements."""
    NAME = "locate TNs (isfinder)"

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

        self.blast_output = Path(output_dir) / "isfinder_blast.txt"

        super().__init__(
            genome=genome,
            output_dir=output_dir,
            source="isfinder",
            input_files=[genome.fasta_path, self.isdb_path],
            artifact_files=[self.blast_output],
            force=force,
        )

    def _generate_artifacts(self) -> None:
        """Run BLAST to produce artifact files."""
        ensure_dir(self.blast_output.parent)

        # Create BLAST DB if needed (with lock to prevent parallel creation)
        with locked_resource(self.isdb_path, "blast_db", timeout=300):
            if not self.isdb_path.with_suffix(".nhr").exists():
                logger.info("Creating ISfinder BLAST database...")
                make_blast_db(self.isdb_path, self.isdb_path)

        # Run BLAST
        run_blastn(
            query=self.genome.fasta_path,
            db=self.isdb_path,
            output=self.blast_output,
            evalue=self.evalue,
        )

    def report_output_message(self, output: RecordTypedDf[RefTn]) -> Optional[str]:
        return f"ISfinder: found {len(output)} TN elements"

    def _calculate_output(self) -> RecordTypedDf[RefTn]:
        """Parse BLAST results, handling multiple hits per IS (following MATLAB blast2ISloc.m)."""
        blast_hits = parse_blast_csv(self.blast_output)
        df = blast_hits.df.copy()

        if df.empty:
            return self._build_ref_tns_dict([])

        # qstart and qend are 1-based from BLAST output, qstart < qend always
        assert (df["qstart"] < df["qend"]).all(), "BLAST query coordinates must have qstart < qend"

        # Remove duplicate alignments (keep first by qstart, qend)
        df = df.drop_duplicates(subset=["qstart", "qend"], keep="first")

        # Add subject length
        tn_lengths = read_fasta_lengths(self.isdb_path)
        df["subject_length"] = df["subject"].map(tn_lengths).fillna(0).astype(int)

        # Sort by query, qstart, bitscore
        df = df.sort_values(["query", "qstart", "bitscore"]).reset_index(drop=True)

        # Find overlapping hit groups
        # Gap between consecutive hits: qstart[i] - qend[i-1]
        df.loc[1:, "gap"] = df["qstart"].iloc[1:].values - df["qend"].iloc[:-1].values

        # Identify group boundaries (positive gaps)
        group_starts = [0] + [i for i in range(1, len(df)) if df.loc[i, "gap"] > 0] + [len(df)]

        # For each group, select best hit
        keep_mask = []
        for i in range(len(group_starts) - 1):
            group_idx = list(range(group_starts[i], group_starts[i+1]))
            group = df.iloc[group_idx]

            # Find hit with alignment length closest to subject length
            alignment_len = group["qend"] - group["qstart"]
            length_diff = abs(alignment_len - group["subject_length"])
            best_idx = length_diff.idxmin()
            keep_mask.append(best_idx)

        df = df.loc[keep_mask]

        # Filter by coverage
        df["coverage"] = abs(df["sstart"] - df["send"]) / df["subject_length"]
        df = df[df["coverage"] > self.critical_coverage]

        # Build ref_tns table
        items = []
        for _, row in df.iterrows():
            is_complement = row["sstart"] > row["send"]
            orientation = Orientation.REVERSE if is_complement else Orientation.FORWARD
            items.append((row["query"], row["qstart"], row["qend"], orientation, row["subject"]))

        return self._build_ref_tns_dict(items)

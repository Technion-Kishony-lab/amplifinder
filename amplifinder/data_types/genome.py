"""Reference genome data model."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict

from Bio import Entrez, SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord

from amplifinder.data_types.record_types import Junction, Orientation
from amplifinder.logger import info
from amplifinder.utils.file_utils import ensure_dir

# Set email for NCBI Entrez (required)
Entrez.email = "amplifinder@example.com"


@dataclass
class Genome:
    """Reference genome - can be from GenBank, FASTA, or both.

    Scaffold names are LOCUS names (GenBank) or sequence IDs (FASTA).
    """

    name: str
    genbank_path: Optional[Path] = None
    fasta_path: Optional[Path] = None

    # Cached records
    _gb_records: Optional[List[SeqRecord]] = field(default=None, repr=False, compare=False)
    _fa_records: Optional[List[SeqRecord]] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.genbank_path is None and self.fasta_path is None:
            raise ValueError("At least one of genbank_path or fasta_path required")

    def __len__(self) -> int:
        """Get total genome length (sum of all scaffold lengths)."""
        return sum(self.scaffold_lengths.values())

    @property
    def gb_records(self) -> Optional[List[SeqRecord]]:
        """Load GenBank records (cached). Returns None if no genbank_path."""
        if not self.genbank_path:
            return None
        if self._gb_records is None:
            with open(self.genbank_path) as f:
                self._gb_records = list(SeqIO.parse(f, "genbank"))
        return self._gb_records

    @property
    def fa_records(self) -> Optional[List[SeqRecord]]:
        """Load FASTA records (cached). Returns None if no fasta_path."""
        if not self.fasta_path:
            return None
        if self._fa_records is None:
            with open(self.fasta_path) as f:
                self._fa_records = list(SeqIO.parse(f, "fasta"))
        return self._fa_records

    @property
    def records(self) -> List[SeqRecord]:
        """Load records from genome file (GenBank preferred, FASTA fallback)."""
        if self.gb_records is not None:
            return self.gb_records
        if self.fa_records is not None:
            return self.fa_records
        raise ValueError(f"No records available for {self.name}")

    @property
    def scaffold_circularities(self) -> Dict[str, bool]:
        """Get circularity status per scaffold (from GenBank annotations).

        Returns:
            Dict mapping scaffold name to circularity (True if circular, False if linear).
            For FASTA-only genomes with a single contig, assumes circular=True.
        """
        records = self.records
        if self.genbank_path:
            return {
                rec.name: rec.annotations.get("topology", "linear").lower() == "circular"
                for rec in records
            }
        # FASTA-only: assume circular if single contig (common for bacterial genomes)
        # TODO: Logic for FASTA-only genomes ok?
        if len(records) == 1:
            return {records[0].name: True}
        return {rec.name: False for rec in records}

    @property
    def scaffold_sequences(self) -> Dict[str, str]:
        """Get all sequences as {scaffold_name: sequence} dict."""
        return {rec.name: str(rec.seq) for rec in self.records}

    @property
    def scaffold_lengths(self) -> Dict[str, int]:
        """Get scaffold lengths as {scaffold_name: length} dict."""
        return {name: len(seq) for name, seq in self.scaffold_sequences.items()}

    @property
    def scaffold_ranges(self) -> Dict[str, tuple[int, int]]:
        """Get scaffold ranges as {scaffold_name: (start, end)} dict.

        Returns:
            Dict mapping scaffold name to (start, end) tuple of 0-based positions
            in concatenated coverage array.
        """
        ranges = {}
        cumulative = 0
        for scaffold_name, length in self.scaffold_lengths.items():
            ranges[scaffold_name] = (cumulative, cumulative + length)
            cumulative += length
        return ranges

    def get_fowrard_sequence_in_range(self, scaf: str, start: int, end: int) -> str:
        """
        Get sequence in a range. 1-based inclusive.
        Uses modulo arithmetic to handle circular genome.
        """
        scaf_circular = self.scaffold_circularities[scaf]
        scaf_seq = self.scaffold_sequences[scaf]
        if not scaf_circular:
            return scaf_seq[start - 1:end]
        start = (start - 1) % len(scaf_seq) + 1
        end = (end - 1) % len(scaf_seq) + 1
        if start <= end:
            return scaf_seq[start - 1:end]
        else:
            return scaf_seq[start - 1:] + scaf_seq[:end]

    def get_sequence_in_range(self, scaf: str, start: int, end: int, direction: Orientation) -> str:
        """Get sequence in a range."""
        if direction == Orientation.FORWARD:
            return self.get_fowrard_sequence_in_range(scaf, start, end)
        else:
            return reverse_complement(self.get_fowrard_sequence_in_range(scaf, end, start))

    def get_junction_arm_sequence(self, jc: Junction, arm: int) -> str:
        """Get sequence for a junction arm."""
        scaf, pos, direction, flank_len = jc.get_scaf_pos_dir_flank(arm)
        return self.get_sequence_in_range(scaf, pos, pos + flank_len * direction, direction)

    def get_junction_sequence(self, jc: Junction) -> str:
        """Get sequence for a junction."""
        return reverse_complement(self.get_junction_arm_sequence(jc, 1)) + self.get_junction_arm_sequence(jc, 2)

class GenomeRegistry:
    """Registry for fetching and caching genomes.

    Handles:
    - NCBI fetching
    - Local caching with accession → locus name mapping
    - File organization

    File structure:
        ref_path/
        ├── {accession}.json    # mapping: accession → locus name
        ├── genbank/
        │   └── {locus}.gb      # GenBank files by locus name
        └── fasta/
            └── {locus}.fasta   # FASTA files by locus name
    """

    def __init__(self, ref_path: Path):
        """Initialize registry with base path for genome storage."""
        self.ref_path = Path(ref_path)
        self.genbank_dir = self.ref_path / "genbank"
        self.fasta_dir = self.ref_path / "fasta"

    def _mapping_file(self, accession: str) -> Path:
        """Get path to mapping file for accession."""
        return self.ref_path / f"{accession}.json"

    def get(self, ref_name: str, ncbi: bool = True) -> Genome:
        """Get genome by name (from cache or NCBI).

        Args:
            ref_name: Reference name or NCBI accession
            ncbi: If True, fetch from NCBI if not cached
        """
        # Try to load from cache
        genome = self._load_cached(ref_name)
        if genome:
            return genome

        # Fetch from NCBI
        if ncbi:
            return self._fetch_from_ncbi(ref_name)

        raise FileNotFoundError(f"Genome {ref_name} not found in {self.ref_path}")

    def _load_cached(self, ref_name: str) -> Optional[Genome]:
        """Try to load genome from cache (GenBank or FASTA)."""
        # Try: mapping file → locus name, or ref_name directly
        locus_name = self._read_mapping(ref_name) or ref_name

        genbank_file = self.genbank_dir / f"{locus_name}.gb"
        fasta_file = self.fasta_dir / f"{locus_name}.fasta"

        has_genbank = genbank_file.exists()
        has_fasta = fasta_file.exists()

        if has_genbank or has_fasta:
            info(f"Reference {ref_name} loaded from cache")
            # Ensure mapping file exists
            if not self._mapping_file(ref_name).exists():
                self._write_mapping(ref_name, locus_name)
            return Genome(
                name=locus_name,
                genbank_path=genbank_file if has_genbank else None,
                fasta_path=fasta_file if has_fasta else None,
            )

        return None

    def _read_mapping(self, accession: str) -> Optional[str]:
        """Read locus name from mapping file, or None if not found."""
        mapping_file = self._mapping_file(accession)
        if not mapping_file.exists():
            return None
        with open(mapping_file) as f:
            return json.load(f)["name"]

    def _write_mapping(self, accession: str, locus_name: str) -> None:
        """Write accession → locus name mapping."""
        mapping_file = self._mapping_file(accession)
        with open(mapping_file, "w") as f:
            json.dump({"accession": accession, "name": locus_name}, f, indent=2)

    def _fetch_from_ncbi(self, ref_name: str) -> Genome:
        """Fetch genome from NCBI and cache locally."""
        info(f"Fetching {ref_name} from NCBI...")
        ensure_dir(self.genbank_dir)

        # Download GenBank
        genbank_file = self.genbank_dir / f"{ref_name}.gb"
        with Entrez.efetch(db="nuccore", id=ref_name, rettype="gb", retmode="text") as handle:
            genbank_text = handle.read()
        with open(genbank_file, "w") as f:
            f.write(genbank_text)

        record = SeqIO.read(genbank_file, "genbank")
        locus_name = record.name
        info(f"Fetched {locus_name}: {len(record.seq):,} bp")

        # Rename to locus name if different
        if locus_name != ref_name:
            info(f"Reference locus name: {locus_name}")
            actual_genbank = self.genbank_dir / f"{locus_name}.gb"
            if not actual_genbank.exists():
                genbank_file.rename(actual_genbank)
            genbank_file = actual_genbank

        # Always write mapping file
        self._write_mapping(ref_name, locus_name)

        # Create FASTA
        fasta_file = self.fasta_dir / f"{locus_name}.fasta"
        if not fasta_file.exists():
            info(f"Creating FASTA for {locus_name}")
            ensure_dir(self.fasta_dir)
            record.id = locus_name
            record.description = ""
            SeqIO.write(record, fasta_file, "fasta")

        return Genome(
            name=locus_name,
            genbank_path=genbank_file.resolve(),
            fasta_path=fasta_file.resolve(),
            _gb_records=[record],
        )


def get_genome(ref_name: str, ref_path: Path, ncbi: bool = True) -> Genome:
    """Get genome by name (convenience wrapper for GenomeRegistry)."""
    registry = GenomeRegistry(ref_path)
    return registry.get(ref_name, ncbi=ncbi)

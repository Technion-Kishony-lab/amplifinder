"""Reference genome data model."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict

from Bio import Entrez, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature

from amplifinder.logger import info

# Set email for NCBI Entrez (required)
Entrez.email = "amplifinder@example.com"


@dataclass
class Genome:
    """Reference genome - can be from GenBank, FASTA, or both."""

    name: str
    genbank_path: Optional[Path] = None
    fasta_path: Optional[Path] = None

    # Cached records
    _gb_record: Optional[SeqRecord] = field(default=None, repr=False, compare=False)
    _fa_record: Optional[SeqRecord] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.genbank_path is None and self.fasta_path is None:
            raise ValueError("At least one of genbank_path or fasta_path required")

    @property
    def genbank_record(self) -> Optional[SeqRecord]:
        """Load and cache GenBank record."""
        if self._gb_record is None and self.genbank_path is not None:
            self._gb_record = SeqIO.read(self.genbank_path, "genbank")
        return self._gb_record

    @property
    def fasta_record(self) -> Optional[SeqRecord]:
        """Load and cache FASTA record."""
        if self._fa_record is None and self.fasta_path is not None:
            self._fa_record = SeqIO.read(self.fasta_path, "fasta")
        return self._fa_record

    @property
    def seq(self) -> Seq:
        """Get sequence from available source (validates if both exist)."""
        gb_seq = self.genbank_record.seq if self.genbank_record else None
        fa_seq = self.fasta_record.seq if self.fasta_record else None

        if gb_seq is not None and fa_seq is not None:
            assert str(gb_seq) == str(fa_seq), f"Sequence mismatch for {self.name}!"
            return gb_seq
        elif gb_seq is not None:
            return gb_seq
        elif fa_seq is not None:
            return fa_seq
        else:
            raise ValueError(f"No sequence available for {self.name}")

    @property
    def sequence(self) -> str:
        """Get genome sequence as string."""
        return str(self.seq)

    @property
    def length(self) -> int:
        """Get genome length."""
        return len(self.seq)

    @property
    def circular(self) -> bool:
        """Check if genome is circular (from GenBank annotations)."""
        if self.genbank_record:
            topology = self.genbank_record.annotations.get("topology", "linear")
            return topology.lower() == "circular"
        return False

    def get_subsequence(self, start: int, end: int) -> str:
        """Get subsequence (1-based coordinates, inclusive)."""
        return str(self.seq[start - 1:end])

    @property
    def record(self) -> SeqRecord:
        """Get primary record (GenBank preferred)."""
        if self.genbank_record:
            return self.genbank_record
        if self.fasta_record:
            return self.fasta_record
        raise ValueError(f"No record available for {self.name}")

    @property
    def features(self) -> List[SeqFeature]:
        """Get GenBank features (None if FASTA-only)."""
        if self.genbank_record:
            return self.genbank_record.features
        return []

    @property
    def annotations(self) -> dict:
        """Get GenBank annotations (empty dict if FASTA-only)."""
        if self.genbank_record:
            return self.genbank_record.annotations
        return {}

    @property
    def records(self) -> List[SeqRecord]:
        """Load all records from genome file."""
        if self.genbank_path:
            return list(SeqIO.parse(self.genbank_path, "genbank"))
        if self.fasta_path:
            return list(SeqIO.parse(self.fasta_path, "fasta"))
        return []

    @property
    def sequences(self) -> Dict[str, str]:
        """Get all sequences as {scaffold_name: sequence} dict.

        Uses record.name (locus name) as key to match TN scaffold names.
        """
        return {rec.name: str(rec.seq) for rec in self.records}


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
        self.genbank_dir.mkdir(parents=True, exist_ok=True)

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
            self.fasta_dir.mkdir(parents=True, exist_ok=True)
            record.id = locus_name
            record.description = ""
            SeqIO.write(record, fasta_file, "fasta")

        return Genome(
            name=locus_name,
            genbank_path=genbank_file.resolve(),
            fasta_path=fasta_file.resolve(),
            _gb_record=record,
        )


def get_genome(ref_name: str, ref_path: Path, ncbi: bool = True) -> Genome:
    """Get genome by name (convenience wrapper for GenomeRegistry)."""
    registry = GenomeRegistry(ref_path)
    return registry.get(ref_name, ncbi=ncbi)

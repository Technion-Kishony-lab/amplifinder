"""Reference genome data model."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict, TypeVar

import numpy as np
from Bio import Entrez, SeqIO
from Bio.Seq import reverse_complement
from Bio.SeqRecord import SeqRecord

from amplifinder.data_types.scaffold import SeqScaffold, SeqSegmentScaffold
from amplifinder.data_types.record_types import Junction
from amplifinder.logger import info
from amplifinder.utils.file_utils import ensure_dir

T = TypeVar('T', str, np.ndarray)

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
    _scaffolds: Optional[Dict[str, SeqScaffold]] = field(default=None, repr=False, compare=False)

    def __post_init__(self):
        if self.genbank_path is None and self.fasta_path is None:
            raise ValueError("At least one of genbank_path or fasta_path required")

    def __len__(self) -> int:
        """Get total genome length (sum of all scaffold lengths)."""
        return sum(len(scaf) for scaf in self.scaffolds.values())

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
    def scaffolds(self) -> Dict[str, SeqScaffold]:
        """Get scaffold objects keyed by name."""
        if self._scaffolds is None:
            records = self.records
            if self.genbank_path:
                circularities = {
                    rec.name: rec.annotations.get("topology", "linear").lower() == "circular"
                    for rec in records
                }
            elif len(records) == 1:
                circularities = {records[0].name: True}
            else:
                circularities = {rec.name: False for rec in records}
            self._scaffolds = {
                rec.name: SeqScaffold(seq=str(rec.seq), is_circular=circularities[rec.name])
                for rec in records
            }
        return self._scaffolds

    @property
    def scaffold_ranges(self) -> Dict[str, tuple[int, int]]:
        """Get scaffold ranges as {scaffold_name: (start, end)} dict.

        Returns:
            Dict mapping scaffold name to (start, end) tuple of 0-based positions
            in concatenated coverage array.
        """
        ranges = {}
        cumulative = 0
        for scaffold_name, scaf in self.scaffolds.items():
            ranges[scaffold_name] = (cumulative, cumulative + len(scaf))
            cumulative += len(scaf)
        return ranges

    def get_scaffold(self, scaf: str) -> SeqScaffold:
        """Get scaffold by name."""
        return self.scaffolds[scaf]

    def get_junction_arm_sequence(self, jc: Junction, arm: int) -> str:
        """Get sequence for a junction arm."""
        jc_arm = jc.get_jc_arm(arm)
        scaffold = self.get_scaffold(jc_arm.scaf)
        seg_scaffold = SeqSegmentScaffold.from_scaffold_and_jc_arm(scaffold, jc_arm)
        return seg_scaffold.slice()

    def get_junction_sequence_arm1_to_arm2(self, jc: Junction) -> str:
        """Get sequence for a junction from arm 1 to arm 2."""
        return reverse_complement(self.get_junction_arm_sequence(jc, 1)) + self.get_junction_arm_sequence(jc, 2)

    def get_junction_sequence_arm2_to_arm1(self, jc: Junction) -> str:
        """Get sequence for a junction from arm 2 to arm 1."""
        return reverse_complement(self.get_junction_arm_sequence(jc, 2)) + self.get_junction_arm_sequence(jc, 1)


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
            info(f"Reference {ref_name} loaded from cache")
            return genome

        # Fetch from NCBI
        if ncbi:
            info(f"Fetching {ref_name} from NCBI...")
            return self._fetch_from_ncbi(ref_name)

        raise FileNotFoundError(f"Genome {ref_name} not found in {self.ref_path}")

    def exists(self, ref_name: str) -> bool:
        """Return True if genome is cached locally."""
        return self._load_cached(ref_name) is not None

    def _load_cached(self, ref_name: str) -> Optional[Genome]:
        """Try to load genome from cache (GenBank or FASTA)."""
        # Try: mapping file → locus name, or ref_name directly
        locus_name = self._read_mapping(ref_name) or ref_name

        genbank_file = self.genbank_dir / f"{locus_name}.gb"
        fasta_file = self.fasta_dir / f"{locus_name}.fasta"

        has_genbank = genbank_file.exists()
        has_fasta = fasta_file.exists()

        if has_genbank or has_fasta:
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
    return GenomeRegistry(ref_path).get(ref_name, ncbi=ncbi)


def exists_genome(ref_name: str, ref_path: Path) -> bool:
    """Return True if genome is cached locally."""
    return GenomeRegistry(ref_path).exists(ref_name)

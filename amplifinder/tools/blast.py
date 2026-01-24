"""BLAST runner and parser."""

from pathlib import Path
from typing import ClassVar, Dict, List, Optional

from amplifinder.logger import info
from amplifinder.records.base_records import Record
from amplifinder.data_types import RecordTypedDf
from amplifinder.env import BLAST_PATH
from amplifinder.utils.run_utils import get_tool_path, run_command


class BlastHit(Record):
    """BLAST alignment hit record."""
    NAME: ClassVar[str] = "BLAST hits"
    CSV_FIELD_FORMATS: ClassVar[Dict[str, str]] = {
        'evalue': '.6e',  # Scientific notation for e-values
    }
    query: str
    subject: str
    percent_identical: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: float
    bitscore: float


def make_blast_db(
    fasta: Path,
    db_path: Path,
    dbtype: str = "nucl",
) -> None:
    """Create a BLAST database from a FASTA file."""
    makeblastdb = get_tool_path("makeblastdb", config_path=BLAST_PATH)
    cmd = [
        makeblastdb,
        "-in", str(fasta),
        "-dbtype", dbtype,
        "-out", str(db_path),
    ]
    info(f"Creating BLAST DB: {' '.join(str(c) for c in cmd)}")
    run_command(cmd, check=True)


def run_blastn(
    query: Path,
    db: Path,
    output: Path,
    evalue: float = 1e-4,
    num_alignments: int = 10000,
    outfmt: int = 10,
    extra_args: Optional[List[str]] = None,
) -> None:
    """Run blastn against a database.

    Args:
        query: Query FASTA file
        db: BLAST database path
        output: Output file path
        evalue: E-value threshold
        num_alignments: Max alignments to report
        outfmt: Output format (10 = CSV)
        extra_args: Additional blastn arguments
    """
    blastn = get_tool_path("blastn", config_path=BLAST_PATH)
    cmd = [
        blastn,
        "-db", str(db),
        "-query", str(query),
        "-evalue", str(evalue),
        "-num_alignments", str(num_alignments),
        "-outfmt", str(outfmt),
        "-out", str(output),
    ]
    if extra_args:
        cmd.extend(extra_args)

    info(f"Running BLAST: {' '.join(str(c) for c in cmd)}")
    run_command(cmd, check=True)


def parse_blast_csv(path: Path) -> RecordTypedDf[BlastHit]:
    """Parse BLAST CSV output (format 10).

    Args:
        path: Path to BLAST output file

    Returns:
        RecordDF with BLAST results, empty with correct schema if file missing/empty
    """
    return RecordTypedDf.from_csv(path, BlastHit, headers=False)

"""BLAST runner and parser."""

import subprocess
from pathlib import Path
from typing import List, Optional

import pandas as pd

from amplifinder.logger import info
from amplifinder.data_types.schema import BLAST_SCHEMA
from amplifinder.env import BLAST_PATH


def _get_blast_cmd(cmd_name: str) -> str:
    """Get full path to BLAST command."""
    if BLAST_PATH:
        return str(BLAST_PATH / cmd_name)
    return cmd_name


def make_blast_db(
    fasta: Path,
    db_path: Path,
    dbtype: str = "nucl",
) -> None:
    """Create a BLAST database from a FASTA file."""
    cmd = [
        _get_blast_cmd("makeblastdb"),
        "-in", str(fasta),
        "-dbtype", dbtype,
        "-out", str(db_path),
    ]
    info(f"Creating BLAST DB: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


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
    cmd = [
        _get_blast_cmd("blastn"),
        "-db", str(db),
        "-query", str(query),
        "-evalue", str(evalue),
        "-num_alignments", str(num_alignments),
        "-outfmt", str(outfmt),
        "-out", str(output),
    ]
    if extra_args:
        cmd.extend(extra_args)
    
    info(f"Running BLAST: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def parse_blast_csv(path: Path) -> pd.DataFrame:
    """Parse BLAST CSV output (format 10).
    
    Args:
        path: Path to BLAST output file
    
    Returns:
        DataFrame with BLAST results, empty with correct schema if file missing/empty
    """
    return BLAST_SCHEMA.read_csv(path, headers=False)

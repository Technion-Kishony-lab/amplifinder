"""FASTA/FASTQ file utilities."""

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

from Bio.Seq import Seq

from amplifinder.logger import info, warning


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence.
    
    Args:
        seq: DNA sequence string
    
    Returns:
        Reverse complement sequence
    """
    return str(Seq(seq).reverse_complement())


def read_fasta_lengths(fasta_path: Path, max_num_reads: Optional[int] = None) -> Dict[str, int]:
    """Read FASTA file and return sequence lengths.

    Args:
        fasta_path: Path to FASTA file
        max_num_reads: Maximum number of sequences to read (None = all)

    Returns:
        Dict mapping sequence ID to length
    """
    lengths = {}
    current_header = None
    current_len = 0

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if current_header:
                    lengths[current_header] = current_len
                    if max_num_reads and len(lengths) >= max_num_reads:
                        return lengths
                current_header = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_header and (not max_num_reads or len(lengths) < max_num_reads):
            lengths[current_header] = current_len

    return lengths


def read_fastq_lengths(fastq_path: Path, max_num_reads: Optional[int] = None) -> List[int]:
    """Read FASTQ file and return read lengths.

    Args:
        fastq_path: Path to FASTQ file (.fastq or .fastq.gz)
        max_num_reads: Maximum number of reads to sample (None = all)

    Returns:
        List of read lengths
    """
    lengths = []
    opener = gzip.open if str(fastq_path).endswith('.gz') else open
    with opener(fastq_path, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                lengths.append(len(line.strip()))
                if max_num_reads and len(lengths) >= max_num_reads:
                    break
    return lengths


@dataclass
class ReadLengthStats:
    """Read length statistics from FASTQ files."""
    max_length: int
    mean_length: float
    is_uniform: bool  # All sampled reads within 5% of each other
    num_bases: int
    num_reads: int


def get_read_length_stats(
    fastq_dir: Path,
    sample_size: int = 1000,
) -> ReadLengthStats:
    """Calculate read length statistics from FASTQ files in a directory.

    Args:
        fastq_dir: Directory containing FASTQ files
        sample_size: Number of reads to sample per file for uniformity check

    Returns:
        ReadLengthStats with max_length, mean_length, is_uniform, num_bases, num_reads
    """
    fastq_dir = Path(fastq_dir)

    # Find FASTQ files
    fastq_files = list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fastq"))
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")

    all_lengths: List[int] = []
    total_bases = 0
    total_reads = 0

    for fq_file in fastq_files:
        lengths = read_fastq_lengths(fq_file, max_num_reads=sample_size)
        if lengths:
            all_lengths.extend(lengths)
            # Estimate total bases/reads from file
            file_lengths = read_fastq_lengths(fq_file, max_num_reads=None)
            total_reads += len(file_lengths)
            total_bases += sum(file_lengths)

    if not all_lengths:
        raise ValueError(f"No reads found in FASTQ files in {fastq_dir}")

    max_length = max(all_lengths)
    mean_length = sum(all_lengths) / len(all_lengths)

    # Check uniformity: all sampled reads within 5% of each other
    min_len = min(all_lengths)
    is_uniform = (max_length / min_len) < 1.05 if min_len > 0 else False

    info(f"Read length stats: max={max_length}, mean={mean_length:.1f}, "
         f"uniform={is_uniform}, reads={total_reads:,}, bases={total_bases:,}")

    return ReadLengthStats(
        max_length=max_length,
        mean_length=mean_length,
        is_uniform=is_uniform,
        num_bases=total_bases,
        num_reads=total_reads,
    )


def get_read_length(
    fastq_dir: Optional[Path] = None,
    provided_length: Optional[int] = None,
    sample_size: int = 1000,
) -> int:
    """Get read length from provided value or calculate from FASTQ files.

    Args:
        fastq_dir: Directory containing FASTQ files (used if provided_length is None)
        provided_length: Pre-specified read length (takes precedence)
        sample_size: Number of reads to sample for calculation

    Returns:
        Read length (max length from FASTQ or provided value)

    Raises:
        ValueError: If neither fastq_dir nor provided_length is given
    """
    if provided_length is not None:
        info(f"Using provided read length: {provided_length}")
        return provided_length

    if fastq_dir is None:
        raise ValueError("Either fastq_dir or provided_length must be specified")

    stats = get_read_length_stats(fastq_dir, sample_size)
    if not stats.is_uniform:
        warning("Read length is not uniform - may affect junction coverage accuracy")
    return stats.max_length

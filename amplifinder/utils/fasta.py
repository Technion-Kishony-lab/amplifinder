"""FASTA/FASTQ file utilities."""

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Union

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from amplifinder.logger import info, warning
from amplifinder.utils.file_utils import ensure_parent_dir


def read_fasta_lengths(fasta_path: Path, max_num_reads: Optional[int] = None) -> Dict[str, int]:
    """Read FASTA file and return sequence lengths.

    Args:
        fasta_path: Path to FASTA file
        max_num_reads: Maximum number of sequences to read (None = all)

    Returns:
        Dict mapping sequence ID to length
    """
    lengths = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        lengths[record.id] = len(record.seq)
        if max_num_reads and len(lengths) >= max_num_reads:
            break
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
        for record in SeqIO.parse(f, "fastq"):
            lengths.append(len(record.seq))
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


def write_fasta(
    sequences: Dict[str, str],
    output_path: Path,
    sort_keys: bool = False,
) -> None:
    """Write sequences to FASTA file.

    Args:
        sequences: Dict mapping sequence ID to sequence string
        output_path: Output FASTA file path
        sort_keys: If True, sort sequences by ID before writing
    """
    ensure_parent_dir(output_path)

    items = sequences.items()
    if sort_keys:
        items = sorted(items)
    records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in items]
    SeqIO.write(records, output_path, "fasta")


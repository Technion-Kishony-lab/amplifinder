"""FASTA/FASTQ file utilities."""

import gzip
import subprocess

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from amplifinder.utils.file_utils import ensure_parent_dir


READ_LENGTH_UNIFORMITY_TOLERANCE: float = 0.05


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
    min_length: int
    max_length: int
    mean_length: float

    @property
    def is_uniform(self) -> bool:
        """True if all sampled reads within 5% of each other."""
        return self.min_length > 0 and (self.max_length / self.min_length) < (1 + READ_LENGTH_UNIFORMITY_TOLERANCE)


def get_read_length_stats(
    fastq_dir: Path,
    sample_per_file: int = 1000,
) -> ReadLengthStats:
    """Calculate read length statistics from FASTQ files in a directory.

    Samples reads from all FASTQ files (e.g., paired-end R1/R2) to determine
    read length. For paired-end, R1 and R2 typically have the same length.

    Args:
        fastq_dir: Directory containing FASTQ files
        sample_per_file: Number of reads to sample per file

    Returns:
        ReadLengthStats with min_length, max_length, mean_length
    """
    fastq_dir = Path(fastq_dir)

    # Find FASTQ files
    fastq_files = list(fastq_dir.glob("*.fastq.gz")) + list(fastq_dir.glob("*.fastq"))
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_dir}")

    all_lengths: List[int] = []
    for fq_file in fastq_files:
        lengths = read_fastq_lengths(fq_file, max_num_reads=sample_per_file)
        all_lengths.extend(lengths)

    if not all_lengths:
        raise ValueError(f"No reads found in FASTQ files in {fastq_dir}")

    return ReadLengthStats(
        min_length=min(all_lengths),
        max_length=max(all_lengths),
        mean_length=sum(all_lengths) / len(all_lengths),
    )


def write_fasta(
    sequences: List[tuple[str, str]],
    output_path: Path,
) -> None:
    """Write sequences to FASTA file.

    Args:
        sequences: List of (seq_id, sequence) tuples
        output_path: Output FASTA file path
        sort_keys: If True, sort sequences by ID before writing
    """
    ensure_parent_dir(output_path)
    records = [SeqRecord(Seq(seq), id=seq_id, description="") for seq_id, seq in sequences]
    SeqIO.write(records, output_path, "fasta")


def count_fastq_reads(fastq_path: Path) -> int:
    """Count number of reads in FASTQ file (lines / 4).
    
    Fast implementation using line counting, not parsing sequences.
    """
    
    try:
        if str(fastq_path).endswith('.gz'):
            result = subprocess.run(['zcat', str(fastq_path), '|', 'wc', '-l'], 
                                  shell=True, capture_output=True, text=True)
        else:
            result = subprocess.run(['wc', '-l', str(fastq_path)], 
                                  capture_output=True, text=True)
        num_lines = int(result.stdout.split()[0])
        return num_lines // 4
    except (subprocess.CalledProcessError, ValueError, IndexError, OSError):
        # Fallback: count by parsing (slower)
        count = 0
        opener = gzip.open if str(fastq_path).endswith('.gz') else open
        with opener(fastq_path, 'rt') as f:
            for record in SeqIO.parse(f, "fastq"):
                count += 1
        return count


def count_total_bases(fastq_dir: Path, read_length: int) -> int:
    """Count total bases in FASTQ directory.
    
    Args:
        fastq_dir: Directory containing FASTQ files
        read_length: Read length (from get_read_length_stats)
    
    Returns:
        Total number of bases
    """
    fastq_files = list(Path(fastq_dir).glob("*.fastq.gz")) + list(Path(fastq_dir).glob("*.fastq"))
    total_reads = sum(count_fastq_reads(fq_file) for fq_file in fastq_files)
    return total_reads * read_length

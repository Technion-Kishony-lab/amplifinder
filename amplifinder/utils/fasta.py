"""FASTA/FASTQ file utilities."""

from pathlib import Path
from typing import Dict, List, Optional


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
        fastq_path: Path to FASTQ file
        max_num_reads: Maximum number of reads to sample (None = all)
    
    Returns:
        List of read lengths
    """
    lengths = []
    with open(fastq_path) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                lengths.append(len(line.strip()))
                if max_num_reads and len(lengths) >= max_num_reads:
                    break
    return lengths


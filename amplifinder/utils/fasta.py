"""FASTA file utilities."""

from pathlib import Path


def read_fasta_lengths(fasta_path: Path) -> dict:
    """Read FASTA file and return sequence lengths.
    
    Args:
        fasta_path: Path to FASTA file
    
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
                current_header = line[1:].strip().split()[0]
                current_len = 0
            else:
                current_len += len(line.strip())
        if current_header:
            lengths[current_header] = current_len
    
    return lengths


"""Sequence manipulation utilities."""

from typing import TypeVar

import numpy as np

T = TypeVar('T', str, np.ndarray)


def concatenate_sequence(seq1: T, seq2: T) -> T:
    """Concatenate two sequences."""
    if isinstance(seq1, np.ndarray) and isinstance(seq2, np.ndarray):
        return np.concatenate([seq1, seq2])
    if isinstance(seq1, str) and isinstance(seq2, str):
        return seq1 + seq2
    if isinstance(seq1, list) and isinstance(seq2, list):
        return seq1 + seq2
    raise ValueError("Invalid sequence type")

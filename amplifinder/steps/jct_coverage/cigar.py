"""CIGAR string processing utilities."""

from __future__ import annotations
from typing import Iterator


class Cigar(list[tuple[int, int]]):
    """CIGAR operations with reference position tracking."""
    
    def iter_with_ref_pos(self, start_pos: int = 0) -> Iterator[tuple[int, int, int]]:
        """Iterate yielding (op, length, ref_pos) where ref_pos is updated for each operation.
        
        Args:
            start_pos: Initial reference position (default 0)
            
        Yields:
            (op, length, ref_pos) where ref_pos is the position before consuming this operation
        """
        ref_pos = start_pos
        for op, length in self:
            yield op, length, ref_pos
            # Update ref_pos based on operation type
            if op in (0, 2, 7, 8):  # M, D, =, X consume reference
                ref_pos += length
            # I (1) does not consume reference
    
    def has_only_operations(self, allowed_ops: set[int]) -> bool:
        """Check if all operations are in the allowed set."""
        return all(op in allowed_ops for op, _ in self)


def merge_consecutive_cigar_ops(cigar: Cigar) -> Cigar:
    """Merge consecutive identical CIGAR operations."""
    merged = Cigar()
    if not cigar:
        return merged
    
    merged.append(cigar[0])
    for op, length in cigar[1:]:
        if op == merged[-1][0]:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    
    return merged


def resolve_cigar_m_operations(cigar: Cigar, query_seq: str, ref_seq: str) -> Cigar:
    """Convert M operations to = (match) and X (mismatch) by comparing sequences."""
    if not any(op == 0 for op, _ in cigar):
        # No M operations, return as is
        return cigar
    
    query_seq = query_seq.upper()
    ref_seq = ref_seq.upper()
    
    expanded = Cigar()
    query_idx = 0
    ref_idx = 0
    
    for op, length in cigar:
        if op == 0:  # M - expand to =/X
            for i in range(length):
                if query_seq[query_idx + i] == ref_seq[ref_idx + i]:
                    expanded.append((7, 1))  # = (match)
                else:
                    expanded.append((8, 1))  # X (mismatch)
        else:
            expanded.append((op, length))

        if op == 1:  # I - insertion
            query_idx += length
        elif op == 2:  # D - deletion
            ref_idx += length
        elif op in (0, 7, 8):  # =, X - match/mismatch
            query_idx += length
            ref_idx += length
        else:
            assert False, "Unsupported CIGAR operation"
    
    return merge_consecutive_cigar_ops(expanded)

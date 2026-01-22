"""CIGAR string processing utilities."""


def merge_consecutive_cigar_ops(cigar: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge consecutive identical CIGAR operations."""
    if not cigar:
        return []
    
    merged = [cigar[0]]
    for op, length in cigar[1:]:
        if op == merged[-1][0]:
            merged[-1] = (op, merged[-1][1] + length)
        else:
            merged.append((op, length))
    
    return merged


def resolve_cigar_m_operations(cigar: list[tuple[int, int]], query_seq: str, ref_seq: str) -> list[tuple[int, int]]:
    """Convert M operations to = (match) and X (mismatch) by comparing sequences."""
    if not any(op == 0 for op, _ in cigar):
        # No M operations, return as is
        return cigar
    
    query_seq = query_seq.upper()
    ref_seq = ref_seq.upper()
    
    expanded = []
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

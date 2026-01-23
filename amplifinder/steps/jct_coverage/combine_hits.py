import pysam
from pathlib import Path
from collections import defaultdict

from amplifinder.steps.jct_coverage.alignment_data import SingleAlignment, CombinedSingleAlignment
from amplifinder.steps.jct_coverage.cigar import Cigar


def combine_same_id_same_orientation_hits(
    bam_path: Path,
) -> dict[tuple[str, bool, str], SingleAlignment]:
    """Combine hits with same read_id and same orientation into single alignments.
    
    Args:
        bam_path: Path to BAM file
        
    Returns:
        Dictionary mapping (read_id, is_reverse, reference_name) to combined SingleAlignment
    """
    # Group hits by (read_id, is_reverse, reference_name)
    grouped_hits: dict[tuple[str, bool, str], list[tuple[int, pysam.AlignedSegment]]] = defaultdict(list)
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for bam_index, hit in enumerate(bam.fetch(), start=1):
            if hit.is_unmapped or hit.is_supplementary:
                continue
            key = (hit.query_name, hit.is_reverse, hit.reference_name)
            grouped_hits[key].append((bam_index, hit))

    # Combine each group into single alignment
    combined_alignments = {}
    for key, indexed_hits in grouped_hits.items():
        if len(indexed_hits) == 1:
            # Single hit - just extract info
            bam_index, hit = indexed_hits[0]
            cigar = Cigar(hit.cigartuples)
            combined_alignments[key] = SingleAlignment(
                start=hit.reference_start,
                end=hit.reference_end,
                cigar=cigar,
                bam_index=bam_index,
            )
        else:
            # Multiple hits - combine them
            combined_alignments[key] = _combine_multiple_hits(indexed_hits)
    
    return combined_alignments


def _combine_multiple_hits(indexed_hits: list[tuple[int, pysam.AlignedSegment]]) -> CombinedSingleAlignment:
    """Combine multiple hits into a single alignment spanning leftmost to rightmost.
    
    Args:
        indexed_hits: List of (bam_index, hit) with same read_id, orientation, and reference
        
    Returns:
        CombinedSingleAlignment with overall span and individual hit info (separate CIGARs)
    """
    # Sort by reference start position
    sorted_indexed_hits = sorted(indexed_hits, key=lambda x: x[1].reference_start)
    
    # Extract info from each hit
    bam_indices = []
    cigars = []
    starts = []
    ends = []
    
    for bam_index, hit in sorted_indexed_hits:
        bam_indices.append(bam_index)
        cigars.append(Cigar(hit.cigartuples))
        starts.append(hit.reference_start)
        ends.append(hit.reference_end)
    
    # Overall span
    combined_start = min(starts)
    combined_end = max(ends)
    
    return CombinedSingleAlignment(
        start=combined_start,
        end=combined_end,
        bam_indices=tuple(bam_indices),
        cigars=tuple(cigars),
        starts=tuple(starts),
        ends=tuple(ends),
    )

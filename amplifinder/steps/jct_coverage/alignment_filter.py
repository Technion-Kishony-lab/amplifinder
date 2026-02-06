import pysam
from typing import Optional

from amplifinder.config import AlignmentFilterParams
from amplifinder.steps.jct_coverage.cigar import Cigar


# Only allow M/I/D/=/X operations
ALLOWED_CIGAR_OPERATIONS = {0, 1, 2, 7, 8}


def bam_hit_passes_filters(
    hit: pysam.AlignedSegment,
    avg_read_length: int,
    alignment_filter_params: AlignmentFilterParams,
    arm_len: int,
    allow_indels_at_junction_distance: Optional[int]
) -> bool:
    """Filter BAM hit based on CIGAR validation, flags, alignment params, and indels check.

    Returns:
        True if all filters pass, False otherwise
    """
    # Create and validate CIGAR
    cigar = Cigar(hit.cigartuples or [])
    if not cigar or not cigar.has_only_operations(ALLOWED_CIGAR_OPERATIONS):
        return False

    # Skip unmapped, supplementary
    if hit.is_unmapped or hit.is_supplementary:
        return False

    # Check alignment length
    read_length_factor = 1 + alignment_filter_params.filter_len_tolerance
    min_alignment_length = avg_read_length / read_length_factor
    max_alignment_length = avg_read_length * read_length_factor
    alignment_length = hit.reference_end - hit.reference_start
    if not (min_alignment_length <= alignment_length <= max_alignment_length):
        return False

    # Check alignment score filters
    if alignment_filter_params.min_as_score is not None:
        as_score = hit.get_tag("AS")
        if as_score < alignment_filter_params.min_as_score:
            return False

    # Check NM score if specified
    if alignment_filter_params.max_nm_score is not None:
        nm_score = hit.get_tag("NM")
        if nm_score > alignment_filter_params.max_nm_score:
            return False

    # Check indels within limit
    return _indels_within_limit(
        hit.reference_start, cigar, arm_len, allow_indels_at_junction_distance
    )


def _indels_within_limit(
    reference_start: int,
    cigar: Cigar,
    arm_len: int,
    allow_indels_at_junction_distance: Optional[int]
) -> bool:
    # Check for indels (insertions=1 or deletions=2)
    has_indels = any(op in (1, 2) for op, _ in cigar)

    if not has_indels:
        # Only M/=/X operations - acceptable
        return True

    # Has indels - check if they're close enough to junction
    if allow_indels_at_junction_distance is None:
        # No distance limit - allow all indels
        return True

    if allow_indels_at_junction_distance == -1:
        # No indels allowed
        return False

    # Check distance of each indel from junction
    junction_pos = arm_len  # 0-based junction position

    for op, length, ref_pos in cigar.iter_with_ref_pos(reference_start):
        if op == 1:  # Insertion - occurs at a single reference position
            dist_to_junction = abs(ref_pos - junction_pos)  # 0 when precisely at junction
            if dist_to_junction > allow_indels_at_junction_distance:
                return False
        elif op == 2:  # Deletion - spans from ref_pos to ref_pos + length - 1
            dist_start = abs(ref_pos - junction_pos)
            dist_end = abs(ref_pos + length - junction_pos)
            if dist_start > allow_indels_at_junction_distance or dist_end > allow_indels_at_junction_distance:
                return False

    return True

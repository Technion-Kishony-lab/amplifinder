from pathlib import Path
from typing import Optional, Union

from pandas.core.base import np
import pysam

from amplifinder.config import AlignmentFilterParams, AlignmentClassifyParams
from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.data_types.enums import ReadType
from amplifinder.steps.jct_coverage.alignment_data import AlignmentData, SingleAlignment, PairedAlignment
from amplifinder.steps.jct_coverage.hits_container import HitsContainer
from amplifinder.env import IGNORE_DUPLICATES, LINK_PAIRED_END, ALLOW_INDELS_AT_JUNCTION_DISTANCE


if LINK_PAIRED_END:
    assert IGNORE_DUPLICATES, "LINK_PAIRED_END requires IGNORE_DUPLICATES"


def get_jct_read_counts(
    bam_path: Path,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
    alignment_filter_params: AlignmentFilterParams,
) -> tuple[
    dict[JunctionType, JunctionReadCounts],
    dict[JunctionType, list[AlignmentData]],
    dict[JunctionType, int]
]:
    """Parse BAM and get coverage for all 7 junction types.

    Args:
        bam_path: Path to BAM file
        avg_read_length: Read length for filtering
        align_params: Alignment analysis parameters
        alignment_filter_params: Read filtering parameters

    Returns:
        tuple of:
            - dict mapping JunctionType -> JunctionReadCounts
            - dict mapping JunctionType -> list of SingleAlignment or PairedAlignment tuples
            - dict mapping JunctionType -> junction length
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    jctypes_to_readtypes_to_readids_to_hits, jct_types_to_lengths = _classify_bam_hits_by_jc_types_and_read_types(
        bam_path, avg_read_length, alignment_classify_params, alignment_filter_params
    )
    
    if LINK_PAIRED_END:
        _merge_hits_of_paired_end_reads(jctypes_to_readtypes_to_readids_to_hits)
    
    # Count and collect hits by junction type and read type
    counts = {jt: JunctionReadCounts() for jt in JunctionType}
    alignment_data = {jt: [] for jt in JunctionType}
    
    for jct_type in JunctionType:
        readtypes_to_readids_to_hits = jctypes_to_readtypes_to_readids_to_hits[jct_type]
        for read_type in ReadType:
            hits_container = readtypes_to_readids_to_hits[read_type]
            counts[jct_type][read_type] = len(hits_container)
            alignment_data[jct_type].extend(hits_container)
    
    return counts, alignment_data, jct_types_to_lengths


def _classify_bam_hits_by_jc_types_and_read_types(
    bam_path: Path,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
    alignment_filter_params: AlignmentFilterParams,
) -> tuple[dict[JunctionType, dict[ReadType, HitsContainer]], dict[JunctionType, int]]:
    """Parse BAM file and classify hits by junction and read type."""
    jctypes_to_readtypes_to_hits_container: dict[JunctionType, dict[ReadType, HitsContainer]] = {
        jt: {rt: HitsContainer(IGNORE_DUPLICATES) for rt in ReadType} for jt in JunctionType
    }
    
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        jct_names_to_lengths = dict(zip(bam.references, bam.lengths))
        jct_types_to_lengths = {JunctionType[jct_name]: length for jct_name, length in jct_names_to_lengths.items()}

        for bam_index, hit in enumerate(bam.fetch(), start=1):
            jct_type = JunctionType[hit.reference_name]
            arm_len = jct_types_to_lengths[jct_type] // 2

            read_type = analyze_hit(hit, arm_len, avg_read_length, alignment_filter_params, alignment_classify_params)
            if read_type is None:
                continue

            read_id = hit.query_name
            hits_container = jctypes_to_readtypes_to_hits_container[jct_type][read_type]
            alignment = SingleAlignment(read_type, hit.reference_start, hit.reference_end, bam_index)
            hits_container.add_hit(alignment, read_id)

    return jctypes_to_readtypes_to_hits_container, jct_types_to_lengths


def _merge_hits_of_paired_end_reads(
    jctypes_to_readtypes_to_readids_to_hits: dict[JunctionType, dict[ReadType, HitsContainer]]
) -> None:
    """
    Detect paired-end reads (ids appearing in both LEFT and RIGHT) 
    and merge them into PAIRED category.
    """
    for jct_type in JunctionType:
        readtypes_to_readids_to_hits = jctypes_to_readtypes_to_readids_to_hits[jct_type]
        
        readids_to_hits_L = readtypes_to_readids_to_hits[ReadType.LEFT]
        readids_to_hits_R = readtypes_to_readids_to_hits[ReadType.RIGHT]
        paired_read_ids = set(readids_to_hits_L.keys()) & set(readids_to_hits_R.keys())

        paired_readids_to_hits = readtypes_to_readids_to_hits[ReadType.PAIRED]
        assert len(paired_readids_to_hits) == 0
        
        for read_id in paired_read_ids:
            hit_L = readids_to_hits_L.pop(read_id)
            hit_R = readids_to_hits_R.pop(read_id)
            paired_readids_to_hits[read_id] = PairedAlignment(ReadType.PAIRED,
                hit_L.start, hit_L.end, hit_L.bam_index,
                hit_R.start, hit_R.end, hit_R.bam_index,
            )


def analyze_hit(
    hit: pysam.AlignedSegment,
    arm_len: int,
    avg_read_length: int,
    alignment_filter_params: AlignmentFilterParams,
    alignment_classify_params: AlignmentClassifyParams,
) -> Optional[ReadType]:

    # Filter by quality
    if not passes_alignment_filters(hit, alignment_filter_params, avg_read_length):
        return None

    # Filter by CIGAR string
    if not is_hit_cigar_acceptable(hit, arm_len):
        return None

    # Get read type
    start = hit.reference_start
    end = hit.reference_end
    start, end = start + 1, end  # convert to 1-based coordinates, end-inclusive

    return get_hit_type(start, end, arm_len, alignment_classify_params)


def passes_alignment_filters(hit: pysam.AlignedSegment, alignment_filter_params: AlignmentFilterParams, avg_read_length: int) -> bool:

    # Filter by alignment length
    read_length_factor = 1 + alignment_filter_params.length_tolerance
    min_alignment_length = avg_read_length / read_length_factor
    max_alignment_length = avg_read_length * read_length_factor

    if not (min_alignment_length <= hit.query_alignment_length <= max_alignment_length):
        return False

    # Filter by alignment flags
    if hit.is_unmapped or hit.is_secondary or hit.is_supplementary:
        return False

    # Filter by alignment tags
    nm_score = hit.get_tag("NM")
    if alignment_filter_params.max_nm_score is not None and nm_score > alignment_filter_params.max_nm_score:
        return False

    as_score = hit.get_tag("AS")
    if alignment_filter_params.min_as_score is not None and as_score < alignment_filter_params.min_as_score:
        return False

    return True


def get_hit_type(start: int, end: int, arm_len: int, alignment_classify_params: AlignmentClassifyParams) -> Optional[ReadType]:
    """Determine the read type based on the start and end positions."""
    # n = arm_len
    # M = min_overlap_len
    # L = avg_read_len
    # F = min_bp_in_frame

    # index right-arm:      <--neg-----0-----pos-->
    #                                  |
    #                                 n+1   n+1+M                    2n
    #                                  |      |                       |
    # +------------------------+------++------+-----------------------+
    # |                        |      |
    # 1                       n-M     n
    #                                 |
    # index left-arm:      <--pos-----0-----neg-->

    idx_L_1 = arm_len - end
    idx_L_2 = arm_len - start
    assert idx_L_2 >= idx_L_1

    idx_R_1 = start - (arm_len + 1)
    idx_R_2 = end - (arm_len + 1)
    assert idx_R_2 >= idx_R_1

    # remove reads that are too far from junction
    if idx_L_1 >= alignment_classify_params.max_dist_from_junction:
        return None
    if idx_R_1 >= alignment_classify_params.max_dist_from_junction:
        return None

    # classify reads
    if idx_L_1 >= 0:
        return ReadType.LEFT
    elif idx_R_1 >= 0:
        return ReadType.RIGHT
    is_start_left_of_margin = idx_L_2 >= alignment_classify_params.min_overlap_len - 1
    is_end_right_of_margin = idx_R_2 >= alignment_classify_params.min_overlap_len - 1
    if is_start_left_of_margin and is_end_right_of_margin:
        return ReadType.MIDDLE
    elif is_start_left_of_margin:
        return ReadType.LEFT_MARGINAL
    elif is_end_right_of_margin:
        return ReadType.RIGHT_MARGINAL
    assert False


def is_hit_cigar_acceptable(hit: pysam.AlignedSegment, arm_len: int) -> bool:
    # CIGAR operation codes:
    # 0=M (match/mismatch)
    # 1=I (insertion)
    # 2=D (deletion)
    # 3=N (skipped region)
    # 4=S (soft clip)
    # 5=H (hard clip)
    # 6=P (padding)
    # 7= (sequence match)
    # 8=X (sequence mismatch)

    if not hit.cigartuples:
        return False

    # Check for operations other than M (0), indels (1, 2), and match/mismatch (7, 8)
    if any(op not in (0, 1, 2, 7, 8) for op, _ in hit.cigartuples):
        return False
    
    # Check for indels (insertions=1 or deletions=2)
    has_indels = any(op in (1, 2) for op, _ in hit.cigartuples)
    
    if not has_indels:
        # Only M operations - acceptable
        return True
    
    # Has indels - check if they're close enough to junction
    if ALLOW_INDELS_AT_JUNCTION_DISTANCE is None:
        # No distance limit - allow all indels
        return True
    
    if ALLOW_INDELS_AT_JUNCTION_DISTANCE == -1:
        # No indels allowed
        return False
    
    # Check distance of each indel from junction
    junction_pos = arm_len  # 0-based junction position
    ref_pos = hit.reference_start  # Current position on reference (0-based)
    
    for op, length in hit.cigartuples:
        if op == 1:  # Insertion - occurs at a single reference position
            dist_to_junction = abs(ref_pos - junction_pos)  # 0 when precisely at junction
            if dist_to_junction > ALLOW_INDELS_AT_JUNCTION_DISTANCE:
                return False
        elif op == 2:  # Deletion - spans from ref_pos to ref_pos + length - 1
            dist_start = abs(ref_pos - junction_pos)
            dist_end = abs(ref_pos + length - junction_pos)
            if dist_start > ALLOW_INDELS_AT_JUNCTION_DISTANCE or dist_end > ALLOW_INDELS_AT_JUNCTION_DISTANCE:
                return False
        
        # Update reference position (insertions don't consume reference)
        if op in (0, 2, 7, 8):  # M, D, =, X consume reference
            ref_pos += length
    
    return True


def get_expected_counts(avg_read_len: int, arm_len: int, alignment_classify_params: AlignmentClassifyParams) -> JunctionReadCounts:
    """Return the expected counts assuming uniform distribution of reads across the junction."""
    min_overlap_len = alignment_classify_params.min_overlap_len
    max_dist = alignment_classify_params.max_dist_from_junction
    if max_dist is None:
        max_dist = np.inf
    max_dist = min(max_dist, arm_len - avg_read_len)
    return JunctionReadCounts(
        left=max_dist + 1,
        left_marginal=min_overlap_len,
        spanning=avg_read_len - 2 * min_overlap_len - 1,
        right_marginal=min_overlap_len,
        right=max_dist + 1
    )

from pathlib import Path
from typing import Optional

import pysam

from amplifinder.config import AlignmentFilterParams, AlignmentClassifyParams
from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.data_types.enums import ReadType


def get_jct_read_counts(
    bam_path: Path,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
    alignment_filter_params: AlignmentFilterParams,
) -> tuple[
    dict[JunctionType, JunctionReadCounts],
    dict[JunctionType, list[tuple[int, int, str]]],
    dict[JunctionType, int],
    dict[JunctionType, dict[ReadType, list[tuple[str, str]]]]
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
            - dict mapping JunctionType -> list of (start, end, read_type) tuples
            - dict mapping JunctionType -> junction length
            - dict mapping JunctionType -> dict[ReadType -> list of (read_id, seq, bam_index) tuples]
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    counts = {jt: JunctionReadCounts() for jt in JunctionType}
    alignment_data = {jt: [] for jt in JunctionType}

    jcs_to_readtypes_to_hits = {
        jt: {rt: [] for rt in ReadType}
        for jt in JunctionType
    }

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:

        # Store lengths by JunctionType
        jct_names_to_lengths = dict(zip(bam.references, bam.lengths))
        jct_types_to_lengths = {JunctionType[jct_name]: length for jct_name, length in jct_names_to_lengths.items()}

        # Count read-alignments of each type at each junction
        # Track 1-based BAM index to match MATLAB format
        for bam_index, hit in enumerate(bam.fetch(), start=1):
            jct_name = hit.reference_name
            jct_type = JunctionType[jct_name]
            arm_len = jct_types_to_lengths[jct_type] // 2

            read_type = analyze_hit(
                hit, arm_len, avg_read_length,
                alignment_filter_params, alignment_classify_params
            )

            if read_type is None:
                continue

            counts[jct_type].increment(read_type)
            alignment_data[jct_type].append((hit.reference_start, hit.reference_end, read_type))

            read_id = hit.query_name
            seq = hit.query_sequence or ""
            jcs_to_readtypes_to_hits[jct_type][read_type].append((read_id, seq, bam_index))

    return counts, alignment_data, jct_types_to_lengths, jcs_to_readtypes_to_hits


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
    if not is_hit_cigar_acceptable(hit):
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
    if idx_L_1 > alignment_classify_params.max_dist_from_junction:
        return None
    if idx_R_1 > alignment_classify_params.max_dist_from_junction:
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


def is_hit_cigar_acceptable(hit: pysam.AlignedSegment) -> bool:
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

    # MATLAB RdFull: only CIGAR-only-M reads are counted
    if any(op != 0 for op, _ in hit.cigartuples):
        return False

    return True


def get_expected_counts(avg_read_len: int, arm_len: int, alignment_classify_params: AlignmentClassifyParams) -> JunctionReadCounts:
    """Return the expected counts assuming uniform distribution of reads across the junction."""
    arm_len = avg_read_len + alignment_classify_params.max_dist_from_junction
    return JunctionReadCounts(
        left=arm_len - avg_read_len + 1,
        left_marginal=alignment_classify_params.min_overlap_len,
        spanning=avg_read_len - 2 * alignment_classify_params.min_overlap_len - 1,
        right_marginal=alignment_classify_params.min_overlap_len,
        right=arm_len - avg_read_len + 1
    )

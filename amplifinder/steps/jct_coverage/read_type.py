from pathlib import Path
import pysam

from typing import Optional

from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.data_types.enums import ReadType


def get_jct_read_counts(
    bam_path: Path,
    avg_read_length: int,
    min_overlap_len: int,
    alignment_length_tolerance: float = 0.1,
    max_dist_from_junction: int = 10,
    max_nm_score: int = 3,
    min_as_score: int = -25
) -> tuple[
    dict[JunctionType, JunctionReadCounts],
    dict[JunctionType, list[tuple[int, int, str]]],
    dict[JunctionType, int]
]:
    """Parse BAM and get coverage for all 7 junction types.

    Args:
        bam_path: Path to BAM file
        avg_read_length: Read length for filtering
        min_overlap_len: Minimum overlap length required
        alignment_length_tolerance: Tolerance for read length filtering (default 0.1 = 10%)
        max_dist_from_junction: Maximum distance from junction for read classification
        max_nm_score: Maximum NM score threshold (default 3)
        min_as_score: Minimum AS score threshold (default -25)

    Returns:
        tuple of:
            - dict mapping JunctionType -> JunctionReadCounts
            - dict mapping JunctionType -> list of (start, end, read_type) tuples
            - dict mapping JunctionType -> junction length
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    read_length_factor = 1 + alignment_length_tolerance
    min_alignment_length = avg_read_length / read_length_factor
    max_alignment_length = avg_read_length * read_length_factor

    counts = {jt: JunctionReadCounts() for jt in JunctionType}
    alignment_data = {jt: [] for jt in JunctionType}

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:

        # Store lengths by JunctionType
        jct_names_to_lengths = dict(zip(bam.references, bam.lengths))
        jct_types_to_lengths = {JunctionType[jct_name]: length for jct_name, length in jct_names_to_lengths.items()}

        # Count read-alignments of each type at each junction
        for read in bam.fetch():
            jct_name = read.reference_name
            jct_type = JunctionType[jct_name]
            arm_len = jct_types_to_lengths[jct_type] // 2

            read_type = classify_read(read, arm_len, min_alignment_length, max_alignment_length,
                                      max_nm_score, min_as_score, min_overlap_len, max_dist_from_junction)

            if read_type is not None:
                counts[jct_type].increment(read_type)
                alignment_data[jct_type].append((read.reference_start, read.reference_end, read_type))

    return counts, alignment_data, jct_types_to_lengths


def classify_read(read: pysam.AlignedSegment, arm_len: int, min_alignment_length: int,
                  max_alignment_length: int, max_nm_score: int, min_as_score: int,
                  min_overlap_len: int, max_dist_from_junction: int) -> Optional[ReadType]:

    # Filter by quality
    if not passes_quality_filters(read, max_nm_score, min_as_score):
        return None

    # Filter by alignment length
    alignment_length = read.query_alignment_length
    if not (min_alignment_length <= alignment_length <= max_alignment_length):
        return None

    # Filter by CIGAR string
    if not is_read_cigar_acceptable(read):
        return None

    # Get read type
    start = read.reference_start
    end = read.reference_end
    start, end = start + 1, end  # convert to 1-based coordinates, end-inclusive

    return get_read_type(start, end, arm_len,
                         min_overlap_len=min_overlap_len,
                         max_dist_from_junction=max_dist_from_junction)


def passes_quality_filters(read: pysam.AlignedSegment, max_nm_score: int, min_as_score: int) -> bool:
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return False
    try:
        as_score = read.get_tag("AS")
    except KeyError:
        as_score = None
    try:
        nm_score = read.get_tag("NM")
    except KeyError:
        try:
            nm_score = read.get_tag("nM")
        except KeyError:
            nm_score = None

    if nm_score is not None and nm_score > max_nm_score:
        return False
    if as_score is not None and as_score < min_as_score:
        return False
    return True


def get_read_type(start: int, end: int, arm_len: int, min_overlap_len: int,
                  max_dist_from_junction: int) -> Optional[ReadType]:
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
    if idx_L_1 > max_dist_from_junction:
        return None
    if idx_R_1 > max_dist_from_junction:
        return None

    # classify reads
    if idx_L_1 >= 0:
        return ReadType.LEFT
    elif idx_R_1 >= 0:
        return ReadType.RIGHT
    is_start_left_of_margin = idx_L_2 >= min_overlap_len - 1
    is_end_right_of_margin = idx_R_2 >= min_overlap_len - 1
    if is_start_left_of_margin and is_end_right_of_margin:
        return ReadType.MIDDLE
    elif is_start_left_of_margin:
        return ReadType.LEFT_MARGINAL
    elif is_end_right_of_margin:
        return ReadType.RIGHT_MARGINAL
    assert False


def is_read_cigar_acceptable(read: pysam.AlignedSegment) -> bool:
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

    if not read.cigartuples:
        return False

    # MATLAB RdFull: only CIGAR-only-M reads are counted
    if any(op != 0 for op, _ in read.cigartuples):
        return False

    return True


def get_expected_counts(arm_len: int, min_overlap_len: int, max_dist_from_junction, read_len: int) -> JunctionReadCounts:
    """Return the expected counts assuming uniform distribution of reads across the junction."""
    arm_len = read_len + max_dist_from_junction
    return JunctionReadCounts(
        left=arm_len - read_len + 1,
        left_marginal=min_overlap_len,
        spanning=read_len - 2 * min_overlap_len - 1,
        right_marginal=min_overlap_len,
        right=arm_len - read_len + 1
    )

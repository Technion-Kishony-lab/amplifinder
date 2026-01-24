"""Classify read types from alignments."""
from pathlib import Path
from typing import Optional
from functools import partial

from amplifinder.config import AlignmentFilterParams, AlignmentClassifyParams
from amplifinder.data_types import JunctionReadCounts, JunctionType, Side, ReadType
from amplifinder.steps.jct_coverage.alignment_data import (
    AlignmentData, BaseSingleAlignment, PairedAlignment
)
from amplifinder.steps.jct_coverage.read_bam import \
    read_bam_and_group_single_alignments
from amplifinder.steps.jct_coverage.combine_hits import (
    flatten_combined_alignments,
    select_or_combine_single_alignments,
    combine_same_id_different_orientation_hits,
)
from amplifinder.steps.jct_coverage.alignment_filter import bam_hit_passes_filters
from amplifinder.env import LINK_PAIRED_END, SELECT_BEST_BY_SCORE, ALLOW_INDELS_AT_JUNCTION_DISTANCE
from amplifinder.env import DEBUG


def get_jct_read_counts(
    bam_path: Path,
    arm_len: int,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
    alignment_filter_params: AlignmentFilterParams,
) -> tuple[
    dict[JunctionType, JunctionReadCounts],
    dict[JunctionType, list[AlignmentData]]
]:
    """Parse BAM and get coverage for all 7 junction types.

    Args:
        bam_path: Path to BAM file
        arm_len: Length of each junction arm
        avg_read_length: Read length for filtering
        alignment_classify_params: Classification parameters
        alignment_filter_params: Filtering parameters

    Returns:
        tuple of:
            - dict mapping JunctionType -> JunctionReadCounts
            - dict mapping JunctionType -> list[AlignmentData]
    """
    # Read BAM and group hits by read_id and orientation, with filtering
    filter_func = partial(
        bam_hit_passes_filters,
        avg_read_length=avg_read_length,
        alignment_filter_params=alignment_filter_params,
        arm_len=arm_len,
        allow_indels_at_junction_distance=ALLOW_INDELS_AT_JUNCTION_DISTANCE
    )

    refnames_to_readids_to_orientations_to_alignments = read_bam_and_group_single_alignments(
        bam_path,
        filter_func=filter_func
    )

    # Select best or combine hits with same read_id and orientation (in place)
    select_or_combine_single_alignments(
        refnames_to_readids_to_orientations_to_alignments, select_best_by_score=SELECT_BEST_BY_SCORE
    )

    # Combine hits with same read_id but different orientations (paired-end) in place
    if LINK_PAIRED_END:
        combine_same_id_different_orientation_hits(refnames_to_readids_to_orientations_to_alignments)

    # Flatten combined alignments and classify by read type
    refnames_to_alignments = flatten_combined_alignments(refnames_to_readids_to_orientations_to_alignments)

    # Classify alignments by read type (also counts as it goes)
    return _classify_alignments_by_read_types(
        refnames_to_alignments, arm_len, avg_read_length, alignment_classify_params
    )


def _classify_alignments_by_read_types(
    refnames_to_alignments: dict[str, list[AlignmentData]],
    arm_len: int,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
) -> tuple[dict[JunctionType, JunctionReadCounts], dict[JunctionType, list[AlignmentData]]]:
    """Classify alignments by junction and read type, returning counts + flat structure."""
    jctypes_to_alignments: dict[JunctionType, list[AlignmentData]] = {jt: [] for jt in JunctionType}
    counts = {jt: JunctionReadCounts() for jt in JunctionType}

    invalid_side_combination_count = 0

    for ref_name, alignments in refnames_to_alignments.items():
        jct_type = JunctionType[ref_name]

        for alignment in alignments:
            classified_alignments = _classify_alignment(
                alignment, arm_len, avg_read_length, alignment_classify_params
            )
            if not classified_alignments:
                invalid_side_combination_count += 1
            for alignment in classified_alignments:
                jctypes_to_alignments[jct_type].append(alignment)
                counts[jct_type].increment(alignment.read_type)

    if invalid_side_combination_count > 0:
        print(f"Total invalid side combinations: {invalid_side_combination_count}", flush=True)

    return counts, jctypes_to_alignments


def _classify_alignment(
    alignment: AlignmentData,
    arm_len: int,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
) -> list[AlignmentData]:
    """Classify alignment by read type.

    Args:
        alignment: BaseSingleAlignment or PairedAlignment
        arm_len: Half the junction length
        avg_read_length: Average read length (unused, kept for consistency)
        alignment_classify_params: Classification parameters

    Returns:
        list[AlignmentData] (empty list if invalid side combination)
    """
    max_debug_examples = 1

    # Check indels within limit for all nested single alignments
    if isinstance(alignment, BaseSingleAlignment):
        classify_single_alignment(alignment, arm_len, alignment_classify_params)
        return [alignment]

    assert isinstance(alignment, PairedAlignment)

    # Classify forward alignment
    fwd_alignment = alignment.forward_alignment
    classify_single_alignment(fwd_alignment, arm_len, alignment_classify_params)
    fwd_read_type = fwd_alignment.read_type
    fwd_side = fwd_read_type.get_side()

    # Classify reverse alignment
    rev_alignment = alignment.reverse_alignment
    classify_single_alignment(rev_alignment, arm_len, alignment_classify_params)
    rev_read_type = rev_alignment.read_type
    rev_side = rev_read_type.get_side()

    # Classify paired alignment
    if fwd_read_type == rev_read_type:
        alignment.read_type = fwd_read_type
        return [alignment]
    if fwd_side == rev_side:
        if fwd_side == Side.LEFT:
            return [rev_alignment]
        if fwd_side == Side.RIGHT:
            return [fwd_alignment]
        assert False, (
            "fwd_side == rev_side == Side.MIDDLE, "
            "so fwd_read_type == rev_read_type == ReadType.SPANNING"
        )
    if fwd_side == Side.LEFT and rev_side == Side.RIGHT:
        alignment.read_type = ReadType.PAIRED
        # We have evience for left, right and spanning
        return [fwd_alignment, rev_alignment, alignment]

    ok = (fwd_side == Side.LEFT and rev_side == Side.MIDDLE) \
        or (fwd_side == Side.MIDDLE and rev_side == Side.RIGHT)
    if not ok:
        if DEBUG and max_debug_examples > 0:
            max_debug_examples -= 1
            print(f"Example of invalid side combination: fwd_side={fwd_side}, rev_side={rev_side}", flush=True)
            print(alignment, flush=True)
        return []
    return [fwd_alignment, rev_alignment]


def classify_single_alignment(
    alignment: BaseSingleAlignment,
    arm_len: int,
    alignment_classify_params: AlignmentClassifyParams,
):
    """Classify single alignment by read type."""
    start = alignment.start + 1  # convert to 1-based coordinates
    end = alignment.end
    read_type = get_hit_type(start, end, arm_len, alignment_classify_params)
    alignment.read_type = read_type


def get_hit_type(
    start: int, end: int, arm_len: int,
    alignment_classify_params: AlignmentClassifyParams
) -> Optional[ReadType]:
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
    max_dist = alignment_classify_params.get_max_dist_from_junction(arm_len)
    min_overlap_len = alignment_classify_params.min_overlap_len

    # classify reads
    if idx_L_1 >= 0:
        return ReadType.LEFT if idx_L_2 < max_dist else ReadType.LEFT_FAR
    elif idx_R_1 >= 0:
        return ReadType.RIGHT if idx_R_2 < max_dist else ReadType.RIGHT_FAR
    is_start_left_of_margin = idx_L_2 >= min_overlap_len - 1
    is_end_right_of_margin = idx_R_2 >= min_overlap_len - 1
    if is_start_left_of_margin and is_end_right_of_margin:
        return ReadType.SPANNING
    elif is_start_left_of_margin:
        return ReadType.LEFT_MARGINAL
    elif is_end_right_of_margin:
        return ReadType.RIGHT_MARGINAL
    assert False


def get_expected_counts(
    avg_read_len: int, arm_len: int,
    alignment_classify_params: AlignmentClassifyParams
) -> JunctionReadCounts:
    """Return the expected counts assuming uniform distribution of reads across the junction."""
    min_overlap_len = alignment_classify_params.min_overlap_len
    max_dist = alignment_classify_params.get_max_dist_from_junction(arm_len)
    counts = JunctionReadCounts(
        left_far=arm_len - max_dist,
        left=max_dist - avg_read_len + 1,
        left_marginal=min_overlap_len,
        spanning=avg_read_len - 2 * min_overlap_len - 1,
        paired=0,  # Not part of expected distribution
        right_marginal=min_overlap_len,
        right=max_dist - avg_read_len + 1,
        right_far=arm_len - max_dist,
    )
    assert counts.total == 2 * arm_len - avg_read_len + 1
    return counts

"""Combine aligned segments."""
from collections import defaultdict

from amplifinder.logger import logger
from amplifinder.steps.jct_coverage.alignment_data import \
    AlignmentData, BaseSingleAlignment, CombinedSingleAlignment, SingleAlignment, PairedAlignment
from amplifinder.utils.file_utils import fmt_count


def select_or_combine_single_alignments(
    grouped_hits: dict[str, dict[str, dict[bool, list[SingleAlignment]]]],
    select_best_by_score: bool = True,
) -> None:
    """Process grouped hits in place: select best by score or combine all.

    Modifies grouped_hits:
    - If select_best_by_score=True: replaces list with best scoring SingleAlignment
    - If select_best_by_score=False: replaces list with CombinedSingleAlignment

    Args:
        grouped_hits: Dictionary to modify in place
        select_best_by_score: If True, keep best scoring hit. If False, combine all hits.
    """
    total_alignments = defaultdict(int)
    single_alignments = defaultdict(int)

    for ref_name, ref_reads in grouped_hits.items():
        for read_id, orientations in ref_reads.items():
            for is_reverse, alignments in orientations.items():
                total_alignments[is_reverse] += 1
                if len(alignments) == 1:
                    # There is only one alignment for this read_id and orientation, so we can use it directly
                    alignment = alignments[0]
                    single_alignments[is_reverse] += 1
                elif select_best_by_score:
                    # Take the alignment with the best (highest) score
                    alignment = max(alignments, key=lambda a: a.alignment_score)
                else:
                    # Create a combined alignment from all alignments
                    alignment = CombinedSingleAlignment.from_alignments(alignments)
                orientations[is_reverse] = alignment

    action = "SELECT_BEST" if select_best_by_score else "COMBINE_ALL"
    logger.info(
        f"More than one alignment per read (action: {action:12s}) " + " "*7 +
        f"FWD={fmt_count(total_alignments[False] - single_alignments[False], total_alignments[False])}, "
        f"REV={fmt_count(total_alignments[True] - single_alignments[True], total_alignments[True])}",
        timestamp=False,
    )


def combine_same_id_different_orientation_hits(
    combined_alignments: dict[str, dict[str, dict[bool, BaseSingleAlignment]]]
) -> None:
    """Combine alignments with same read_id but different orientations into PairedAlignment.

    Modifies in place:
    - Paired reads: removes is_reverse=True/False, adds is_reverse=None -> PairedAlignment
    - Unpaired reads: unchanged (is_reverse=True/False -> BaseSingleAlignment)

    Args:
        combined_alignments: Dictionary to modify in place
    """
    total_reads = 0
    total_normal_pairs = 0
    total_swapped_pairs = 0
    sum_distances_normal = 0
    sum_distances_swapped = 0

    for ref_name, ref_reads in combined_alignments.items():
        for read_id, orientations in ref_reads.items():
            total_reads += 1
            if len(orientations) == 2:  # Has both forward and reverse - create paired
                fwd = orientations.pop(False)
                rev = orientations.pop(True)
                if fwd.left < rev.right:
                    total_normal_pairs += 1
                    paired = PairedAlignment(
                        forward_alignment=fwd,
                        reverse_alignment=rev,
                        is_swapped=False
                    )
                    sum_distances_normal += paired.overlapping_length
                else:
                    total_swapped_pairs += 1
                    paired = PairedAlignment(
                        forward_alignment=rev,
                        reverse_alignment=fwd,
                        is_swapped=True
                    )
                    sum_distances_swapped += paired.overlapping_length
                    logger.debug_message(
                        f"Example of a swapped paired alignment:\n{paired}",
                        category="swapped_paired_alignment",
                        max_prints=1
                    )
                orientations[None] = paired

    total_pairs = total_normal_pairs + total_swapped_pairs
    avg_distance_normal = sum_distances_normal / total_normal_pairs if total_normal_pairs > 0 else 0
    avg_distance_swapped = sum_distances_swapped / total_swapped_pairs if total_swapped_pairs > 0 else 0
    total_singletons = total_reads - total_pairs
    total_hits = total_singletons + total_pairs * 2
    logger.info(f"Total hits: {total_hits:6} = ({total_singletons:6} Singletons) + 2 * ({total_pairs:6} Pairs)\n"
          f"\tNormal  pairs: {total_normal_pairs:6}, Avg overlap (bp): {avg_distance_normal:6.1f}\n"
          f"\tSwapped pairs: {total_swapped_pairs:6}, Avg overlap (bp): {avg_distance_swapped:6.1f}",
          timestamp=False,
    )


def flatten_combined_alignments(
    combined_alignments: dict[str, dict[str, dict[bool | None, BaseSingleAlignment | PairedAlignment]]]
) -> dict[str, list[AlignmentData]]:
    """Flatten alignments into reference_name -> list[AlignmentData].

    Args:
        combined_alignments: reference_name -> read_id -> is_reverse/None -> AlignmentData

    Returns:
        reference_name -> list[AlignmentData]
    """
    flattened_alignments: dict[str, list[AlignmentData]] = defaultdict(list)
    for ref_name, ref_reads in combined_alignments.items():
        for read_id, orientations in ref_reads.items():
            for is_reverse, alignment in orientations.items():
                flattened_alignments[ref_name].append(alignment)
    return flattened_alignments

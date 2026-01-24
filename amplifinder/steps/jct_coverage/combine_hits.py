from pathlib import Path
from collections import defaultdict

from amplifinder.steps.jct_coverage.alignment_data import \
    AlignmentData, BaseSingleAlignment, SingleAlignment, PairedAlignment, \
    create_single_or_combined_alignment


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
    for ref_name, ref_reads in grouped_hits.items():
        for read_id, orientations in ref_reads.items():
            for is_reverse, alignments in orientations.items():
                if select_best_by_score:
                    # Take the alignment with the best (highest) score
                    orientations[is_reverse] = max(alignments, key=lambda a: a.alignment_score)
                else:
                    # Combine all alignments
                    orientations[is_reverse] = create_single_or_combined_alignment(alignments)


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
    for ref_name, ref_reads in combined_alignments.items():
        for read_id, orientations in ref_reads.items():
            if len(orientations) == 2:  # Has both forward and reverse - create paired
                fwd = orientations.pop(False)
                rev = orientations.pop(True)
                assert fwd.left < rev.right
                orientations[None] = PairedAlignment(
                    forward_alignment=fwd,
                    reverse_alignment=rev
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

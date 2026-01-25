"""Read BAM file and group alignments."""
import pysam
from pathlib import Path
from collections import defaultdict
from typing import Optional, Callable

from amplifinder.steps.jct_coverage.alignment_data import SingleAlignment
from amplifinder.steps.jct_coverage.cigar import Cigar, resolve_cigar_m_operations
from amplifinder.env import DEBUG, RESOLVE_CIGAR_MATCHES_VS_MISMATCHES
from amplifinder.utils.file_utils import fmt_count


def read_bam_and_group_single_alignments(
    bam_path: Path,
    filter_func: Optional[Callable[[pysam.AlignedSegment], bool]] = None,
) -> dict[str, dict[str, dict[bool, list[SingleAlignment]]]]:
    """Read BAM file into SingleAlignment objects grouped by reference, read_id, and orientation.

    Args:
        bam_path: Path to BAM file
        filter_func: Optional function to filter BAM hits (returns True to keep)

    Returns:
        Dictionary mapping reference_name -> read_id -> is_reverse -> list[SingleAlignment]
    """
    grouped_hits: dict[str, dict[str, dict[bool, list[SingleAlignment]]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(list))
    )

    # Statistics counters
    total_hits = 0
    pass_filter_total = 0
    pass_filter_forward = 0
    pass_filter_reverse = 0

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for bam_index, hit in enumerate(bam.fetch(), start=1):
            total_hits += 1
            
            # Apply filter if provided
            if filter_func is not None and not filter_func(hit):
                continue

            pass_filter_total += 1
            if hit.is_reverse:
                pass_filter_reverse += 1
            else:
                pass_filter_forward += 1

            # Create CIGAR (filter already validated it)
            cigar = Cigar(hit.cigartuples or [])

            # Resolve M operations
            if RESOLVE_CIGAR_MATCHES_VS_MISMATCHES and any(op == 0 for op, _ in cigar):
                query_seq = hit.query_alignment_sequence
                ref_seq = hit.get_reference_sequence()
                cigar = resolve_cigar_m_operations(cigar, query_seq, ref_seq)

            single_alignment = SingleAlignment(
                read_type=None,  # read type will be set later
                start=hit.reference_start,
                end=hit.reference_end,
                is_reverse=hit.is_reverse,
                cigar=cigar,
                bam_index=bam_index,
                alignment_score=hit.get_tag("AS") if hit.has_tag("AS") else 0,
                read_id=hit.query_name,  # TODO: not really needed, remove?
            )

            grouped_hits[hit.reference_name][hit.query_name][hit.is_reverse].append(single_alignment)

    print(
        f"BAM statistics: total hits={fmt_count(total_hits)}, "
        f"pass_filter={fmt_count(pass_filter_total, total_hits)}, "
        f"FWD={fmt_count(pass_filter_forward, pass_filter_total)}, "
        f"REV={fmt_count(pass_filter_reverse, pass_filter_total)}",
        flush=True
    )

    return grouped_hits

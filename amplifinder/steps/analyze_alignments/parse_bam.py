"""Parsing junction coverage from BAM files."""

import pysam

from pathlib import Path
from typing import Dict

from amplifinder.data_types import JunctionReadCounts


def get_junction_coverage(
    bam_path: Path,
    read_length: int,
    read_length_tolerance: float = 0.1,
    min_overlap: int = 12,
) -> Dict[str, JunctionReadCounts]:
    """Parse BAM and get coverage for all 7 junction types.

    Args:
        bam_path: Path to sorted BAM file
        read_length: Read length for filtering
        read_length_tolerance: Tolerance for read length filtering (default 0.1 = 10%)
        min_overlap: Minimum overlap to count as spanning

    Returns:
        dict mapping reference name -> JunctionReadCounts
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    read_length_factor = 1 + read_length_tolerance
    min_len = read_length / read_length_factor
    max_len = read_length * read_length_factor

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))
        junction_points = {name: length // 2 for name, length in ref_lengths.items()}

        # [left, right, spanning] for each junction type
        counts = {name: JunctionReadCounts() for name in ref_lengths.keys()}

        for read in bam.fetch():
            if read.is_unmapped:
                continue

            ref_name = read.reference_name
            assert ref_name in junction_points, f"Reference name {ref_name} not found in junction points"
            length = read.query_alignment_length
            if not (min_len <= length <= max_len):
                continue

            start = read.reference_start
            end = read.reference_end
            junction_point = junction_points[ref_name]

            if end <= junction_point:
                counts[ref_name].left += 1      # left of the junction
            elif start > junction_point:
                counts[ref_name].right += 1     # right of the junction
            elif start <= junction_point - min_overlap and end >= junction_point + min_overlap:
                counts[ref_name].spanning += 1  # spanning the junction
            else:
                counts[ref_name].other += 1     # other (partial overlap)

    return counts

"""BAM file parsing utilities using pysam."""

from pathlib import Path
from typing import List, NamedTuple

import pysam


class AlignedRead(NamedTuple):
    """Parsed alignment information."""
    ref_name: str      # reference name (junction ID: 1-7)
    start: int         # 0-based start position
    end: int           # 0-based end position (exclusive)
    length: int        # alignment length
    is_mapped: bool    # True if read is mapped
    mapq: int          # mapping quality


class JunctionReadCounts(NamedTuple):
    """Read counts at a junction."""
    left: int      # reads on left side of junction
    right: int     # reads on right side of junction
    spanning: int  # reads spanning the junction


def parse_bam_reads(bam_path: Path) -> List[AlignedRead]:
    """Parse BAM file and extract alignment information.
    
    Args:
        bam_path: Path to sorted BAM file
    
    Returns:
        List of AlignedRead records
    
    Raises:
        FileNotFoundError: If BAM file doesn't exist
    """
    
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")
    
    reads = []
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch():
            if read.is_unmapped:
                continue
            
            reads.append(AlignedRead(
                ref_name=read.reference_name,
                start=read.reference_start,
                end=read.reference_end,
                length=read.query_alignment_length,
                is_mapped=not read.is_unmapped,
                mapq=read.mapping_quality,
            ))
    
    return reads


def count_junction_reads(
    reads: List[AlignedRead],
    junction_length: int,
    read_length: int,
    req_overlap: int = 12,
) -> dict[str, JunctionReadCounts]:
    """Count reads at each junction position.
    
    Junction structure: [left_side | junction_point | right_side]
    Junction point is at position junction_length // 2
    
    A read is classified as:
    - spanning: crosses the junction point with sufficient overlap on both sides
    - left: ends before or at the junction point
    - right: starts at or after the junction point
    
    Args:
        reads: List of aligned reads
        junction_length: Length of synthetic junction sequence (typically 2 * read_length)
        read_length: Read length for calculating overlap requirements
        req_overlap: Minimum overlap on each side to count as spanning
    
    Returns:
        Dict mapping junction ID (ref_name) to JunctionReadCounts
    """
    # Junction point is in the middle
    junction_point = junction_length // 2
    
    # Group reads by reference (junction ID)
    reads_by_ref: dict[str, List[AlignedRead]] = {}
    for read in reads:
        if read.ref_name not in reads_by_ref:
            reads_by_ref[read.ref_name] = []
        reads_by_ref[read.ref_name].append(read)
    
    results = {}
    for ref_name, ref_reads in reads_by_ref.items():
        left_count = 0
        right_count = 0
        spanning_count = 0
        
        for read in ref_reads:
            # Check if read spans junction with sufficient overlap
            left_overlap = junction_point - read.start
            right_overlap = read.end - junction_point
            
            if left_overlap >= req_overlap and right_overlap >= req_overlap:
                spanning_count += 1
            elif read.end <= junction_point:
                left_count += 1
            elif read.start >= junction_point:
                right_count += 1
            else:
                # Partial overlap - classify based on center
                read_center = (read.start + read.end) // 2
                if read_center < junction_point:
                    left_count += 1
                else:
                    right_count += 1
        
        results[ref_name] = JunctionReadCounts(
            left=left_count,
            right=right_count,
            spanning=spanning_count,
        )
    
    return results


def get_junction_coverage(
    bam_path: Path,
    junction_length: int,
    read_length: int,
    req_overlap: int = 12,
) -> List[JunctionReadCounts]:
    """Parse BAM and get coverage for all 7 junction types.
    
    Args:
        bam_path: Path to sorted BAM file
        junction_length: Length of synthetic junction sequences
        read_length: Read length
        req_overlap: Minimum overlap to count as spanning
    
    Returns:
        List of 7 JunctionReadCounts (indexed 0-6 for junction types 1-7)
    """
    reads = parse_bam_reads(bam_path)
    counts = count_junction_reads(reads, junction_length, read_length, req_overlap)
    
    # Build result list for junction types 1-7
    result = []
    for jc_type in range(1, 8):
        ref_name = str(jc_type)
        if ref_name in counts:
            result.append(counts[ref_name])
        else:
            result.append(JunctionReadCounts(left=0, right=0, spanning=0))
    
    return result

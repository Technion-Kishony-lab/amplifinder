#!/usr/bin/env python3
"""Test script for combine_same_id_same_orientation_hits function."""

from pathlib import Path
from amplifinder.steps.jct_coverage.combine_hits import combine_same_id_same_orientation_hits
from amplifinder.steps.jct_coverage.alignment_data import CombinedSingleAlignment

# Test with a real BAM file
bam_path = Path("path/to/your/test.bam")  # Update this path

if bam_path.exists():
    print(f"Processing: {bam_path}")
    combined = combine_same_id_same_orientation_hits(bam_path)
    
    print(f"\nTotal combined alignments: {len(combined)}")
    
    # Find combined alignments (those that merged multiple hits)
    combined_alignments = {k: v for k, v in combined.items() if isinstance(v, CombinedSingleAlignment) and len(v.bam_indices) > 1}
    print(f"Alignments combined from multiple hits: {len(combined_alignments)}")
    
    # Show examples of combined alignments
    print("\n=== Examples of Combined Alignments ===")
    for i, (key, alignment) in enumerate(list(combined_alignments.items())[:5]):
        read_id, is_reverse, ref_name = key
        strand = "-" if is_reverse else "+"
        print(f"\n{i+1}. {read_id} [{strand}] on {ref_name}")
        print(f"   Combined span: {alignment.start} - {alignment.end} ({alignment.end - alignment.start} bp)")
        print(f"   Number of original hits: {len(alignment.bam_indices)}")
        print(f"   Original BAM indices: {alignment.bam_indices}")
        print(f"   Original segments:")
        for j, (start, end, cigar) in enumerate(zip(alignment.starts, alignment.ends, alignment.cigars)):
            print(f"      {j+1}. {start}-{end} ({end-start} bp) CIGAR: {cigar}")
    
    # Check for alignments with gaps between segments
    with_gaps = {}
    for k, v in combined_alignments.items():
        for j in range(1, len(v.starts)):
            if v.starts[j] > v.ends[j-1]:  # Gap exists
                with_gaps[k] = v
                break
    
    print(f"\n\n=== Alignments with Gaps (non-overlapping hits) ===")
    print(f"Count: {len(with_gaps)}")
    
    for i, (key, alignment) in enumerate(list(with_gaps.items())[:3]):
        read_id, is_reverse, ref_name = key
        strand = "-" if is_reverse else "+"
        print(f"\n{i+1}. {read_id} [{strand}] on {ref_name}")
        print(f"   Combined span: {alignment.start} - {alignment.end} ({alignment.end - alignment.start} bp)")
        print(f"   Original segments:")
        for j, (start, end) in enumerate(zip(alignment.starts, alignment.ends)):
            gap_info = ""
            if j > 0:
                gap = start - alignment.ends[j-1]
                if gap > 0:
                    gap_info = f" (gap: {gap} bp)"
            print(f"      {j+1}. {start}-{end} ({end-start} bp){gap_info}")
else:
    print(f"BAM file not found: {bam_path}")
    print("\nUpdate the bam_path variable in this script to test.")

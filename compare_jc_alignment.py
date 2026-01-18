#!/usr/bin/env python3
"""
Compare MATLAB and Python Junction Alignment Results

This script performs a comprehensive comparison between MATLAB and Python pipeline
outputs for junction alignment analysis, identifying differences in BAM file
alignments and read-type classification (left/right/green).

Usage:
    python compare_jc_alignment.py \\
        --matlab-bam <path> \\
        --python-bam <path> \\
        --matlab-fasta <path> \\
        --python-fasta <path> \\
        --matlab-data-dir <path> \\
        --output-dir <path>
"""

import argparse
import csv
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any

import pandas as pd
import pysam
from Bio import SeqIO

# Import pipeline utilities for samtools path detection
from amplifinder.env import SAMTOOLS_PATH
from amplifinder.utils.run_utils import get_tool_path

# Junction type mapping (from plan)
MATLAB_TO_PYTHON_JUNCTION_MAP = {
    '1': 'CHR_TO_AMP_LEFT',
    '2': 'CHR_TO_TN_LEFT',
    '3': 'AMP_RIGHT_TO_TN_LEFT',
    '4': 'AMP_RIGHT_TO_AMP_LEFT',
    '5': 'TN_RIGHT_TO_AMP_LEFT',
    '6': 'TN_RIGHT_TO_CHR',
    '7': 'AMP_RIGHT_TO_CHR'
}

PYTHON_TO_MATLAB_JUNCTION_MAP = {v: k for k, v in MATLAB_TO_PYTHON_JUNCTION_MAP.items()}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description='Compare MATLAB and Python junction alignment results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument(
        '--matlab-bam',
        type=Path,
        required=True,
        help='Path to MATLAB BAM file'
    )
    parser.add_argument(
        '--python-bam',
        type=Path,
        required=True,
        help='Path to Python BAM file'
    )
    parser.add_argument(
        '--matlab-fasta',
        type=Path,
        required=True,
        help='Path to MATLAB FASTA file'
    )
    parser.add_argument(
        '--python-fasta',
        type=Path,
        required=True,
        help='Path to Python FASTA file'
    )
    parser.add_argument(
        '--matlab-data-dir',
        type=Path,
        required=True,
        help='Directory containing MATLAB CSV files (exported from .mat files)'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path('./comparison_outputs'),
        help='Output directory for comparison results (default: ./comparison_outputs)'
    )
    parser.add_argument(
        '--junction-length',
        type=int,
        default=None,
        help='Junction sequence length (auto-detected from FASTA if not provided)'
    )
    parser.add_argument(
        '--read-length',
        type=int,
        default=150,
        help='Average read length (default: 150)'
    )
    parser.add_argument(
        '--min-overlap-len',
        type=int,
        default=12,
        help='Minimum overlap length (default: 12)'
    )
    parser.add_argument(
        '--read-length-tolerance',
        type=float,
        default=0.1,
        help='Read length tolerance (default: 0.1)'
    )
    parser.add_argument(
        '--min-bp-in-frame',
        type=int,
        default=10,
        help='Minimum bp in frame (default: 10)'
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    parser.add_argument(
        '--skip-bam-comparison',
        action='store_true',
        help='Skip Part 1 (BAM comparison)'
    )
    parser.add_argument(
        '--skip-read-counts',
        action='store_true',
        help='Skip Part 2 (read count comparison)'
    )
    return parser.parse_args()


def log(message: str, verbose: bool = True) -> None:
    """Print log message if verbose."""
    if verbose:
        print(message, flush=True)


def ensure_dir(path: Path) -> None:
    """Ensure directory exists."""
    path.mkdir(parents=True, exist_ok=True)


def get_samtools_path() -> Path:
    """Get samtools path using pipeline utilities."""
    try:
        return get_tool_path("samtools", config_path=SAMTOOLS_PATH, required=True)
    except FileNotFoundError:
        # Fallback: try to find in PATH
        try:
            result = subprocess.run(['which', 'samtools'], capture_output=True, text=True, check=True)
            return Path(result.stdout.strip())
        except (subprocess.CalledProcessError, FileNotFoundError):
            raise FileNotFoundError("samtools not found. Please install samtools or configure samtools_path in config.yaml")


def get_junction_length(fasta_path: Path) -> int:
    """Get junction length from FASTA file."""
    for record in SeqIO.parse(fasta_path, 'fasta'):
        return len(record.seq)
    raise ValueError(f"No sequences found in {fasta_path}")


def is_cigar_only_m(read: pysam.AlignedSegment) -> bool:
    """Check if read has CIGAR-only-M (all match operations)."""
    if not read.cigartuples:
        return False
    return all(op == 0 for op, _ in read.cigartuples)


def classify_read_matlab(
    start_1: int,
    alignment_length: int,
    jct_len: int,
    read_length: int,
    min_overlap_len: int,
    min_bp_in_frame: int,
) -> Optional[str]:
    """Classify read using MATLAB logic.
    
    Args:
        start_1: 1-based start position (inclusive)
        alignment_length: Length of alignment
        jct_len: Junction length
        read_length: Expected read length
        min_overlap_len: Minimum overlap length
        min_bp_in_frame: Minimum bp in frame
    """
    jct_point = jct_len / 2
    read_end = start_1 + alignment_length - 1  # 1-based inclusive end
    
    # Check read length
    okz_length = (alignment_length > read_length * 0.9) and (alignment_length < read_length * 1.1)
    if not okz_length:
        return None
    
    # Green classification: start before junction - overlap, end after junction + overlap
    okz_green = (start_1 < jct_point - min_overlap_len) and (read_end > jct_point + min_overlap_len)
    
    # Left classification: extends into left region but ends before junction + overlap
    okz_left = (start_1 + alignment_length > jct_point - read_length + min_bp_in_frame) and \
               (read_end < jct_point + min_overlap_len)
    
    # Right classification: starts after junction - overlap, but before junction + read_length - min_bp
    okz_right = (start_1 > jct_point - min_overlap_len) and \
                (start_1 < jct_point + read_length - min_bp_in_frame)
    
    if okz_green:
        return "green"
    elif okz_left:
        return "left"
    elif okz_right:
        return "right"
    return None


def classify_read_python(
    read: pysam.AlignedSegment,
    jct_len: int,
    avg_read_length: int,
    min_overlap_len: int,
    read_length_tolerance: float = 0.1,
    max_dist_from_junction: int = 10,
) -> Optional[str]:
    """Classify read using Python logic (from JunctionReadCounts.get_read_type)."""
    alignment_length = read.query_alignment_length
    read_length_factor = 1 + read_length_tolerance
    min_alignment_length = avg_read_length / read_length_factor
    max_alignment_length = avg_read_length * read_length_factor
    
    if not (min_alignment_length <= alignment_length <= max_alignment_length):
        return None
    
    if not is_cigar_only_m(read):
        return None
    
    # Filter by alignment tags
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
    
    if (nm_score is not None and nm_score > 3) or (as_score is not None and as_score < -25):
        return None
    
    # Convert to 1-based coordinates
    start_1 = read.reference_start + 1
    end_1 = read.reference_end
    arm_len = jct_len // 2
    
    # Use Python classification logic (from JunctionReadCounts.get_read_type)
    idx_L_1 = arm_len - end_1
    idx_L_2 = arm_len - start_1
    idx_R_1 = start_1 - (arm_len + 1)
    idx_R_2 = end_1 - (arm_len + 1)
    
    if idx_L_1 > max_dist_from_junction:
        return None  # too far from junction on left
    if idx_R_1 > max_dist_from_junction:
        return None  # too far from junction on right
    
    if idx_L_1 >= 0:
        return "left"
    elif idx_R_1 >= 0:
        return "right"
    
    is_start_left_of_margin = idx_L_2 >= min_overlap_len - 1
    is_end_right_of_margin = idx_R_2 >= min_overlap_len - 1
    
    if is_start_left_of_margin and is_end_right_of_margin:
        return "green"  # MIDDLE/spanning
    elif is_start_left_of_margin:
        return "left"  # LEFT_MARGINAL -> left
    elif is_end_right_of_margin:
        return "right"  # RIGHT_MARGINAL -> right
    
    return None


def extract_all_tags(read: pysam.AlignedSegment) -> Dict[str, Any]:
    """Extract all tags from a read."""
    tags = {}
    for tag, value in read.get_tags():
        tags[tag] = value
    return tags


def extract_alignment_record(read: pysam.AlignedSegment, source: str) -> Dict[str, Any]:
    """Extract all alignment fields and tags from a read."""
    # Get all tags first
    all_tags = extract_all_tags(read)

    ref_name = read.reference_name if read.reference_name else '*'
    pos_1 = read.reference_start + 1 if read.reference_start is not None and read.reference_start >= 0 else 0
    rnext = read.next_reference_name if read.next_reference_name else '*'
    pnext = read.next_reference_start + 1 if read.next_reference_start is not None else 0
    seq = read.query_sequence if read.query_sequence else ''
    qual = ''.join(chr(q + 33) for q in read.query_qualities) if read.query_qualities else ''

    # Standard SAM fields
    record = {
        'read_id': read.query_name,
        'source': source,
        'flag': read.flag,
        'ref_name': ref_name,
        'pos': pos_1,  # 1-based
        'mapq': read.mapping_quality,
        'cigar': read.cigarstring,
        'rnext': rnext,
        'pnext': pnext,
        'tlen': read.template_length,
        'sequence': seq,
        'quality': qual,
    }
    
    # Flag-derived boolean fields
    record.update({
        'is_paired': read.is_paired,
        'is_proper_pair': read.is_proper_pair,
        'is_unmapped': read.is_unmapped,
        'mate_unmapped': read.mate_is_unmapped,
        'is_reverse': read.is_reverse,
        'mate_reverse': read.mate_is_reverse,
        'is_read1': read.is_read1,
        'is_read2': read.is_read2,
        'is_secondary': read.is_secondary,
        'is_qcfail': read.is_qcfail,
        'is_duplicate': read.is_duplicate,
        'is_supplementary': read.is_supplementary,
    })
    
    # Coordinate fields
    record.update({
        'ref_start': read.reference_start,  # 0-based
        'ref_end': read.reference_end,  # 0-based exclusive
        'query_alignment_start': read.query_alignment_start,
        'query_alignment_end': read.query_alignment_end,
        'query_length': read.query_length,
    })
    
    # Add all tags (common ones explicitly, others dynamically)
    common_tags = ['AS', 'NM', 'nM', 'NH', 'HI', 'XS', 'MD', 'MC', 'MQ', 'SA', 'RG', 'BC', 'OC', 'OP', 'OA']
    for tag in common_tags:
        record[f'tag_{tag}'] = all_tags.get(tag, 'N/A')
    
    # Add any other tags found
    for tag, value in all_tags.items():
        if tag not in common_tags:
            record[f'tag_{tag}'] = value
    
    return record


def convert_bam_to_sam(bam_path: Path, sam_path: Path, verbose: bool = True) -> None:
    """Convert BAM to SAM using samtools."""
    log(f"Converting {bam_path.name} to SAM...", verbose)
    samtools = get_samtools_path()
    cmd = [str(samtools), 'view', '-h', str(bam_path)]
    with open(sam_path, 'w') as f:
        subprocess.run(cmd, stdout=f, check=True)
    log(f"  Saved to {sam_path}", verbose)


def split_unmapped_reads(
    bam_path: Path,
    mapped_bam_path: Path,
    unmapped_bam_path: Path,
    verbose: bool = True,
) -> Tuple[int, int]:
    """Split BAM into mapped and unmapped reads."""
    log(f"Splitting unmapped reads in {bam_path.name}...", verbose)
    mapped_count = 0
    unmapped_count = 0
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam_in:
        with pysam.AlignmentFile(str(mapped_bam_path), 'wb', header=bam_in.header) as bam_mapped:
            with pysam.AlignmentFile(str(unmapped_bam_path), 'wb', header=bam_in.header) as bam_unmapped:
                for read in bam_in:
                    if read.is_unmapped:
                        bam_unmapped.write(read)
                        unmapped_count += 1
                    else:
                        bam_mapped.write(read)
                        mapped_count += 1
    log(f"  Mapped: {mapped_count}, Unmapped: {unmapped_count}", verbose)
    log(f"  Mapped BAM: {mapped_bam_path}", verbose)
    log(f"  Unmapped BAM: {unmapped_bam_path}", verbose)
    return mapped_count, unmapped_count


def write_fastq_from_records(records: List[Dict[str, Any]], fastq_path: Path, verbose: bool = True) -> int:
    """Write FASTQ from alignment records (uses first record per read_id)."""
    seen = set()
    written = 0
    with open(fastq_path, 'w') as f:
        for record in records:
            read_id = record.get('read_id')
            if not read_id or read_id in seen:
                continue
            seq = record.get('sequence') or ''
            if not seq:
                continue
            qual = record.get('quality') or ''
            if not qual:
                qual = 'I' * len(seq)
            f.write(f"@{read_id}\n{seq}\n+\n{qual}\n")
            seen.add(read_id)
            written += 1
    if verbose:
        log(f"  Wrote {written} reads to {fastq_path}", verbose)
    return written


def compare_bam_files(
    matlab_bam: Path,
    python_bam: Path,
    output_dir: Path,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Part 1: Compare BAM files read-by-read."""
    log("\n=== Part 1: BAM File Comparison ===", verbose)
    
    bam_output_dir = output_dir / 'bam_comparison'
    ensure_dir(bam_output_dir)
    
    # Split unmapped reads to separate files (compare only mapped)
    matlab_mapped_bam = bam_output_dir / 'matlab_alignment.mapped.bam'
    matlab_unmapped_bam = bam_output_dir / 'matlab_alignment.unmapped.bam'
    python_mapped_bam = bam_output_dir / 'python_alignment.mapped.bam'
    python_unmapped_bam = bam_output_dir / 'python_alignment.unmapped.bam'

    split_unmapped_reads(matlab_bam, matlab_mapped_bam, matlab_unmapped_bam, verbose)
    split_unmapped_reads(python_bam, python_mapped_bam, python_unmapped_bam, verbose)

    # Convert BAM to SAM
    matlab_sam = bam_output_dir / 'matlab_alignment.sam'
    python_sam = bam_output_dir / 'python_alignment.sam'
    
    convert_bam_to_sam(matlab_mapped_bam, matlab_sam, verbose)
    convert_bam_to_sam(python_mapped_bam, python_sam, verbose)
    
    # Parse BAM files and extract alignments
    log("Parsing MATLAB BAM file...", verbose)
    matlab_alignments = defaultdict(list)  # read_id -> list of alignments
    matlab_all_tags = set()
    
    with pysam.AlignmentFile(str(matlab_mapped_bam), 'rb') as bam:
        for read in bam:
            read_id = read.query_name
            record = extract_alignment_record(read, 'matlab')
            matlab_alignments[read_id].append(record)
            matlab_all_tags.update(extract_all_tags(read).keys())
    
    log(f"  Found {len(matlab_alignments)} unique reads with {sum(len(v) for v in matlab_alignments.values())} alignments", verbose)
    
    log("Parsing Python BAM file...", verbose)
    python_alignments = defaultdict(list)
    python_all_tags = set()
    
    with pysam.AlignmentFile(str(python_mapped_bam), 'rb') as bam:
        for read in bam:
            read_id = read.query_name
            record = extract_alignment_record(read, 'python')
            python_alignments[read_id].append(record)
            python_all_tags.update(extract_all_tags(read).keys())
    
    log(f"  Found {len(python_alignments)} unique reads with {sum(len(v) for v in python_alignments.values())} alignments", verbose)
    
    # Collect all unique tag names
    all_tags = sorted(matlab_all_tags | python_all_tags)
    common_tags = ['AS', 'NM', 'nM', 'NH', 'HI', 'XS', 'MD', 'MC', 'MQ', 'SA', 'RG', 'BC', 'OC', 'OP', 'OA']
    other_tags = sorted(set(all_tags) - set(common_tags))
    
    # Compare alignments
    log("Comparing alignments...", verbose)
    matched = []
    matlab_only = []
    python_only = []
    differing = []
    
    all_read_ids = set(matlab_alignments.keys()) | set(python_alignments.keys())
    
    for read_id in sorted(all_read_ids):
        matlab_reads = matlab_alignments.get(read_id, [])
        python_reads = python_alignments.get(read_id, [])
        
        if not matlab_reads:
            python_only.extend(python_reads)
            continue
        if not python_reads:
            matlab_only.extend(matlab_reads)
            continue
        
        # Try to match alignments
        matched_reads = []
        unmatched_matlab = list(matlab_reads)
        unmatched_python = list(python_reads)
        
        for m_read in matlab_reads:
            for p_read in python_reads:
                # Convert MATLAB ref_name (numeric) to Python ref_name (descriptive) for comparison
                matlab_ref_mapped = MATLAB_TO_PYTHON_JUNCTION_MAP.get(m_read['ref_name'], m_read['ref_name'])
                python_ref = p_read['ref_name']
                
                # Match based on ref_name (after mapping), position (within 1 bp), and CIGAR
                if (matlab_ref_mapped == python_ref and
                    abs(m_read['pos'] - p_read['pos']) <= 1 and
                    m_read['cigar'] == p_read['cigar'] and
                    m_read['flag'] == p_read['flag']):
                    matched_reads.append((m_read, p_read))
                    if m_read in unmatched_matlab:
                        unmatched_matlab.remove(m_read)
                    if p_read in unmatched_python:
                        unmatched_python.remove(p_read)
                    break
        
        if matched_reads:
            matched.extend([m for m, p in matched_reads])
            matched.extend([p for m, p in matched_reads])
        
        if unmatched_matlab:
            differing.extend(unmatched_matlab)
        if unmatched_python:
            differing.extend(unmatched_python)
    
    log(f"  Matched: {len(matched) // 2} alignments", verbose)
    log(f"  MATLAB only: {len(matlab_only)} alignments", verbose)
    log(f"  Python only: {len(python_only)} alignments", verbose)
    log(f"  Differing: {len(differing)} alignments", verbose)
    
    # Write differences CSV
    if differing or matlab_only or python_only:
        log("Writing differences CSV...", verbose)
        all_differences = differing + matlab_only + python_only
        
        # Build column list
        base_columns = [
            'read_id', 'source', 'flag', 'ref_name', 'pos', 'mapq', 'cigar',
            'rnext', 'pnext', 'tlen', 'sequence', 'quality',
            'is_paired', 'is_proper_pair', 'is_unmapped', 'mate_unmapped',
            'is_reverse', 'mate_reverse', 'is_read1', 'is_read2',
            'is_secondary', 'is_qcfail', 'is_duplicate', 'is_supplementary',
            'ref_start', 'ref_end', 'query_alignment_start', 'query_alignment_end', 'query_length'
        ]
        tag_columns = [f'tag_{tag}' for tag in common_tags + other_tags]
        columns = base_columns + tag_columns
        
        # Ensure all records have all columns
        for record in all_differences:
            for col in columns:
                if col not in record:
                    record[col] = 'N/A'
        
        df = pd.DataFrame(all_differences)
        # Sort by read_id, then source, then ref_name, then pos
        df = df.sort_values(['read_id', 'source', 'ref_name', 'pos'])
        diff_csv_path = bam_output_dir / 'bam_comparison_differences.csv'
        df.to_csv(diff_csv_path, index=False)
        log(f"  Saved to {diff_csv_path}", verbose)

        # Write FASTQ for differing reads
        fastq_path = bam_output_dir / 'green_mismatch_reads.fastq'
        write_fastq_from_records(all_differences, fastq_path, verbose)
    
    # Write summary
    summary_path = bam_output_dir / 'comparison_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("BAM File Comparison Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"MATLAB BAM: {matlab_bam}\n")
        f.write(f"Python BAM: {python_bam}\n\n")
        f.write(f"Total unique reads in MATLAB: {len(matlab_alignments)}\n")
        f.write(f"Total unique reads in Python: {len(python_alignments)}\n")
        f.write(f"Total alignments in MATLAB: {sum(len(v) for v in matlab_alignments.values())}\n")
        f.write(f"Total alignments in Python: {sum(len(v) for v in python_alignments.values())}\n\n")
        f.write(f"Matched alignments: {len(matched) // 2}\n")
        f.write(f"MATLAB-only alignments: {len(matlab_only)}\n")
        f.write(f"Python-only alignments: {len(python_only)}\n")
        f.write(f"Differing alignments: {len(differing)}\n")
    
    log(f"Summary saved to {summary_path}", verbose)
    
    return {
        'matlab_reads': len(matlab_alignments),
        'python_reads': len(python_alignments),
        'matched': len(matched) // 2,
        'matlab_only': len(matlab_only),
        'python_only': len(python_only),
        'differing': len(differing),
    }


def load_matlab_csv(matlab_data_dir: Path, filename: str) -> pd.Series:
    """Load MATLAB CSV file."""
    path = matlab_data_dir / filename
    if not path.exists():
        raise FileNotFoundError(f"MATLAB CSV file not found: {path}")
    df = pd.read_csv(path, header=None)
    values = df.values.flatten()
    return pd.Series(values)


def compare_read_counts(
    matlab_bam: Path,
    python_bam: Path,
    matlab_data_dir: Path,
    output_dir: Path,
    junction_length: int,
    read_length: int,
    min_overlap_len: int,
    read_length_tolerance: float,
    min_bp_in_frame: int,
    verbose: bool = True,
) -> Dict[str, Any]:
    """Part 2: Compare read counts and classifications."""
    log("\n=== Part 2: Read Count Comparison ===", verbose)
    
    count_output_dir = output_dir / 'read_count_comparison'
    ensure_dir(count_output_dir)
    
    # Load MATLAB data
    log("Loading MATLAB data...", verbose)
    try:
        RdStart = load_matlab_csv(matlab_data_dir, 'biomap__RdStart.csv')
        RdRef = load_matlab_csv(matlab_data_dir, 'biomap__RdRef.csv')
        RdLength = load_matlab_csv(matlab_data_dir, 'biomap__RdLength.csv')
        RdFlag = load_matlab_csv(matlab_data_dir, 'biomap__RdFlag.csv')
        RdFull = load_matlab_csv(matlab_data_dir, 'biomap__RdFull.csv')
        nmbr_green = load_matlab_csv(matlab_data_dir, 'bamreads__nmbr_green_reads.csv')
        nmbr_left = load_matlab_csv(matlab_data_dir, 'bamreads__nmbr_left_reads.csv')
        nmbr_right = load_matlab_csv(matlab_data_dir, 'bamreads__nmbr_right_reads.csv')
    except FileNotFoundError as e:
        raise FileNotFoundError(f"Missing MATLAB CSV file: {e}. Please export .mat files to CSV first.")
    
    log(f"  Loaded {len(RdStart)} MATLAB alignment records", verbose)
    
    # Build MATLAB index-to-readID mapping
    log("Building MATLAB index-to-readID mapping...", verbose)
    matlab_index_to_readid = {}
    matlab_alignments_by_index = {}
    
    with pysam.AlignmentFile(str(matlab_bam), 'rb') as bam:
        for idx, read in enumerate(bam):
            read_id = read.query_name
            matlab_index_to_readid[idx] = read_id
            matlab_alignments_by_index[idx] = {
                'read_id': read_id,
                'ref_name': read.reference_name,
                'start': read.reference_start + 1,
                'end': read.reference_end,
                'cigar': read.cigarstring,
                'flag': read.flag,
                'mapq': read.mapping_quality,
                'is_paired': read.is_paired,
                'is_proper_pair': read.is_proper_pair,
                'is_unmapped': read.is_unmapped,
                'mate_unmapped': read.mate_is_unmapped,
                'is_reverse': read.is_reverse,
                'mate_reverse': read.mate_is_reverse,
                'is_read1': read.is_read1,
                'is_read2': read.is_read2,
                'is_secondary': read.is_secondary,
                'is_qcfail': read.is_qcfail,
                'is_duplicate': read.is_duplicate,
                'is_supplementary': read.is_supplementary,
                'tags': extract_all_tags(read),
            }
    
    log(f"  Mapped {len(matlab_index_to_readid)} indices to read IDs", verbose)
    
    # Extract MATLAB classifications
    log("Extracting MATLAB classifications...", verbose)
    matlab_classifications = {}  # (read_id, jct_type) -> classification
    matlab_classification_details = {}
    
    for i in range(len(RdStart)):
        if i not in matlab_index_to_readid:
            continue  # Skip if index not in BAM
        
        if not RdFull.iloc[i]:  # Skip if not CIGAR-only-M
            continue
        
        alignment_length = int(RdLength.iloc[i])
        okz_length = (alignment_length > read_length * 0.9) and (alignment_length < read_length * 1.1)
        if not okz_length:
            continue
        
        jct_type = str(RdRef.iloc[i])  # '1'..'7'
        start_1 = int(RdStart.iloc[i])
        read_end = start_1 + alignment_length - 1
        
        classification = classify_read_matlab(
            start_1, alignment_length,
            junction_length, read_length, min_overlap_len, min_bp_in_frame
        )
        
        read_id = matlab_index_to_readid[i]
        key = (read_id, jct_type)
        matlab_classifications.setdefault(key, []).append(classification)
        matlab_details = matlab_alignments_by_index.get(i, {})
        matlab_classification_details[key] = {
            'classification': classification,
            'start': start_1,
            'end': read_end,
            'length': alignment_length,
            'flag': int(RdFlag.iloc[i]),
            'matlab_index': i,
            'cigar': matlab_details.get('cigar'),
            'mapq': matlab_details.get('mapq'),
            'is_paired': matlab_details.get('is_paired'),
            'is_proper_pair': matlab_details.get('is_proper_pair'),
            'is_unmapped': matlab_details.get('is_unmapped'),
            'mate_unmapped': matlab_details.get('mate_unmapped'),
            'is_reverse': matlab_details.get('is_reverse'),
            'mate_reverse': matlab_details.get('mate_reverse'),
            'is_read1': matlab_details.get('is_read1'),
            'is_read2': matlab_details.get('is_read2'),
            'is_secondary': matlab_details.get('is_secondary'),
            'is_qcfail': matlab_details.get('is_qcfail'),
            'is_duplicate': matlab_details.get('is_duplicate'),
            'is_supplementary': matlab_details.get('is_supplementary'),
            'tags': matlab_details.get('tags', {}),
        }
    
    log(f"  Classified {len(matlab_classifications)} MATLAB read-junction pairs", verbose)
    
    # Extract Python classifications
    log("Extracting Python classifications...", verbose)
    python_classifications = {}
    python_classification_details = {}
    jct_lengths = {}
    
    with pysam.AlignmentFile(str(python_bam), 'rb') as bam:
        jct_lengths = dict(zip(bam.references, bam.lengths))
        
        for read in bam:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            
            read_id = read.query_name
            jct_name = read.reference_name
            jct_type = PYTHON_TO_MATLAB_JUNCTION_MAP.get(jct_name)
            
            if jct_type is None:
                continue  # Unknown junction type
            
            jct_len = jct_lengths[jct_name]
            
            classification = classify_read_python(
                read, jct_len, read_length, min_overlap_len,
                read_length_tolerance, min_bp_in_frame
            )
            
            key = (read_id, jct_type)
            python_classifications.setdefault(key, []).append(classification)
            
            start_1 = read.reference_start + 1
            end_1 = read.reference_end
            
            python_classification_details[key] = {
                'classification': classification,
                'start': start_1,
                'end': end_1,
                'length': read.query_alignment_length,
                'cigar': read.cigarstring,
                'flag': read.flag,
                'mapq': read.mapping_quality,
                'is_paired': read.is_paired,
                'is_proper_pair': read.is_proper_pair,
                'is_unmapped': read.is_unmapped,
                'mate_unmapped': read.mate_is_unmapped,
                'is_reverse': read.is_reverse,
                'mate_reverse': read.mate_is_reverse,
                'is_read1': read.is_read1,
                'is_read2': read.is_read2,
                'is_secondary': read.is_secondary,
                'is_qcfail': read.is_qcfail,
                'is_duplicate': read.is_duplicate,
                'is_supplementary': read.is_supplementary,
                'tags': extract_all_tags(read),
            }
    
    log(f"  Classified {len(python_classifications)} Python read-junction pairs", verbose)
    
    # Compare summary counts
    log("Comparing summary counts...", verbose)
    summary_data = []
    
    for jct_type in sorted(MATLAB_TO_PYTHON_JUNCTION_MAP.keys()):
        # MATLAB counts (from bamreads.mat)
        matlab_green = int(nmbr_green.iloc[int(jct_type) - 1]) if int(jct_type) <= len(nmbr_green) else 0
        matlab_left = int(nmbr_left.iloc[int(jct_type) - 1]) if int(jct_type) <= len(nmbr_left) else 0
        matlab_right = int(nmbr_right.iloc[int(jct_type) - 1]) if int(jct_type) <= len(nmbr_right) else 0
        
        # Python counts (count from classifications)
        python_green = sum(1 for (rid, jt), cls_list in python_classifications.items()
                           if jt == jct_type for cls in cls_list if cls == 'green')
        python_left = sum(1 for (rid, jt), cls_list in python_classifications.items()
                          if jt == jct_type for cls in cls_list if cls == 'left')
        python_right = sum(1 for (rid, jt), cls_list in python_classifications.items()
                           if jt == jct_type for cls in cls_list if cls == 'right')
        
        summary_data.append({
            'junction_type': jct_type,
            'matlab_green': matlab_green,
            'python_green': python_green,
            'green_diff': python_green - matlab_green,
            'matlab_left': matlab_left,
            'python_left': python_left,
            'left_diff': python_left - matlab_left,
            'matlab_right': matlab_right,
            'python_right': python_right,
            'right_diff': python_right - matlab_right,
        })
    
    df_summary = pd.DataFrame(summary_data)
    df_summary.to_csv(count_output_dir / 'read_count_comparison_summary.csv', index=False)
    log(f"  Saved to {count_output_dir / 'read_count_comparison_summary.csv'}", verbose)
    
    # Find differentially labelled reads
    log("Finding differentially labelled reads...", verbose)
    common_keys = set(matlab_classifications.keys()) & set(python_classifications.keys())
    matlab_only = set(matlab_classifications.keys()) - set(python_classifications.keys())
    python_only = set(python_classifications.keys()) - set(matlab_classifications.keys())

    matlab_all_tags = set()
    python_all_tags = set()
    for details in matlab_classification_details.values():
        matlab_all_tags.update((details.get('tags') or {}).keys())
    for details in python_classification_details.values():
        python_all_tags.update((details.get('tags') or {}).keys())
    all_tags = sorted(matlab_all_tags | python_all_tags)
    
    diff_reads = []
    
    for key in sorted(common_keys):
        matlab_cls_list = matlab_classifications[key]
        python_cls_list = python_classifications[key]
        matlab_cls = matlab_cls_list[0] if matlab_cls_list else None
        python_cls = python_cls_list[0] if python_cls_list else None
        if matlab_cls != python_cls:
            read_id, jct_type = key
            matlab_details = matlab_classification_details.get(key, {})
            python_details = python_classification_details.get(key, {})
            
            reason = f"MATLAB: {matlab_cls}, Python: {python_cls}"
            
            row = {
                'read_id': read_id,
                'junction_type': jct_type,
                'matlab_classification': matlab_cls if matlab_cls else 'N/A',
                'python_classification': python_cls if python_cls else 'N/A',
                'difference_reason': reason,
            }
            
            # Add MATLAB details
            row.update({
                'matlab_start': matlab_details.get('start', 'N/A'),
                'matlab_end': matlab_details.get('end', 'N/A'),
                'matlab_length': matlab_details.get('length', 'N/A'),
                'matlab_flag': matlab_details.get('flag', 'N/A'),
                'matlab_cigar': matlab_details.get('cigar', 'N/A'),
                'matlab_mapq': matlab_details.get('mapq', 'N/A'),
                'matlab_is_paired': matlab_details.get('is_paired', 'N/A'),
                'matlab_is_proper_pair': matlab_details.get('is_proper_pair', 'N/A'),
                'matlab_is_unmapped': matlab_details.get('is_unmapped', 'N/A'),
                'matlab_mate_unmapped': matlab_details.get('mate_unmapped', 'N/A'),
                'matlab_is_reverse': matlab_details.get('is_reverse', 'N/A'),
                'matlab_mate_reverse': matlab_details.get('mate_reverse', 'N/A'),
                'matlab_is_read1': matlab_details.get('is_read1', 'N/A'),
                'matlab_is_read2': matlab_details.get('is_read2', 'N/A'),
                'matlab_is_secondary': matlab_details.get('is_secondary', 'N/A'),
                'matlab_is_qcfail': matlab_details.get('is_qcfail', 'N/A'),
                'matlab_is_duplicate': matlab_details.get('is_duplicate', 'N/A'),
                'matlab_is_supplementary': matlab_details.get('is_supplementary', 'N/A'),
            })
            
            # Add Python details
            row.update({
                'python_start': python_details.get('start', 'N/A'),
                'python_end': python_details.get('end', 'N/A'),
                'python_length': python_details.get('length', 'N/A'),
                'python_cigar': python_details.get('cigar', 'N/A'),
                'python_flag': python_details.get('flag', 'N/A'),
                'python_mapq': python_details.get('mapq', 'N/A'),
                'python_is_paired': python_details.get('is_paired', 'N/A'),
                'python_is_proper_pair': python_details.get('is_proper_pair', 'N/A'),
                'python_is_unmapped': python_details.get('is_unmapped', 'N/A'),
                'python_mate_unmapped': python_details.get('mate_unmapped', 'N/A'),
                'python_is_reverse': python_details.get('is_reverse', 'N/A'),
                'python_mate_reverse': python_details.get('mate_reverse', 'N/A'),
                'python_is_read1': python_details.get('is_read1', 'N/A'),
                'python_is_read2': python_details.get('is_read2', 'N/A'),
                'python_is_secondary': python_details.get('is_secondary', 'N/A'),
                'python_is_qcfail': python_details.get('is_qcfail', 'N/A'),
                'python_is_duplicate': python_details.get('is_duplicate', 'N/A'),
                'python_is_supplementary': python_details.get('is_supplementary', 'N/A'),
            })

            for tag in all_tags:
                row[f'matlab_tag_{tag}'] = (matlab_details.get('tags') or {}).get(tag, 'N/A')
                row[f'python_tag_{tag}'] = (python_details.get('tags') or {}).get(tag, 'N/A')
            
            diff_reads.append(row)
    
    # Add MATLAB-only and Python-only reads
    for key in sorted(matlab_only):
        read_id, jct_type = key
        matlab_details = matlab_classification_details.get(key, {})
        row = {
            'read_id': read_id,
            'junction_type': jct_type,
            'matlab_classification': matlab_classifications[key][0] if matlab_classifications[key] else 'N/A',
            'python_classification': 'N/A',
            'difference_reason': 'Only in MATLAB',
        }
        row.update({
            'matlab_start': matlab_details.get('start', 'N/A'),
            'matlab_end': matlab_details.get('end', 'N/A'),
            'matlab_length': matlab_details.get('length', 'N/A'),
            'matlab_flag': matlab_details.get('flag', 'N/A'),
            'matlab_cigar': matlab_details.get('cigar', 'N/A'),
            'matlab_mapq': matlab_details.get('mapq', 'N/A'),
            'matlab_is_paired': matlab_details.get('is_paired', 'N/A'),
            'matlab_is_proper_pair': matlab_details.get('is_proper_pair', 'N/A'),
            'matlab_is_unmapped': matlab_details.get('is_unmapped', 'N/A'),
            'matlab_mate_unmapped': matlab_details.get('mate_unmapped', 'N/A'),
            'matlab_is_reverse': matlab_details.get('is_reverse', 'N/A'),
            'matlab_mate_reverse': matlab_details.get('mate_reverse', 'N/A'),
            'matlab_is_read1': matlab_details.get('is_read1', 'N/A'),
            'matlab_is_read2': matlab_details.get('is_read2', 'N/A'),
            'matlab_is_secondary': matlab_details.get('is_secondary', 'N/A'),
            'matlab_is_qcfail': matlab_details.get('is_qcfail', 'N/A'),
            'matlab_is_duplicate': matlab_details.get('is_duplicate', 'N/A'),
            'matlab_is_supplementary': matlab_details.get('is_supplementary', 'N/A'),
        })
        row.update({f'python_{k}': 'N/A' for k in ['start', 'end', 'length', 'cigar', 'flag', 'mapq',
                                                   'is_paired', 'is_proper_pair', 'is_unmapped',
                                                   'mate_unmapped', 'is_reverse', 'mate_reverse',
                                                   'is_read1', 'is_read2', 'is_secondary',
                                                   'is_qcfail', 'is_duplicate', 'is_supplementary']})
        for tag in all_tags:
            row[f'matlab_tag_{tag}'] = (matlab_details.get('tags') or {}).get(tag, 'N/A')
            row[f'python_tag_{tag}'] = 'N/A'
        diff_reads.append(row)
    
    for key in sorted(python_only):
        read_id, jct_type = key
        python_details = python_classification_details.get(key, {})
        row = {
            'read_id': read_id,
            'junction_type': jct_type,
            'matlab_classification': 'N/A',
            'python_classification': python_classifications[key][0] if python_classifications[key] else 'N/A',
            'difference_reason': 'Only in Python',
        }
        row.update({f'matlab_{k}': 'N/A' for k in ['start', 'end', 'length', 'flag', 'cigar', 'mapq',
                                                   'is_paired', 'is_proper_pair', 'is_unmapped',
                                                   'mate_unmapped', 'is_reverse', 'mate_reverse',
                                                   'is_read1', 'is_read2', 'is_secondary',
                                                   'is_qcfail', 'is_duplicate', 'is_supplementary']})
        row.update({
            'python_start': python_details.get('start', 'N/A'),
            'python_end': python_details.get('end', 'N/A'),
            'python_length': python_details.get('length', 'N/A'),
            'python_cigar': python_details.get('cigar', 'N/A'),
            'python_flag': python_details.get('flag', 'N/A'),
            'python_mapq': python_details.get('mapq', 'N/A'),
            'python_is_paired': python_details.get('is_paired', 'N/A'),
            'python_is_proper_pair': python_details.get('is_proper_pair', 'N/A'),
            'python_is_unmapped': python_details.get('is_unmapped', 'N/A'),
            'python_mate_unmapped': python_details.get('mate_unmapped', 'N/A'),
            'python_is_reverse': python_details.get('is_reverse', 'N/A'),
            'python_mate_reverse': python_details.get('mate_reverse', 'N/A'),
            'python_is_read1': python_details.get('is_read1', 'N/A'),
            'python_is_read2': python_details.get('is_read2', 'N/A'),
            'python_is_secondary': python_details.get('is_secondary', 'N/A'),
            'python_is_qcfail': python_details.get('is_qcfail', 'N/A'),
            'python_is_duplicate': python_details.get('is_duplicate', 'N/A'),
            'python_is_supplementary': python_details.get('is_supplementary', 'N/A'),
        })
        for tag in all_tags:
            row[f'matlab_tag_{tag}'] = 'N/A'
            row[f'python_tag_{tag}'] = (python_details.get('tags') or {}).get(tag, 'N/A')
        diff_reads.append(row)
    
    if diff_reads:
        df_diff = pd.DataFrame(diff_reads)
        df_diff = df_diff.sort_values(['read_id', 'junction_type'])
        df_diff.to_csv(count_output_dir / 'differentially_labelled_reads.csv', index=False)
        log(f"  Saved {len(diff_reads)} differentially labelled reads to {count_output_dir / 'differentially_labelled_reads.csv'}", verbose)
    else:
        log("  No differentially labelled reads found!", verbose)

    # Per-junction read lists
    per_junction_rows = []
    for (read_id, jct_type), cls_list in matlab_classifications.items():
        details = matlab_classification_details.get((read_id, jct_type), {})
        for cls in cls_list:
            per_junction_rows.append({
                'junction_type': jct_type,
                'source': 'matlab',
                'read_id': read_id,
                'classification': cls if cls else 'N/A',
                'start': details.get('start', 'N/A'),
                'end': details.get('end', 'N/A'),
                'length': details.get('length', 'N/A'),
                'cigar': details.get('cigar', 'N/A'),
                'flag': details.get('flag', 'N/A'),
            })
    for (read_id, jct_type), cls_list in python_classifications.items():
        details = python_classification_details.get((read_id, jct_type), {})
        for cls in cls_list:
            per_junction_rows.append({
                'junction_type': jct_type,
                'source': 'python',
                'read_id': read_id,
                'classification': cls if cls else 'N/A',
                'start': details.get('start', 'N/A'),
                'end': details.get('end', 'N/A'),
                'length': details.get('length', 'N/A'),
                'cigar': details.get('cigar', 'N/A'),
                'flag': details.get('flag', 'N/A'),
            })
    if per_junction_rows:
        df_per_junction = pd.DataFrame(per_junction_rows)
        df_per_junction = df_per_junction.sort_values(['junction_type', 'source', 'read_id'])
        per_junction_path = count_output_dir / 'per_junction_read_lists.csv'
        df_per_junction.to_csv(per_junction_path, index=False)
        log(f"  Saved per-junction read lists to {per_junction_path}", verbose)
    
    return {
        'total_differences': len(diff_reads),
        'common_reads': len(common_keys),
        'matlab_only': len(matlab_only),
        'python_only': len(python_only),
    }


def generate_readme(output_dir: Path, args: argparse.Namespace) -> None:
    """Generate README.md in output directory."""
    readme_path = output_dir / 'README.md'
    
    with open(readme_path, 'w') as f:
        f.write("# Junction Alignment Comparison Results\n\n")
        f.write("This directory contains the results of comparing MATLAB and Python ")
        f.write("junction alignment outputs.\n\n")
        f.write("## Directory Structure\n\n")
        f.write("- `bam_comparison/`: Part 1 results (BAM file comparison)\n")
        f.write("- `read_count_comparison/`: Part 2 results (read count comparison)\n")
        f.write("- `comparison_summary.txt`: Overall summary\n")
        f.write("- `README.md`: This file\n\n")
        f.write("## Key Files\n\n")
        f.write("### Quick Start\n")
        f.write("1. Start with `comparison_summary.txt` for overall differences\n")
        f.write("2. Check `read_count_comparison/read_count_comparison_summary.csv` for count differences\n")
        f.write("3. See `read_count_comparison/differentially_labelled_reads.csv` for specific reads\n\n")
        f.write("### Part 1: BAM Comparison\n")
        f.write("- `bam_comparison/comparison_summary.txt`: Summary of alignment differences\n")
        f.write("- `bam_comparison/bam_comparison_differences.csv`: Detailed CSV of all differing alignments\n")
        f.write("- `bam_comparison/matlab_alignment.sam`: SAM conversion of MATLAB BAM\n")
        f.write("- `bam_comparison/python_alignment.sam`: SAM conversion of Python BAM\n\n")
        f.write("### Part 2: Read Count Comparison\n")
        f.write("- `read_count_comparison/read_count_comparison_summary.csv`: Side-by-side count comparison\n")
        f.write("- `read_count_comparison/differentially_labelled_reads.csv`: Reads with different classifications\n\n")
        f.write("## Interpretation Guide\n\n")
        f.write("### Count Differences\n")
        f.write("- Look at `*_diff` columns in `read_count_comparison_summary.csv`\n")
        f.write("- Zero = perfect match\n")
        f.write("- Non-zero = difference (positive = Python has more, negative = MATLAB has more)\n\n")
        f.write("### Junction Type Mapping\n")
        f.write("MATLAB uses numeric IDs (1-7), Python uses descriptive names:\n")
        for matlab_id, python_name in sorted(MATLAB_TO_PYTHON_JUNCTION_MAP.items()):
            f.write(f"- MATLAB {matlab_id} <-> Python {python_name}\n")
        f.write("\n## Command Used\n\n")
        f.write("```bash\n")
        f.write("python compare_jc_alignment.py \\\n")
        f.write(f"    --matlab-bam {args.matlab_bam} \\\n")
        f.write(f"    --python-bam {args.python_bam} \\\n")
        f.write(f"    --matlab-fasta {args.matlab_fasta} \\\n")
        f.write(f"    --python-fasta {args.python_fasta} \\\n")
        f.write(f"    --matlab-data-dir {args.matlab_data_dir} \\\n")
        f.write(f"    --output-dir {args.output_dir}\n")
        f.write("```\n")


def main():
    """Main entry point."""
    args = parse_args()
    
    # Validate inputs
    for path, name in [
        (args.matlab_bam, 'MATLAB BAM'),
        (args.python_bam, 'Python BAM'),
        (args.matlab_fasta, 'MATLAB FASTA'),
        (args.python_fasta, 'Python FASTA'),
        (args.matlab_data_dir, 'MATLAB data directory'),
    ]:
        if not path.exists():
            print(f"ERROR: {name} not found: {path}", file=sys.stderr)
            sys.exit(1)
    
    # Create output directory
    ensure_dir(args.output_dir)
    
    # Get junction length
    if args.junction_length is None:
        log("Auto-detecting junction length from FASTA...", args.verbose)
        args.junction_length = get_junction_length(args.matlab_fasta)
        log(f"  Detected junction length: {args.junction_length} bp", args.verbose)
    
    # Run comparisons
    bam_results = {}
    count_results = {}
    
    if not args.skip_bam_comparison:
        try:
            bam_results = compare_bam_files(
                args.matlab_bam, args.python_bam, args.output_dir, args.verbose
            )
        except Exception as e:
            print(f"ERROR in BAM comparison: {e}", file=sys.stderr)
            if args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
    
    if not args.skip_read_counts:
        try:
            count_results = compare_read_counts(
                args.matlab_bam, args.python_bam, args.matlab_data_dir,
                args.output_dir, args.junction_length, args.read_length,
                args.min_overlap_len, args.read_length_tolerance, args.min_bp_in_frame,
                args.verbose
            )
        except Exception as e:
            print(f"ERROR in read count comparison: {e}", file=sys.stderr)
            if args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
    
    # Generate overall summary
    log("\n=== Generating Summary ===", args.verbose)
    summary_path = args.output_dir / 'comparison_summary.txt'
    with open(summary_path, 'w') as f:
        f.write("Junction Alignment Comparison Summary\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"MATLAB BAM: {args.matlab_bam}\n")
        f.write(f"Python BAM: {args.python_bam}\n")
        f.write(f"Junction Length: {args.junction_length} bp\n")
        f.write(f"Read Length: {args.read_length} bp\n")
        f.write(f"Min Overlap Length: {args.min_overlap_len} bp\n\n")
        
        if bam_results:
            f.write("Part 1: BAM Comparison\n")
            f.write("-" * 30 + "\n")
            f.write(f"MATLAB unique reads: {bam_results.get('matlab_reads', 0)}\n")
            f.write(f"Python unique reads: {bam_results.get('python_reads', 0)}\n")
            f.write(f"Matched alignments: {bam_results.get('matched', 0)}\n")
            f.write(f"MATLAB-only: {bam_results.get('matlab_only', 0)}\n")
            f.write(f"Python-only: {bam_results.get('python_only', 0)}\n")
            f.write(f"Differing: {bam_results.get('differing', 0)}\n\n")
        
        if count_results:
            f.write("Part 2: Read Count Comparison\n")
            f.write("-" * 30 + "\n")
            f.write(f"Total differences: {count_results.get('total_differences', 0)}\n")
            f.write(f"Common read-junction pairs: {count_results.get('common_reads', 0)}\n")
            f.write(f"MATLAB-only: {count_results.get('matlab_only', 0)}\n")
            f.write(f"Python-only: {count_results.get('python_only', 0)}\n\n")
        
        # Check for differences
        has_differences = (
            (bam_results.get('matlab_only', 0) > 0 or
             bam_results.get('python_only', 0) > 0 or
             bam_results.get('differing', 0) > 0) or
            count_results.get('total_differences', 0) > 0
        )
        
        if has_differences:
            f.write("⚠️  DIFFERENCES FOUND\n")
            f.write("See detailed files in subdirectories for more information.\n")
        else:
            f.write("✅ NO DIFFERENCES FOUND\n")
            f.write("MATLAB and Python results are identical!\n")
    
    log(f"Summary saved to {summary_path}", args.verbose)
    
    # Generate README
    generate_readme(args.output_dir, args)
    log(f"README generated at {args.output_dir / 'README.md'}", args.verbose)
    
    log("\n=== Comparison Complete ===", args.verbose)
    log(f"Results saved to: {args.output_dir}", args.verbose)
    
    # Exit with appropriate code
    has_differences = (
        (bam_results.get('matlab_only', 0) > 0 or
         bam_results.get('python_only', 0) > 0 or
         bam_results.get('differing', 0) > 0) or
        count_results.get('total_differences', 0) > 0
    )
    sys.exit(1 if has_differences else 0)


if __name__ == '__main__':
    main()

"""Utilities for comparing junction coverage between MATLAB and Python."""

import pandas as pd
import pysam
from pathlib import Path


def load_read_buckets_from_fasta(alignment_dir: Path, file_ext: str = "fasta") -> dict[str, tuple[int, str]]:
    """Load FASTQ/FASTA files and map read_id/index to (jct_num, side).

    Args:
        alignment_dir: Directory containing reads_jct_*_*.fasta/fastq files
        file_ext: File extension ('fasta' or 'fastq')

    Returns dict mapping "read_name/index" -> (jct_num, side)
    where side is "left", "right", or "span"
    """
    buckets = {}
    for j in range(1, 8):
        for side in ("left", "right", "green"):
            fasta_path = alignment_dir / f"reads_jct_{j}_{side}.{file_ext}"
            with open(fasta_path, "r") as f:
                for line in f:
                    if line.startswith(">"):  # FASTA format
                        # Format: >read_name/index
                        read_id = line.strip()[1:]  # Remove >
                        buckets[read_id] = (j, "span" if side == "green" else side)
                    elif line.startswith("@"):  # FASTQ format
                        # Format: @read_name/index
                        read_id = line.strip()[1:]  # Remove @
                        buckets[read_id] = (j, "span" if side == "green" else side)
    return buckets


def load_matlab_read_buckets(matlab_alignment_dir: Path) -> dict[str, tuple[int, str]]:
    """Load MATLAB FASTQ files and map read_id/index to (jct_num, side).

    Returns dict mapping "read_name/index" -> (jct_num, side)
    where side is "left", "right", or "span"
    """
    return load_read_buckets_from_fasta(matlab_alignment_dir, file_ext="fastq")


def load_python_read_buckets(python_alignment_dir: Path) -> dict[str, tuple[int, str]]:
    """Load Python .bam_indices files and map read_id to (jct_num, side).

    Reads files written by write_junction_read_bam_indices (DEBUG mode).
    Format per line: read_id\\tstart\\tend\\tis_reverse\\tbam_index

    Returns dict mapping "read_name/bam_index" -> (jct_num, side)
    where side is "left", "right", or "span"
    """
    # Map ReadType value names to the bucket side names used by comparison
    _SIDE_MAP = {
        "left": "left",
        "left_marginal": "left",
        "right": "right",
        "right_marginal": "right",
        "spanning": "span",
        "paired": "span",
    }
    buckets: dict[str, tuple[int, str]] = {}
    for j in range(1, 8):
        for read_type_value, side in _SIDE_MAP.items():
            indices_path = python_alignment_dir / f"reads_jct_{j}_{read_type_value}.bam_indices"
            if not indices_path.exists():
                continue
            with open(indices_path, "r") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    read_id = parts[0]
                    bam_index = parts[4]
                    key = f"{read_id}/{bam_index}"
                    buckets[key] = (j, side)
    return buckets


def create_jct_comparison_table(
    python_bam: Path,
    matlab_buckets: dict[str, tuple[int, str]],
    python_buckets: dict[str, tuple[int, str]],
    jct_num: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create comparison table for a single junction.

    Args:
        python_bam: Path to Python BAM file
        matlab_buckets: MATLAB read buckets from FASTQ files
        python_buckets: Python read buckets from FASTA files
        jct_num: Junction number (1-7)

    Returns tuple of:
    - DataFrame with columns (only rows where matlab_type != python_type):
      - read_index: 1-based BAM index
      - read_name: Read name
      - start: Start position (1-based)
      - end: End position (1-based)
      - matlab_type: Classification from MATLAB ('left', 'right', 'span', or None)
      - python_type: Classification from Python ('left', 'right', 'span', or None)
    - Confusion matrix DataFrame (4x4) with rows/cols: left, span, right, none
    """
    rows = []

    # Initialize confusion matrix counts
    categories = ['left', 'span', 'right', 'none']
    confusion_counts = {matlab: {python: 0 for python in categories} for matlab in categories}

    with pysam.AlignmentFile(str(python_bam), "rb") as bam:
        # Iterate through BAM once, tracking 1-based index
        for read_index, read in enumerate(bam.fetch(until_eof=True), start=1):
            read_name = read.query_name
            start = read.reference_start + 1
            end = read.reference_end

            read_id_key = f"{read_name}/{read_index}"

            # Check MATLAB classification
            matlab_type = None
            if read_id_key in matlab_buckets:
                matlab_jct, matlab_side = matlab_buckets[read_id_key]
                if matlab_jct == jct_num:
                    matlab_type = matlab_side

            # Check Python classification
            python_type = None
            if read_id_key in python_buckets:
                python_jct, python_side = python_buckets[read_id_key]
                if python_jct == jct_num:
                    python_type = python_side

            # Normalize None to "none" for confusion matrix
            matlab_key = matlab_type if matlab_type is not None else "none"
            python_key = python_type if python_type is not None else "none"

            # Update confusion matrix
            confusion_counts[matlab_key][python_key] += 1

            # NOTE: Change this condition to indicate the disagreement that you want to see.
            if matlab_type != python_type and (matlab_type == "span" or python_type == "span"):
                # Only include rows where matlab_type and python_type differ
                rows.append({
                    "read_index": read_index,
                    "read_name": read_name,
                    "start": start,
                    "end": end,
                    "matlab_type": matlab_type,
                    "python_type": python_type,
                })

    # Create confusion matrix DataFrame
    confusion_matrix = pd.DataFrame(
        [[confusion_counts[matlab][python] for python in categories] for matlab in categories],
        index=categories,
        columns=categories
    )

    return pd.DataFrame(rows), confusion_matrix

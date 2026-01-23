"""Bowtie2 wrapper for read alignment."""

from pathlib import Path
from typing import Optional, Union, Tuple

from amplifinder.env import SAMTOOLS_PATH, BOWTIE2_PATH
from amplifinder.utils.run_utils import get_tool_path, run_command
from amplifinder.utils.file_utils import ensure_parent_dir
from amplifinder.utils.timing import print_timer


def run_bowtie2_build(ref_fasta: Path, index_prefix: Path) -> None:
    """Build bowtie2 index from FASTA file.

    Args:
        ref_fasta: Path to reference FASTA file
        index_prefix: Prefix for output index files (e.g., /path/to/index)

    Raises:
        FileNotFoundError: If bowtie2-build not found or ref_fasta doesn't exist
        RuntimeError: If bowtie2-build fails
    """
    bowtie2_build = get_tool_path("bowtie2-build", config_path=BOWTIE2_PATH, required=True)

    if not ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")

    # Create output directory if needed
    ensure_parent_dir(index_prefix)

    cmd = [
        bowtie2_build,
        "--quiet",
        str(ref_fasta),
        str(index_prefix),
    ]

    run_command(cmd, check=True, capture_output=True, text=True)


def _format_score_min(score_min: Union[str, Tuple[float, float]], local: bool) -> str:
    """Format score_min for bowtie2."""
    if isinstance(score_min, str):
        return score_min
    mode = "G" if local else "L"
    return f"{mode},{score_min[0]},{score_min[1]}"


def _format_mismatch_penalty(mismatch_penalty: Union[str, Tuple[int, int]]) -> str:
    """Format mismatch penalty tuple or string as a single bowtie2 argument."""
    if isinstance(mismatch_penalty, str):
        return mismatch_penalty
    return f"{mismatch_penalty[0]},{mismatch_penalty[1]}"


def run_bowtie2_align(
    index_prefix: Path,
    fastq_path: Path,
    output_sam: Path,
    # bowtie scoring/behavior
    score_min: Tuple[float, float] = (0, -0.25),
    mismatch_penalty: Union[str, Tuple[int, int]] = (5, 5),
    local: bool = True,
    num_alignments: int = 100,
    # resources
    threads: int = 4,
) -> None:
    """Align reads to index using bowtie2.

    Args:
        index_prefix: Path to bowtie2 index prefix
        fastq_path: Path to FASTQ file (or directory with FASTQ files)
        output_sam: Path to output SAM file
        score_min: Minimum alignment score function. If None, uses:
                   - "G,0,-0.25" for local alignment (matches MATLAB)
                   - "L,0,-0.25" for end-to-end alignment (matches MATLAB)
        num_alignments: Maximum alignments to report per read (-k)
        threads: Number of threads to use
        local: Use local alignment (default True for junction alignment)

    Raises:
        FileNotFoundError: If bowtie2 not found or inputs don't exist
        RuntimeError: If bowtie2 fails
    """
    bowtie2 = get_tool_path("bowtie2", config_path=BOWTIE2_PATH, required=True)

    fastq_path = Path(fastq_path)
    if not fastq_path.exists():
        raise FileNotFoundError(f"FASTQ not found: {fastq_path}")

    # Handle directory or file input
    if fastq_path.is_dir():
        # Find FASTQ files in directory
        fastq_files = list(fastq_path.glob("*.fastq")) + list(fastq_path.glob("*.fastq.gz"))
        fastq_files += list(fastq_path.glob("*.fq")) + list(fastq_path.glob("*.fq.gz"))
        if not fastq_files:
            raise FileNotFoundError(f"No FASTQ files found in: {fastq_path}")
        reads_arg = ",".join(str(f) for f in sorted(fastq_files))
    else:
        reads_arg = str(fastq_path)

    # Create output directory if needed
    ensure_parent_dir(output_sam)

    # Set score_min based on alignment mode (or override)
    score_min_str = _format_score_min(score_min, local)
    mp = _format_mismatch_penalty(mismatch_penalty)
    cmd = [
        bowtie2,
        "-x", str(index_prefix),
        "-U", reads_arg,
        "-S", str(output_sam),
        "-k", str(num_alignments),
        "--score-min", score_min_str,
        "--mp", mp,
        "-p", str(threads),
        "--reorder",
    ]

    if local:
        cmd.append("--local")
    else:
        cmd.append("--end-to-end")  # Explicitly match MATLAB's --end-to-end flag

    run_command(cmd, check=True, capture_output=True, text=True)


def sam_to_sorted_bam(
    sam_path: Path,
    bam_path: Path,
    samtools_path: Optional[Path] = None,
    threads: int = 1,
    min_qlen: Optional[int] = None,
) -> None:
    """Convert SAM to sorted BAM and index.

    Args:
        sam_path: Input SAM file
        bam_path: Output BAM file (will also create .bai index)
        samtools_path: Path to samtools (auto-detect if None)
        threads: Number of threads for sorting

    Raises:
        FileNotFoundError: If samtools not found or SAM doesn't exist
        RuntimeError: If samtools fails
    """
    samtools = get_tool_path("samtools", config_path=samtools_path or SAMTOOLS_PATH)

    if not sam_path.exists():
        raise FileNotFoundError(f"SAM file not found: {sam_path}")

    # Convert and sort (optionally filter by min_qlen)
    input_for_sort = sam_path
    if min_qlen is not None:
        temp_bam = bam_path.with_suffix(".filtered.bam")
        cmd_view = [
            samtools, "view",
            "-bS",
            "--min-qlen", str(min_qlen),
            "-@", str(threads),
            "-o", str(temp_bam),
            str(sam_path),
        ]
        run_command(cmd_view, check=True, capture_output=True, text=True)
        input_for_sort = temp_bam

    cmd_sort = [
        samtools, "sort",
        "-@", str(threads),
        "-o", str(bam_path),
        str(input_for_sort),
    ]
    run_command(cmd_sort, check=True, capture_output=True, text=True)

    if min_qlen is not None and temp_bam.exists():
        temp_bam.unlink()

    # Index
    cmd_index = [samtools, "index", str(bam_path)]

    run_command(cmd_index, check=True, capture_output=True, text=True)


def align_reads_to_fasta(
    ref_fasta: Path,
    fastq_path: Path,
    output_bam: Path,
    threads: int = 1,
    score_min: Union[str, Tuple[float, float]] = (0, -0.25),
    mismatch_penalty: Union[str, Tuple[int, int]] = (5, 5),
    num_alignments: int = 100,
    local: bool = False,  # Default: end-to-end (matches MATLAB)
    min_qlen: Optional[int] = None,
    keep_sam: bool = False,
) -> None:
    """Full alignment pipeline: build index, align, convert to sorted BAM.

    Args:
        ref_fasta: Reference FASTA file
        fastq_path: FASTQ file or directory
        output_bam: Output sorted BAM file
        threads: Number of threads
        score_min: Alignment score function (tuple uses mode G/L based on local).
        mismatch_penalty: Bowtie2 --mp value (tuple or string).
        num_alignments: Max alignments per read (-k).
        local: Use local alignment (default False = end-to-end).
        keep_sam: Keep intermediate SAM file (default False).
    """
    # Paths
    work_dir = output_bam.parent
    index_prefix = work_dir / "bowtie2_index"
    sam_path = work_dir / "alignments.sam"

    # Build index
    with print_timer("Indexing... ", end_msg=", ", seperate_prints=True):
        run_bowtie2_build(ref_fasta, index_prefix)

    # Align
    with print_timer("Aligning... ", end_msg=", ", seperate_prints=True):
        run_bowtie2_align(
            index_prefix=index_prefix,
            fastq_path=fastq_path,
            output_sam=sam_path,
            score_min=score_min,
            mismatch_penalty=mismatch_penalty,
            local=local,
            num_alignments=num_alignments,
            threads=threads,
        )

    # Convert to sorted BAM
    with print_timer("BAM... ", end_msg="", seperate_prints=True):
        sam_to_sorted_bam(sam_path, output_bam, threads=threads, min_qlen=min_qlen)
    print()  # Newline after final step

    # Cleanup
    if not keep_sam:
        sam_path.unlink()

    # Cleanup index files
    for idx_file in work_dir.glob("bowtie2_index.*"):
        idx_file.unlink()

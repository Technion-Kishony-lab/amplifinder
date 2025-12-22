"""Bowtie2 wrapper for read alignment."""

import subprocess
import shutil
from pathlib import Path
from typing import Optional

from amplifinder.logger import info, debug
from amplifinder.env import SAMTOOLS_PATH


def get_bowtie2_path() -> Optional[Path]:
    """Find bowtie2 executable in PATH."""
    path = shutil.which("bowtie2")
    return Path(path) if path else None


def get_bowtie2_build_path() -> Optional[Path]:
    """Find bowtie2-build executable in PATH."""
    path = shutil.which("bowtie2-build")
    return Path(path) if path else None


def check_bowtie2_available() -> bool:
    """Check if bowtie2 and bowtie2-build are available."""
    return get_bowtie2_path() is not None and get_bowtie2_build_path() is not None


def run_bowtie2_build(ref_fasta: Path, index_prefix: Path) -> None:
    """Build bowtie2 index from FASTA file.
    
    Args:
        ref_fasta: Path to reference FASTA file
        index_prefix: Prefix for output index files (e.g., /path/to/index)
    
    Raises:
        FileNotFoundError: If bowtie2-build not found or ref_fasta doesn't exist
        RuntimeError: If bowtie2-build fails
    """
    bowtie2_build = get_bowtie2_build_path()
    if not bowtie2_build:
        raise FileNotFoundError("bowtie2-build not found in PATH")
    
    if not ref_fasta.exists():
        raise FileNotFoundError(f"Reference FASTA not found: {ref_fasta}")
    
    # Create output directory if needed
    index_prefix.parent.mkdir(parents=True, exist_ok=True)
    
    cmd = [
        str(bowtie2_build),
        "--quiet",
        str(ref_fasta),
        str(index_prefix),
    ]
    
    debug(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"bowtie2-build failed: {result.stderr}")
    
    info(f"Built bowtie2 index: {index_prefix}")


def run_bowtie2_align(
    index_prefix: Path,
    fastq_path: Path,
    output_sam: Path,
    score_min: Optional[str] = None,
    num_alignments: int = 10,
    threads: int = 1,
    local: bool = True,
) -> None:
    """Align reads to index using bowtie2.
    
    Args:
        index_prefix: Path to bowtie2 index prefix
        fastq_path: Path to FASTQ file (or directory with FASTQ files)
        output_sam: Path to output SAM file
        score_min: Minimum alignment score function. If None, uses:
                   - "G,0,-0.25" for local alignment (matches MATLAB)
                   - "L,0,-0.6" for end-to-end alignment
        num_alignments: Maximum alignments to report per read (-k)
        threads: Number of threads to use
        local: Use local alignment (default True for junction alignment)
    
    Raises:
        FileNotFoundError: If bowtie2 not found or inputs don't exist
        RuntimeError: If bowtie2 fails
    """
    bowtie2 = get_bowtie2_path()
    if not bowtie2:
        raise FileNotFoundError("bowtie2 not found in PATH")
    
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
    output_sam.parent.mkdir(parents=True, exist_ok=True)
    
    # Set default score_min based on alignment mode
    # MATLAB uses end-to-end with L,0,-0.6, local with G,0,-0.25
    if score_min is None:
        score_min = "G,0,-0.25" if local else "L,0,-0.6"
    
    cmd = [
        str(bowtie2),
        "-x", str(index_prefix),
        "-U", reads_arg,
        "-S", str(output_sam),
        "-k", str(num_alignments),
        "--score-min", score_min,
        "-p", str(threads),
    ]
    
    # Add mismatch penalty to match MATLAB (--mp 5,5)
    # MATLAB uses this for both end-to-end and local modes
    cmd.extend(["--mp", "5,5"])
    
    if local:
        cmd.append("--local")
    
    debug(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"bowtie2 failed: {result.stderr}")
    
    info(f"Aligned reads to: {output_sam}")


def sam_to_sorted_bam(
    sam_path: Path,
    bam_path: Path,
    samtools_path: Optional[Path] = None,
    threads: int = 1,
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
    if samtools_path is None:
        # Try config path first
        if SAMTOOLS_PATH:
            samtools_path = SAMTOOLS_PATH
        else:
            # Try PATH as last resort
            samtools_path = shutil.which("samtools")
            if not samtools_path:
                raise FileNotFoundError(
                    "samtools not found. Set samtools_path in config.yaml or ensure samtools is in PATH. "
                    f"Config path: {SAMTOOLS_PATH}"
                )
    
    samtools_path = Path(samtools_path)
    
    # If path is a directory, append "samtools" to the path
    if samtools_path.is_dir():
        samtools_path = samtools_path / "samtools"
    
    # Verify the final path exists
    if not samtools_path.exists():
        raise FileNotFoundError(f"samtools not found at: {samtools_path}")
    
    if not sam_path.exists():
        raise FileNotFoundError(f"SAM file not found: {sam_path}")
    
    # Convert and sort
    cmd_sort = [
        str(samtools_path), "sort",
        "-@", str(threads),
        "-o", str(bam_path),
        str(sam_path),
    ]
    
    debug(f"Running: {' '.join(cmd_sort)}")
    result = subprocess.run(cmd_sort, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"samtools sort failed: {result.stderr}")
    
    # Index
    cmd_index = [str(samtools_path), "index", str(bam_path)]
    
    debug(f"Running: {' '.join(cmd_index)}")
    result = subprocess.run(cmd_index, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"samtools index failed: {result.stderr}")
    
    info(f"Created sorted BAM: {bam_path}")


def align_reads_to_fasta(
    ref_fasta: Path,
    fastq_path: Path,
    output_bam: Path,
    threads: int = 1,
    score_min: Optional[str] = None,
    num_alignments: int = 10,
    local: bool = False,  # MATLAB uses --end-to-end by default
    keep_sam: bool = False,
) -> None:
    """Full alignment pipeline: build index, align, convert to sorted BAM.
    
    Args:
        ref_fasta: Reference FASTA file
        fastq_path: FASTQ file or directory
        output_bam: Output sorted BAM file
        threads: Number of threads
        score_min: Minimum alignment score function. If None, uses appropriate default.
        num_alignments: Max alignments per read
        local: Use local alignment (default True for junction alignment)
        keep_sam: Keep intermediate SAM file (default False)
    """
    # Paths
    work_dir = output_bam.parent
    index_prefix = work_dir / "bowtie2_index"
    sam_path = work_dir / "alignments.sam"
    
    # Build index
    run_bowtie2_build(ref_fasta, index_prefix)
    
    # Align
    run_bowtie2_align(
        index_prefix=index_prefix,
        fastq_path=fastq_path,
        output_sam=sam_path,
        score_min=score_min,
        num_alignments=num_alignments,
        threads=threads,
        local=local,
    )
    
    # Convert to sorted BAM
    sam_to_sorted_bam(sam_path, output_bam, threads=threads)
    
    # Cleanup
    if not keep_sam:
        sam_path.unlink()
    
    # Cleanup index files
    for idx_file in work_dir.glob("bowtie2_index.*"):
        idx_file.unlink()


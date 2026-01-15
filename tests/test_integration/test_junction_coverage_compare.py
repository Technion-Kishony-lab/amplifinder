"""Compare junction coverage and read labels between MATLAB and Python."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set, Tuple

import csv
import pytest
import pysam

from amplifinder.data_types import Side, JunctionType
from amplifinder.tools.bowtie2 import align_reads_to_fasta
from amplifinder.utils.fasta import get_read_length_stats
from amplifinder.utils.file_utils import ensure_dir


pytestmark = pytest.mark.integration
REQUIRE_MATLAB_FILES = True

MATLAB_JUNCTION_MAPPING = {
    "1": JunctionType.CHR_TO_AMP_LEFT,
    "2": JunctionType.CHR_TO_TN_LEFT,
    "3": JunctionType.AMP_RIGHT_TO_TN_LEFT,
    "4": JunctionType.AMP_RIGHT_TO_AMP_LEFT,
    "5": JunctionType.TN_RIGHT_TO_AMP_LEFT,
    "6": JunctionType.TN_RIGHT_TO_CHR,
    "7": JunctionType.AMP_RIGHT_TO_CHR,
}


def _normalize_junction_name(name: str) -> str:
    """Normalize MATLAB numeric junction names to Python JunctionType names."""
    if name in MATLAB_JUNCTION_MAPPING:
        return MATLAB_JUNCTION_MAPPING[name].name
    if name in JunctionType.__members__:
        return name
    return name


@dataclass
class ReadSummary:
    green: bool
    alignments: Set[Tuple[str, int, int, str]]
    green_alignments: Set[Tuple[str, int, int, str]]
    flags: Set[int]
    cigars: Set[str]
    seq: Optional[str] = None
    qual: Optional[List[int]] = None


def _is_cigar_only_m(read: pysam.AlignedSegment) -> bool:
    if not read.cigartuples:
        return False
    return all(op == 0 for op, _ in read.cigartuples)


def _classify_read(
    read: pysam.AlignedSegment,
    jct_len: int,
    avg_read_length: int,
    min_overlap_len: int,
    read_length_tolerance: float = 0.1,
    min_bp_in_frame: int = 10,
) -> Optional[Side]:
    """MATLAB-equivalent read classification (left/right/green/undetermined)."""
    alignment_length = read.query_alignment_length
    read_length_factor = 1 + read_length_tolerance
    min_alignment_length = avg_read_length / read_length_factor
    max_alignment_length = avg_read_length * read_length_factor
    if not (min_alignment_length <= alignment_length <= max_alignment_length):
        return None

    if not _is_cigar_only_m(read):
        return None

    # Filter by alignment tags (match pipeline)
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

    start_1 = read.reference_start + 1
    end_1 = read.reference_end
    jct_point = jct_len / 2

    okz_left = (start_1 + alignment_length > jct_point - avg_read_length + min_bp_in_frame) and \
               (end_1 < jct_point + min_overlap_len)
    okz_right = (start_1 > jct_point - min_overlap_len) and \
                (start_1 < jct_point + avg_read_length - min_bp_in_frame)
    okz_green = (start_1 < jct_point - min_overlap_len) and \
                (end_1 > jct_point + min_overlap_len)

    if okz_green:
        return Side.MIDDLE
    if okz_left:
        return Side.LEFT
    if okz_right:
        return Side.RIGHT
    return None


def _collect_counts_and_reads(
    bam_path: Path,
    avg_read_length: int,
    min_overlap_len: int,
) -> Tuple[Dict[str, Dict[str, int]], Dict[str, ReadSummary]]:
    """Collect per-junction counts and per-read summaries from BAM."""
    counts: Dict[str, Dict[str, int]] = {}
    reads: Dict[str, ReadSummary] = {}

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        jct_lengths = dict(zip(bam.references, bam.lengths))
        for read in bam.fetch():
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            jct_name_raw = read.reference_name
            jct_len = jct_lengths[jct_name_raw]
            jct_name = _normalize_junction_name(jct_name_raw)
            read_type = _classify_read(
                read=read,
                jct_len=jct_len,
                avg_read_length=avg_read_length,
                min_overlap_len=min_overlap_len,
            )

            if jct_name not in counts:
                counts[jct_name] = {"left": 0, "right": 0, "spanning": 0}
            if read_type == Side.LEFT:
                counts[jct_name]["left"] += 1
            elif read_type == Side.RIGHT:
                counts[jct_name]["right"] += 1
            elif read_type == Side.MIDDLE:
                counts[jct_name]["spanning"] += 1

            read_name = read.query_name
            if read_name not in reads:
                reads[read_name] = ReadSummary(
                    green=False,
                    alignments=set(),
                    green_alignments=set(),
                    flags=set(),
                    cigars=set(),
                    seq=read.query_sequence,
                    qual=read.query_qualities,
                )
            summary = reads[read_name]
            alignment_tuple = (
                jct_name,
                read.reference_start,
                read.reference_end,
                read.cigarstring or "",
            )
            if read_type == Side.MIDDLE:
                summary.green = True
                summary.green_alignments.add(alignment_tuple)
            summary.alignments.add(alignment_tuple)
            summary.flags.add(read.flag)
            summary.cigars.add(read.cigarstring or "")

    return counts, reads


def _ensure_python_bam(
    jc_fasta: Path,
    fastq_dir: Path,
    output_dir: Path,
) -> Path:
    """Run Python alignment against MATLAB jc.fasta if needed."""
    output_dir = ensure_dir(output_dir)
    bam_path = output_dir / "sorted.bam"
    if bam_path.exists():
        return bam_path

    align_reads_to_fasta(
        ref_fasta=jc_fasta,
        fastq_path=fastq_dir,
        output_bam=bam_path,
        threads=4,
        score_min=None,  # use default (matches MATLAB)
        num_alignments=100,
        local=False,
        keep_sam=False,
    )
    return bam_path


def _write_fastq(
    path: Path,
    python_reads: Dict[str, ReadSummary],
    matlab_reads: Dict[str, ReadSummary],
    names: Iterable[str],
) -> int:
    """Write reads to FASTQ. Returns number of reads written."""
    count = 0
    with path.open("w") as handle:
        for name in names:
            summary = python_reads.get(name) or matlab_reads.get(name)
            if summary is None or summary.seq is None or summary.qual is None:
                continue
            qual_str = "".join(chr(q + 33) for q in summary.qual)
            handle.write(f"@{name}\n{summary.seq}\n+\n{qual_str}\n")
            count += 1
    return count


@pytest.mark.slow
def test_compare_junction_coverage_matlab_python(isolate_srr25242877, tmp_path):
    """Compare MATLAB vs Python junction coverage and read labels."""
    matlab_candidate_dir = Path(
        "/zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/"
        "SRR25242877/chr_U00096_873161F_890745R_IS_41R"
    )
    matlab_bam = matlab_candidate_dir / "alignment" / "alignment.sorted.bam"
    matlab_jc_fasta = matlab_candidate_dir / "jc.fasta"

    if not matlab_bam.exists() or not matlab_jc_fasta.exists():
        if REQUIRE_MATLAB_FILES:
            pytest.fail(f"Missing MATLAB alignment files in: {matlab_candidate_dir}")
        pytest.skip("MATLAB alignment files not available")

    fastq_dir = isolate_srr25242877["fastq_path"]
    read_stats = get_read_length_stats(fastq_dir, sample_per_file=500)
    avg_read_length = int(round(read_stats.mean_length))
    min_overlap_len = 12

    python_output_dir = ensure_dir(
        Path(__file__).parent / "test_outputs" / "integration" /
        "matlab_compare" / isolate_srr25242877["iso_name"]
    )
    python_bam = _ensure_python_bam(
        jc_fasta=matlab_jc_fasta,
        fastq_dir=fastq_dir,
        output_dir=python_output_dir,
    )

    matlab_counts, matlab_reads = _collect_counts_and_reads(
        matlab_bam, avg_read_length=avg_read_length, min_overlap_len=min_overlap_len
    )
    python_counts, python_reads = _collect_counts_and_reads(
        python_bam, avg_read_length=avg_read_length, min_overlap_len=min_overlap_len
    )

    # Compare per-junction coverage
    all_jcts = sorted(set(matlab_counts) | set(python_counts))
    mismatches = []
    for jct in all_jcts:
        m = matlab_counts.get(jct, {"left": 0, "right": 0, "spanning": 0})
        p = python_counts.get(jct, {"left": 0, "right": 0, "spanning": 0})
        if m != p:
            mismatches.append((jct, m, p))

    # Identify reads that are green in one but not the other, and alignments differ
    read_names = set(matlab_reads) | set(python_reads)
    diff_green_reads: List[str] = []
    for name in sorted(read_names):
        m = matlab_reads.get(name)
        p = python_reads.get(name)
        m_green = m.green if m else False
        p_green = p.green if p else False
        if m_green == p_green:
            continue
        m_align = m.alignments if m else set()
        p_align = p.alignments if p else set()
        if m_align != p_align:
            diff_green_reads.append(name)

    # Write FASTQ for differing green reads
    fastq_out = python_output_dir / "green_mismatch_reads.fastq"
    fastq_written = _write_fastq(
        fastq_out, python_reads, matlab_reads, diff_green_reads
    )

    # Write table if <100 reads in FASTQ
    if 0 < fastq_written < 100:
        table_path = python_output_dir / "green_mismatch_reads.csv"
        with table_path.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow([
                "read_header",
                "python_flags",
                "matlab_flags",
                "python_cigars",
                "matlab_cigars",
                "python_alignments",
                "matlab_alignments",
                "python_green_alignments",
                "matlab_green_alignments",
            ])
            for name in diff_green_reads:
                m = matlab_reads.get(name)
                p = python_reads.get(name)
                p_flags = ";".join(str(f) for f in sorted(p.flags)) if p else ""
                m_flags = ";".join(str(f) for f in sorted(m.flags)) if m else ""
                p_cigars = ";".join(sorted(p.cigars)) if p else ""
                m_cigars = ";".join(sorted(m.cigars)) if m else ""
                p_alignments = ";".join(
                    f"{scaf}:{start}-{end}:{cigar}"
                    for scaf, start, end, cigar in sorted(p.alignments)
                ) if p else ""
                m_alignments = ";".join(
                    f"{scaf}:{start}-{end}:{cigar}"
                    for scaf, start, end, cigar in sorted(m.alignments)
                ) if m else ""
                p_green_alignments = ";".join(
                    f"{scaf}:{start}-{end}:{cigar}"
                    for scaf, start, end, cigar in sorted(p.green_alignments)
                ) if p else ""
                m_green_alignments = ";".join(
                    f"{scaf}:{start}-{end}:{cigar}"
                    for scaf, start, end, cigar in sorted(m.green_alignments)
                ) if m else ""
                writer.writerow([
                    name,
                    p_flags,
                    m_flags,
                    p_cigars,
                    m_cigars,
                    p_alignments,
                    m_alignments,
                    p_green_alignments,
                    m_green_alignments,
                ])

    if mismatches:
        mismatch_str = "\n".join(
            f"{jct}: matlab={m} python={p}" for jct, m, p in mismatches
        )
        pytest.fail(f"Junction coverage mismatches:\n{mismatch_str}")

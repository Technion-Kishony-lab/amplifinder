"""Export BAM indices for reads by junction and read type."""
from pathlib import Path

from amplifinder.steps.jct_coverage.alignment_data import AlignmentData
from amplifinder.data_types import JunctionType, ReadType


def load_junction_read_bam_indices(
    indices_dir: Path
) -> dict[JunctionType, dict[ReadType, list[int]]]:
    """Load BAM indices from files (reflection of write_junction_read_bam_indices).

    Args:
        indices_dir: Directory containing reads_jct_*_*.bam_indices files

    Returns:
        dict[JunctionType, dict[ReadType, list[int]]]
    """
    result: dict[JunctionType, dict[ReadType, list[int]]] = {}

    for jct_type in JunctionType:
        result[jct_type] = {}

        for read_type in ReadType:
            indices_path = indices_dir / f"reads_jct_{jct_type.num}_{read_type.value}.bam_indices"
            indices = []
            with open(indices_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    # Parse tab-separated: read_id, start, end, is_reverse, bam_index
                    parts = line.split('\t')
                    bam_idx = int(parts[4])
                    indices.append(bam_idx)

            result[jct_type][read_type] = indices

    return result


def extract_bam_indices_by_read_type(
    alignments_by_read_type: dict[ReadType, list[AlignmentData]]
) -> dict[ReadType, list[tuple[str, int, int, bool, int]]]:
    """Extract (read_id, start, end, is_reverse, bam_index) tuples by read type."""
    
    indices_by_read_type: dict[ReadType, list[tuple[str, int, int, bool, int]]] = {rt: [] for rt in ReadType}

    for read_type, alignments_list in alignments_by_read_type.items():
        for alignment in alignments_list:
            # Get all single alignments (handles both SingleAlignment and CombinedSingleAlignment)
            single_alignments = alignment.get_all_single_alignments()
            for single_aln in single_alignments:
                indices_by_read_type[read_type].append((
                    single_aln.read_id,
                    single_aln.start,
                    single_aln.end,
                    single_aln.is_reverse,
                    single_aln.bam_index
                ))

    return indices_by_read_type


def write_junction_read_bam_indices(
    alignment_data: dict[JunctionType, list[AlignmentData]],
    output_dir: Path
) -> None:
    """Write BAM indices with read_id, start, end, is_reverse for each read type by junction (1-7)."""

    output_dir = Path(output_dir)
    for jc_type, alignments in alignment_data.items():
        # Group alignments by read_type
        alignments_by_read_type = {rt: [a for a in alignments if a.read_type == rt] for rt in ReadType}
        indices_by_read_type = extract_bam_indices_by_read_type(alignments_by_read_type)

        # Write indices to text files (tab-separated: read_id, start, end, is_reverse, bam_index)
        for read_type, indices in indices_by_read_type.items():

            # Convert JunctionType and ReadType enum values to filename-safe string
            output_path = output_dir / f"reads_jct_{jc_type.num}_{read_type.value}.bam_indices"

            with open(output_path, 'w') as f:
                for read_id, start, end, is_reverse, bam_idx in indices:
                    f.write(f"{read_id}\t{start}\t{end}\t{is_reverse}\t{bam_idx}\n")

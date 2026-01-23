"""Export BAM indices for reads by junction and read type."""
from pathlib import Path

from amplifinder.steps.jct_coverage.alignment_data import AlignmentData
from amplifinder.data_types.enums import JunctionType, ReadType


def load_junction_read_bam_indices(
    indices_dir: Path
) -> dict[JunctionType, dict[ReadType, list[int] | list[tuple[int, int]]]]:
    """Load BAM indices from files (reflection of write_junction_read_bam_indices).

    Args:
        indices_dir: Directory containing reads_jct_*_*.bam_indices files

    Returns:
        dict[JunctionType, dict[ReadType, list[int] | list[tuple[int, int]]]]
    """
    result: dict[JunctionType, dict[ReadType, list[int] | list[tuple[int, int]]]] = {}

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

                    if read_type == ReadType.PAIRED:
                        parts = line.split()
                        indices.append((int(parts[0]), int(parts[1])))
                    else:
                        indices.append(int(line))

            result[jct_type][read_type] = indices

    return result


def extract_bam_indices_by_read_type(
    alignments: list[AlignmentData]
) -> dict[ReadType, list[int] | list[tuple[int, int]]]:

    indices_by_read_type: dict[ReadType, list[int] | list[tuple[int, int]]] = {rt: [] for rt in ReadType}

    for alignment in alignments:
        indices = alignment.get_bam_indices()
        indices = indices[0] if len(indices) == 1 else indices
        indices_by_read_type[alignment.read_type].append(indices)

    return indices_by_read_type


def write_junction_read_bam_indices(
    alignment_data: dict[JunctionType, list[AlignmentData]],
    output_dir: Path
) -> None:
    """Write BAM indices for each read type by junction (1-7)."""

    output_dir = Path(output_dir)
    for jc_type, alignments in alignment_data.items():
        indices_by_read_type = extract_bam_indices_by_read_type(alignments)

        # Write indices to text files (one per line, pairs for PAIRED)
        for read_type, indices in indices_by_read_type.items():

            # Convert JunctionType and ReadType enum values to filename-safe string
            output_path = output_dir / f"reads_jct_{jc_type.num}_{read_type.value}.bam_indices"

            with open(output_path, 'w') as f:
                for item in indices:
                    # Handle both single int and tuple[int, int] uniformly
                    if isinstance(item, tuple):
                        f.write(f"{item[0]} {item[1]}\n")
                    else:
                        f.write(f"{item}\n")

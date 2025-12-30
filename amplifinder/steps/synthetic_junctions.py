"""Step 10: Create synthetic junction sequences for alignment."""

from pathlib import Path
from typing import Optional

from Bio.Seq import reverse_complement

from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, Genome, JunctionType, RefTnLoc,
)
from amplifinder.steps.base import Step
from amplifinder.utils.file_utils import ensure_parent_dir


def create_synthetic_junctions(
    candidate: FilteredTnJc2,
    chr_seq: str,
    tn_loc: RefTnLoc,
    tn_seq: str,
    read_length: int,
    genome: 'Genome',
) -> dict[JunctionType, str]:
    """Create 7 synthetic junction sequences for a candidate amplicon.

    Junction types (amplicon: ~~~>>>======>>>======>>>~~~):
    1. ~~==  left reference (chromosome-cassette)
    2. ~~>>  left IS transposition (chromosome-IS)
    3. ==>>  left of mid IS (cassette-IS)
    4. ====  lost IS (cassette-cassette, no IS)
    5. >>==  right of mid IS (IS-cassette)
    6. >>~~  right IS transposition (IS-chromosome)
    7. ==~~  right reference (cassette-chromosome)

    Args:
        candidate: CandidateTnJc2 record
        chr_seq: Chromosome sequence
        tn_loc: TN location record for the chosen TN
        tn_seq: TN element sequence
        read_length: Read length for junction width

    Returns:
        Dict mapping JunctionType to sequence string
    """
    WID = read_length * 2  # junction width

    # Convert positions from 1-based (genomic coordinates) to 0-based (array indexing)
    # start and end are 1-based inclusive (from BLAST/junction positions)
    # Get segment scaffold to access normalized left/right positions
    seg_scaf = candidate.get_segment_scaffold(genome)
    pos_L = seg_scaf.left - 1  # Convert to 0-based start
    pos_R = seg_scaf.right - 1  # Convert to 0-based start

    # Outside positions (flanking the amplicon)
    pos_out_L = pos_L - 1
    pos_out_R = pos_R + 1

    # TN boundaries (convert 1-based inclusive to 0-based for array indexing)
    # loc_left and loc_right are 1-based inclusive (from BLAST/GenBank)
    tn_left = tn_loc.loc_left - 1  # Convert to 0-based start
    tn_right = tn_loc.loc_right  # 1-based inclusive -> 0-based exclusive end
    # TN orientation (from chosen TN)
    # Get orientation from candidate's tn_orientations list
    # MATLAB: IS_orientation = isjc2.chosen_IS{1}.orientation
    # Orientation computed in create_tnjc2: side_i * dir2(i), accounting for span_origin
    tn_orientation = 1  # default to forward
    if candidate.chosen_tn_id is not None:
        try:
            idx = candidate.tn_ids.index(candidate.chosen_tn_id)
            tn_orientation = candidate.tn_orientations[idx].value
        except (ValueError, IndexError):
            # Fallback if chosen_tn_id not found in tn_ids
            pass

    junctions = {}

    # Helper to get sequence with bounds checking
    def get_seq(seq: str, start: int, end: int, circular: bool = False) -> str:
        seq_len = len(seq)
        if circular and (start < 0 or end > seq_len):
            # Handle circular genome
            result = []
            for i in range(start, end):
                result.append(seq[i % seq_len])
            return ''.join(result)
        else:
            start = max(0, start)
            end = min(seq_len, end)
            return seq[start:end]

    # (1) Left reference: ~~==
    seq1 = get_seq(chr_seq, pos_out_L - WID + 1, pos_out_L + 1, circular=True)
    seq2 = get_seq(chr_seq, pos_L, pos_L + WID, circular=True)
    junctions[JunctionType.LEFT_REF] = seq1 + seq2

    # (2) Left IS transposition: ~~>>
    seq1 = get_seq(chr_seq, pos_out_L - WID + 1, pos_out_L + 1, circular=True)
    if tn_orientation == 1:
        seq2 = get_seq(tn_seq, tn_left, tn_left + WID)
    else:
        seq2 = reverse_complement(get_seq(tn_seq, tn_right - WID, tn_right))
    junctions[JunctionType.LEFT_IS_TRANS] = seq1 + seq2

    # (3) Left of mid IS: ==>>
    seq1 = get_seq(chr_seq, pos_R - WID + 1, pos_R + 1, circular=True)
    if tn_orientation == 1:
        seq2 = get_seq(tn_seq, tn_left, tn_left + WID)
    else:
        seq2 = reverse_complement(get_seq(tn_seq, tn_right - WID, tn_right))
    junctions[JunctionType.LEFT_MID_IS] = seq1 + seq2

    # (4) Lost IS: ====
    seq1 = get_seq(chr_seq, pos_R - WID + 1, pos_R + 1, circular=True)
    seq2 = get_seq(chr_seq, pos_L, pos_L + WID, circular=True)
    junctions[JunctionType.LOST_IS] = seq1 + seq2

    # (5) Right of mid IS: >>==
    if tn_orientation == 1:
        seq1 = get_seq(tn_seq, tn_right - WID, tn_right)
    else:
        seq1 = reverse_complement(get_seq(tn_seq, tn_left, tn_left + WID))
    seq2 = get_seq(chr_seq, pos_L, pos_L + WID, circular=True)
    junctions[JunctionType.RIGHT_MID_IS] = seq1 + seq2

    # (6) Right IS transposition: >>~~
    if tn_orientation == 1:
        seq1 = get_seq(tn_seq, tn_right - WID, tn_right)
    else:
        seq1 = reverse_complement(get_seq(tn_seq, tn_left, tn_left + WID))
    seq2 = get_seq(chr_seq, pos_out_R, pos_out_R + WID, circular=True)
    junctions[JunctionType.RIGHT_IS_TRANS] = seq1 + seq2

    # (7) Right reference: ==~~
    seq1 = get_seq(chr_seq, pos_R - WID + 1, pos_R + 1, circular=True)
    seq2 = get_seq(chr_seq, pos_out_R, pos_out_R + WID, circular=True)
    junctions[JunctionType.RIGHT_REF] = seq1 + seq2

    return junctions


def write_junctions_fasta(
    junctions: dict[JunctionType, str],
    output_path: Path,
) -> None:
    """Write junction sequences to FASTA file.

    Args:
        junctions: Dict mapping JunctionType to sequence
        output_path: Output FASTA file path
    """
    ensure_parent_dir(output_path)

    with open(output_path, 'w') as f:
        for jtype, seq in sorted(junctions.items(), key=lambda x: x[0].value):
            f.write(f">{jtype.value}\n{seq}\n")


class CreateSyntheticJunctionsStep(Step[RecordTypedDf[FilteredTnJc2]]):
    """Create synthetic junction FASTA files for each candidate.

    Creates 7 junction sequences per candidate for read alignment analysis.
    """

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        genome: Genome,
        tn_locs: RecordTypedDf[RefTnLoc],
        output_dir: Path,
        read_length: int = 150,
        force: Optional[bool] = None,
    ):
        self.filtered_tnjc2s = filtered_tnjc2s
        self.genome = genome
        self.tn_locs = tn_locs
        self.output_dir = Path(output_dir)
        self.read_length = read_length

        # Output files are per-candidate analysis directories
        self.analysis_dirs = [
            output_dir / filtered_tnjc2.analysis_dir
            for filtered_tnjc2 in filtered_tnjc2s
        ]

        super().__init__(
            input_files=[genome.fasta_path],
            output_files=[d / "junctions.fasta" for d in self.analysis_dirs],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[FilteredTnJc2]:
        """Create synthetic junctions for each candidate."""
        # Build TN lookup
        tn_lookup = {tn.tn_id: tn for tn in self.tn_locs}

        for filtered_tnjc2 in self.filtered_tnjc2s:
            analysis_dir = self.output_dir / filtered_tnjc2.analysis_dir

            # Get chosen TN
            chosen_tn_id = filtered_tnjc2.chosen_tn_id
            if chosen_tn_id is None or chosen_tn_id not in tn_lookup:
                # Skip candidates without valid chosen TN
                continue

            tn_loc = tn_lookup[chosen_tn_id]

            # Get scaffold sequence
            try:
                scaf = self.genome.get_seq_scaffold(tn_loc.tn_scaf)
            except KeyError:
                continue
            chr_seq = scaf.seq

            # Get TN sequence (loc_left and loc_right are 1-based inclusive, convert to 0-based for slicing)
            tn_seq = chr_seq[tn_loc.loc_left - 1:tn_loc.loc_right]

            # Create junctions
            junctions = create_synthetic_junctions(
                candidate=filtered_tnjc2,
                chr_seq=chr_seq,
                tn_loc=tn_loc,
                tn_seq=tn_seq,
                read_length=self.read_length,
                genome=self.genome,
            )

            # Write FASTA
            write_junctions_fasta(junctions, analysis_dir / "junctions.fasta")

        return self.filtered_tnjc2s

    def _save_output(self, output: RecordTypedDf[FilteredTnJc2]) -> None:
        """Output already saved in _calculate_output."""
        pass

    def load_outputs(self) -> RecordTypedDf[FilteredTnJc2]:
        """Return candidates (junction files are side effects)."""
        return self.filtered_tnjc2s

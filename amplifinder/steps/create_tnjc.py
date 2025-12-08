"""Step: Match junctions to TN elements (TnJc)."""

from pathlib import Path
from typing import Optional, List

from Bio.Seq import Seq

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDF, Junction, TnEndSeq, TnMatch, TnJunction, Orientation
from amplifinder.data_types.genome import Genome


class CreateTnJcStep(Step[RecordTypedDF[TnJunction]]):
    """Match junctions to TN elements by sequence comparison.

    For each junction, extracts flanking sequences and compares them to
    precomputed TN end sequences. Junctions matching a TN end are marked
    as TN-associated (TnJc).

    Based on assign_potential_ISs.m
    """

    def __init__(
        self,
        jc_df: RecordTypedDF[Junction],
        ref_tn_end_seqs: RecordTypedDF[TnEndSeq],
        genome: Genome,
        output_dir: Path,
        max_dist_to_tn: int,
        trim_jc_flanking: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            jc_df: All junctions (breseq + reference TN junctions combined)
            ref_tn_end_seqs: Precomputed TN end sequences (from CreateRefTnEndSeqsStep)
            genome: Reference genome (for junction sequence extraction)
            output_dir: Directory to write output
            max_dist_to_tn: Maximum distance from junction to TN boundary
            trim_jc_flanking: Trim junction edges to avoid misalignment
            force: Force re-run
        """
        self.jc_df = jc_df
        self.ref_tn_end_seqs = ref_tn_end_seqs
        self.genome = genome
        self.output_dir = Path(output_dir)
        self.max_dist_to_tn = max_dist_to_tn
        self.trim_jc_flanking = trim_jc_flanking

        self.output_file = self.output_dir / "TnJc.csv"

        super().__init__(
            input_files=[p for p in [genome.genbank_path, genome.fasta_path] if p],
            output_files=[self.output_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDF[TnJunction]:
        """Match junctions to TN elements."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.ref_seqs = self.genome.sequences
        tnjc_records = []

        for jc in self.jc_df:
            # Get junction flanking sequences
            seq1 = self._get_junction_seq(jc, side=1)
            seq2 = self._get_junction_seq(jc, side=2)

            # Match against TN sequences
            matches1 = self._find_tn_matches(seq1)
            matches2 = self._find_tn_matches(seq2)

            is_side1_tn = len(matches1) > 0
            is_side2_tn = len(matches2) > 0

            # Skip if neither side matches a TN
            if not is_side1_tn and not is_side2_tn:
                continue

            # Normalize: side 1 should be the TN side
            switched = not is_side1_tn and is_side2_tn
            if switched:
                matches = matches2
                jc = jc.switch_sides()
            else:
                matches = matches1

            tnjc_records.append(TnJunction.from_other(jc, matches=matches, switched=switched))

        tnjc = RecordTypedDF.from_records(tnjc_records, TnJunction)
        tnjc = tnjc.pipe(lambda df: df.sort_values(["scaf2", "pos2"]))

        info(f"Found {len(tnjc)} TN-associated junctions (TnJc)")
        return tnjc

    def _save_output(self, output: RecordTypedDF[TnJunction]) -> None:
        output.to_csv(self.output_file)

    def _get_junction_seq(self, jc: Junction, side: int) -> str:
        """Extract sequence at junction side."""
        scaf = jc.scaf1 if side == 1 else jc.scaf2
        pos = jc.pos1 if side == 1 else jc.pos2
        direction = jc.dir1 if side == 1 else jc.dir2
        flank_len = jc.flanking_left if side == 1 else jc.flanking_right

        if scaf not in self.ref_seqs:
            return ""

        ref_seq = self.ref_seqs[scaf]
        pos = pos - 1  # 0-based
        flank_len = flank_len - self.trim_jc_flanking - 1

        if flank_len <= 0:
            flank_len = 20

        if direction == Orientation.FORWARD:
            start = pos
            end = min(len(ref_seq), pos + flank_len)
            seq = ref_seq[start:end]
        else:
            start = max(0, pos - flank_len)
            end = pos + 1
            seq = ref_seq[start:end]
            seq = str(Seq(seq).reverse_complement())

        return seq

    def _find_tn_matches(self, jc_seq: str) -> List[TnMatch]:
        """Find TN elements matching a junction sequence."""
        if not jc_seq:
            return []

        matches = []
        threshold = self.max_dist_to_tn * 2

        for tn in self.ref_tn_end_seqs:
            # Check forward
            pos_fwd = tn.seq_fwd.find(jc_seq)
            if pos_fwd >= 0 and pos_fwd < threshold:
                dist = pos_fwd - self.max_dist_to_tn
                matches.append(TnMatch(tn.tn_id, tn.tn_side, dist))
                continue

            # Check reverse complement
            pos_rc = tn.seq_rc.find(jc_seq)
            if pos_rc >= 0 and pos_rc < threshold:
                dist = pos_rc - self.max_dist_to_tn
                matches.append(TnMatch(tn.tn_id, tn.tn_side, dist))

        return matches

    def load_outputs(self) -> RecordTypedDF[TnJunction]:
        """Load TnJc from output file."""
        return RecordTypedDF.from_csv(self.output_file, TnJunction)

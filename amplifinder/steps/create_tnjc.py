"""Step: Create TN-associated junction candidates (TNJC)."""

from dataclasses import replace
from pathlib import Path
from typing import Optional, List, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types.records import TNJC_SCHEMA, JUNCTION_SCHEMA, TnMatch, Tnjc, Junction


class CreateTNJCStep(Step[pd.DataFrame]):
    """Match junctions to TN elements by sequence comparison.

    For each junction, extracts flanking sequences and compares them to
    TN element end sequences. Junctions matching a TN end are marked
    as TN-associated (TNJC).

    Based on assign_potential_ISs.m
    """

    MAX_DIST_TO_TN = 20  # Max distance from junction to TN boundary
    TRIM_JC_FLANKING = 5  # Trim junction edges to avoid misalignment

    def __init__(
        self,
        jc_df: pd.DataFrame,
        tn_loc: pd.DataFrame,
        ref_fasta: Path,
        output_dir: Path,
        max_dist_to_tn: int = 20,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            jc_df: All junctions (breseq + reference TN junctions combined)
            tn_loc: TN locations DataFrame
            ref_fasta: Path to reference FASTA file
            output_dir: Directory to write output
            max_dist_to_tn: Maximum distance from junction to TN boundary
            force: Force re-run
        """
        self.jc_df = jc_df.copy()
        self.tn_loc = tn_loc
        self.ref_fasta = Path(ref_fasta)
        self.output_dir = Path(output_dir)
        self.max_dist_to_tn = max_dist_to_tn

        self.output_file = self.output_dir / "TNJC.csv"

        super().__init__(
            inputs=[self.ref_fasta],
            outputs=[self.output_file],
            force=force,
        )

    def _run(self) -> None:
        """Match junctions to TN elements."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Load reference sequences
        self.ref_seqs = self._load_reference_seqs()

        # Get TN end sequences (with margins)
        tn_seqs = self._get_tn_end_sequences()

        tnjc_records = []
        junctions = JUNCTION_SCHEMA.df_to_instances(self.jc_df)

        for jc in junctions:
            # Get junction flanking sequences
            seq1 = self._get_junction_seq(jc, side=1)
            seq2 = self._get_junction_seq(jc, side=2)

            # Match against TN sequences
            matches1 = self._find_tn_matches(seq1, tn_seqs)
            matches2 = self._find_tn_matches(seq2, tn_seqs)

            is_side1_tn = len(matches1) > 0
            is_side2_tn = len(matches2) > 0

            # Skip if neither side matches a TN
            if not is_side1_tn and not is_side2_tn:
                continue

            # Normalize: side 1 should be the TN side
            switched = not is_side1_tn and is_side2_tn
            if switched:
                matches = matches2
                jc = self._switch_jc_sides(jc)
            else:
                matches = matches1

            # Create TNJC record
            tnjc_records.append(Tnjc.from_instance(jc, matches=matches, switched=switched))

        if tnjc_records:
            tnjc = TNJC_SCHEMA.instances_to_df(tnjc_records)
            # Sort by chromosome position
            if "scaf2" in tnjc.columns and "pos2" in tnjc.columns:
                tnjc = tnjc.sort_values(["scaf2", "pos2"])
        else:
            tnjc = TNJC_SCHEMA.empty()

        TNJC_SCHEMA.to_csv(tnjc, self.output_file)
        info(f"Found {len(tnjc)} TN-associated junctions (TNJC)")

    def _load_reference_seqs(self) -> dict:
        """Load reference sequences from FASTA."""
        seqs = {}
        for record in SeqIO.parse(self.ref_fasta, "fasta"):
            seqs[record.id] = str(record.seq)
        return seqs

    def _get_tn_end_sequences(self) -> List[Tuple[int, str, str, str]]:
        """Get sequences at TN element ends with margins.

        Returns:
            List of (tn_id, tn_side, left_seq, right_seq) tuples
        """
        tn_seqs = []
        margin = self.max_dist_to_tn

        for idx, tn_row in self.tn_loc.iterrows():
            scaf = tn_row["TN_scaf"]
            if scaf not in self.ref_seqs:
                continue

            ref_seq = self.ref_seqs[scaf]
            left = tn_row["LocLeft"] - 1  # 0-based
            right = tn_row["LocRight"]  # 0-based exclusive

            # Left end sequence (including margin outside TN)
            left_start = max(0, left - margin)
            left_end = min(len(ref_seq), left + margin)
            left_seq = ref_seq[left_start:left_end]

            # Right end sequence (including margin outside TN)
            right_start = max(0, right - margin)
            right_end = min(len(ref_seq), right + margin)
            right_seq = ref_seq[right_start:right_end]

            # Also add reverse complements
            left_seq_rc = str(Seq(left_seq).reverse_complement())
            right_seq_rc = str(Seq(right_seq).reverse_complement())

            tn_seqs.append((idx, "left", left_seq, left_seq_rc))
            tn_seqs.append((idx, "right", right_seq, right_seq_rc))

        return tn_seqs

    def _get_junction_seq(self, jc: Junction, side: int) -> str:
        """Extract sequence at junction side.

        Args:
            jc: Junction record
            side: 1 or 2

        Returns:
            Sequence string at junction
        """
        scaf = jc.scaf1 if side == 1 else jc.scaf2
        pos = jc.pos1 if side == 1 else jc.pos2
        direction = jc.dir1 if side == 1 else jc.dir2
        flank_len = jc.flanking_left if side == 1 else jc.flanking_right

        if scaf not in self.ref_seqs:
            return ""

        ref_seq = self.ref_seqs[scaf]
        pos = pos - 1  # 0-based
        flank_len = flank_len - self.TRIM_JC_FLANKING - 1

        if flank_len <= 0:
            flank_len = 20

        if direction == 1:
            start = pos
            end = min(len(ref_seq), pos + flank_len)
            seq = ref_seq[start:end]
        else:
            start = max(0, pos - flank_len)
            end = pos + 1
            seq = ref_seq[start:end]
            seq = str(Seq(seq).reverse_complement())

        return seq

    def _find_tn_matches(
        self, jc_seq: str, tn_seqs: List[Tuple]
    ) -> List[TnMatch]:
        """Find TN elements matching a junction sequence.

        Args:
            jc_seq: Junction sequence
            tn_seqs: List of TN end sequences

        Returns:
            List of TnMatch(tn_id, side, distance) for matches
        """
        if not jc_seq:
            return []

        matches = []
        threshold = self.max_dist_to_tn * 2

        for tn_id, tn_side, tn_seq_fwd, tn_seq_rc in tn_seqs:
            # Check forward
            pos_fwd = tn_seq_fwd.find(jc_seq)
            if pos_fwd >= 0 and pos_fwd < threshold:
                dist = pos_fwd - self.max_dist_to_tn
                matches.append(TnMatch(tn_id, tn_side, dist))
                continue

            # Check reverse complement
            pos_rc = tn_seq_rc.find(jc_seq)
            if pos_rc >= 0 and pos_rc < threshold:
                dist = pos_rc - self.max_dist_to_tn
                matches.append(TnMatch(tn_id, tn_side, dist))

        return matches

    def _switch_jc_sides(self, jc: Junction) -> Junction:
        """Swap side 1 and side 2 fields in junction record."""
        return replace(
            jc,
            scaf1=jc.scaf2,
            scaf2=jc.scaf1,
            pos1=jc.pos2,
            pos2=jc.pos1,
            dir1=jc.dir2,
            dir2=jc.dir1,
            flanking_left=jc.flanking_right,
            flanking_right=jc.flanking_left,
        )

    def read_outputs(self) -> pd.DataFrame:
        """Load TNJC from output file."""
        return TNJC_SCHEMA.read_csv(self.output_file)

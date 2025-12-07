"""Step: Combine TN junction pairs (TNJC2)."""

from pathlib import Path
from typing import Optional, List, Tuple

from amplifinder.steps.base import Step
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDF, TnJunction, TnJunctionPair, TnMatch, Side, Orientation
from amplifinder.data_types.genome import Genome


class CreateTNJC2Step(Step[RecordTypedDF[TnJunctionPair]]):
    """Combine TN junctions into pairs (candidate amplicons).

    For each pair of junctions, checks:
    - (a) Same scaffold, opposite chromosome directions
    - (b) Same TN element, different sides (left/right)

    Based on MATLAB combine_ISJC_pairs.m and calculate_amplicon_length.m
    """

    def __init__(
        self,
        tnjc: RecordTypedDF[TnJunction],
        genome: Genome,
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tnjc: TN-associated junctions from CreateTNJCStep
            genome: Reference genome (for circularity and length)
            output_dir: Directory to write output
            force: Force re-run
        """
        self.tnjc = tnjc
        self.genome = genome
        self.output_dir = Path(output_dir)

        self.output_file = self.output_dir / "TNJC2.csv"

        super().__init__(
            inputs=[],
            outputs=[self.output_file],
            force=force,
        )

    def _run(self) -> None:
        """Combine junction pairs."""
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Convert to list for O(n^2) pairing
        junctions = list(self.tnjc)
        pairs = self._pair_junctions(junctions)

        if pairs:
            tnjc2 = RecordTypedDF.from_records(pairs, TnJunctionPair)
        else:
            tnjc2 = RecordTypedDF.empty(TnJunctionPair)

        tnjc2.to_csv(self.output_file)
        info(f"Found {len(tnjc2)} junction pairs (TNJC2)")

    def _pair_junctions(self, junctions: List[TnJunction]) -> List[TnJunctionPair]:
        """Find all valid junction pairs.

        Based on MATLAB combine_ISJC_pairs.m
        """
        pairs = []
        n = len(junctions)

        for i in range(n - 1):
            for j in range(i + 1, n):
                jc_i = junctions[i]
                jc_j = junctions[j]

                # (a1) Same scaffold for chromosome side (side 2)
                if jc_i.scaf2 != jc_j.scaf2:
                    continue

                # (a2) Opposing chromosome directions
                if jc_i.dir2 == jc_j.dir2:
                    continue

                # (b) Find matching TN: same ID, different sides
                matching = self._find_matching_tns(jc_i.matches, jc_j.matches)
                if not matching:
                    continue

                # Create pair record
                pair = self._create_pair(jc_i, jc_j, matching)
                pairs.append(pair)

        return pairs

    def _find_matching_tns(
        self,
        matches_i: List[TnMatch],
        matches_j: List[TnMatch],
    ) -> List[Tuple[int, Orientation]]:
        """Find TN elements that match both junctions on different sides.

        Returns list of (tn_id, orientation) tuples.
        """
        result = []

        for mi in matches_i:
            for mj in matches_j:
                # Same TN ID
                if mi.tn_id != mj.tn_id:
                    continue

                # Different sides (left vs right)
                if mi.side == mj.side:
                    continue

                # Orientation based on which side connects to which junction
                # If left side of TN connects to left junction -> same orientation (FORWARD)
                # If right side of TN connects to left junction -> inverted (REVERSE)
                orientation = Orientation.FORWARD if mi.side == Side.LEFT else Orientation.REVERSE

                result.append((mi.tn_id, orientation))

        return result

    def _create_pair(
        self,
        jc_i: TnJunction,
        jc_j: TnJunction,
        matching: List[Tuple[int, Orientation]],
    ) -> TnJunctionPair:
        """Create a junction pair record.

        Normalizes so that L (left) junction has lower chromosome position.
        """
        # Determine which is left/right based on position
        if jc_i.pos2 <= jc_j.pos2:
            jc_L, jc_R = jc_i, jc_j
        else:
            jc_L, jc_R = jc_j, jc_i

        # Extract TN IDs and compute combined orientation
        tn_ids = [m[0] for m in matching]
        orientations = [m[1] for m in matching]

        # Orientation: if all same -> that value, if mixed -> BOTH
        if all(o == Orientation.FORWARD for o in orientations):
            tn_orientation = Orientation.FORWARD
        elif all(o == Orientation.REVERSE for o in orientations):
            tn_orientation = Orientation.REVERSE
        else:
            tn_orientation = Orientation.BOTH

        # Span origin: if left junction points left (REVERSE), amplicon spans origin
        # (for circular genomes)
        span_origin = jc_L.dir2 == Orientation.REVERSE

        # Adjust orientation for span_origin (MATLAB: orientation * ISJC.dir2(i))
        if span_origin:
            tn_orientation = tn_orientation.opposite()

        # Calculate amplicon length
        amplicon_length, complementary_length = self._calculate_amplicon_length(
            jc_L.pos2, jc_R.pos2, jc_L.scaf2, span_origin
        )

        return TnJunctionPair(
            jc_num_L=jc_L.num,
            jc_num_R=jc_R.num,
            scaf_chr=jc_L.scaf2,
            pos_chr_L=jc_L.pos2,
            pos_chr_R=jc_R.pos2,
            pos_tn_L=jc_L.pos1,
            pos_tn_R=jc_R.pos1,
            dir_chr_L=jc_L.dir2,
            dir_chr_R=jc_R.dir2,
            dir_tn_L=jc_L.dir1,
            dir_tn_R=jc_R.dir1,
            tn_ids=tn_ids,
            tn_orientation=tn_orientation,
            span_origin=span_origin,
            amplicon_length=amplicon_length,
            complementary_length=complementary_length,
        )

    def _calculate_amplicon_length(
        self,
        pos_L: int,
        pos_R: int,
        scaf: str,
        span_origin: bool,
    ) -> Tuple[int, int]:
        """Calculate amplicon and complementary lengths.

        Based on MATLAB calculate_amplicon_length.m
        Handles circular genomes and origin-spanning amplicons.
        """
        # Basic length between positions
        raw_length = pos_R - pos_L + 1

        # Get genome properties
        is_circular = self.genome.circular
        scaf_length = self.genome.length

        if is_circular:
            if span_origin:
                # Amplicon spans origin: actual amplicon is the complement
                amplicon_length = scaf_length - raw_length
                complementary_length = raw_length
            else:
                # Normal case
                amplicon_length = raw_length
                complementary_length = scaf_length - raw_length
        else:
            # Linear chromosome
            if span_origin:
                # Can't span origin on linear - mark as infinite
                amplicon_length = float("inf")
                complementary_length = raw_length
            else:
                amplicon_length = raw_length
                complementary_length = float("inf")

        return int(amplicon_length) if amplicon_length != float("inf") else -1, \
               int(complementary_length) if complementary_length != float("inf") else -1

    def read_outputs(self) -> RecordTypedDF[TnJunctionPair]:
        """Load TNJC2 from output file."""
        return RecordTypedDF.from_csv(self.output_file, TnJunctionPair)


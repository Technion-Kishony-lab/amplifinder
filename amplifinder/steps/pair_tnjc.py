"""Step: Pair TN junctions into TnJc2."""

from pathlib import Path
from typing import Optional, List, Tuple

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.data_types import RecordTypedDf, TnJunction, RawTnJc2, RefTnSide, Side, Orientation
from amplifinder.data_types.genome import Genome


class PairTnJcToRawTnJc2Step(RecordTypedDfStep[RawTnJc2]):
    """Combine TN junctions into pairs (candidate amplicons).

    For each pair of junctions, checks:
    - (a) Same scaffold, opposite chromosome directions
    - (b) Same TN element, different sides (left/right)

    Based on MATLAB combine_ISJC_pairs.m and calculate_amplicon_length.m
    """

    def __init__(
        self,
        tnjcs: RecordTypedDf[TnJunction],
        genome: Genome,
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tnjcs: TN-associated junctions from CreateTnJcStep
            genome: Reference genome (for circularity and length per scaffold)
            output_dir: Directory to write output
            force: Force re-run
        """
        self.tnjcs = tnjcs
        self.genome = genome

        # Cache scaffold properties for multi-scaffold support

        super().__init__(
            output_dir=output_dir,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[RawTnJc2]:
        """Combine junction pairs."""

        # Convert to list for O(n^2) pairing
        junctions = list(self.tnjcs)
        pairs = self._pair_junctions(junctions)

        raw_tnjc2s = RecordTypedDf.from_records(pairs, RawTnJc2)

        self.log(f"Found {len(raw_tnjc2s)} junction pairs (RawTnJc2)")
        return raw_tnjc2s

    def _pair_junctions(self, junctions: List[TnJunction]) -> List[RawTnJc2]:
        """Find all valid junction pairs.

        Based on MATLAB combine_ISJC_pairs.m
        """
        pairs = []
        n = len(junctions)

        for i in range(n - 1):
            for j in range(i + 1, n):
                jc_i = junctions[i]
                jc_j = junctions[j]

                # (a1) Same scaffold for chromosome arm (arm 2)
                if jc_i.scaf2 != jc_j.scaf2:
                    continue

                # (a2) Opposing chromosome directions
                if jc_i.dir2 == jc_j.dir2:
                    continue

                # (b) Find matching TN: same ID, different sides
                matching_tns = self._find_matching_tns(jc_i.ref_tn_sides, jc_j.ref_tn_sides)
                if not matching_tns:
                    continue

                # Create pair record
                pair = self._create_pair(jc_i, jc_j, matching_tns)
                pairs.append(pair)

        return pairs

    def _find_matching_tns(
        self,
        matches_i: List[RefTnSide],
        matches_j: List[RefTnSide],
    ) -> List[Tuple[int, Side]]:
        """Find TN elements that match both junctions on different sides.

        Returns list of (tn_id, side_i) tuples where side_i is the TN side
        that junction i connects to.
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

                result.append((mi.tn_id, mi.side))

        return result

    def _create_pair(
        self,
        jc_i: TnJunction,
        jc_j: TnJunction,
        matching_tns: List[Tuple[int, Side]],
    ) -> RawTnJc2:
        """Create a junction pair record.

        Normalizes so that L (left) junction has lower chromosome position.
        """
        # Determine which is left/right based on position
        swapped = jc_i.pos2 > jc_j.pos2
        if swapped:
            jc_L, jc_R = jc_j, jc_i
        else:
            jc_L, jc_R = jc_i, jc_j

        # Extract TN IDs and compute orientations
        # MATLAB: orientation = side_i * dir2(i), where i is L
        tn_ids = [m[0] for m in matching_tns]
        # Use dir2 of the left junction (jc_L) for orientation calculation
        tn_orientations = [Orientation(side_i.value * jc_L.dir2.value * (-1 if swapped else 1))
                           for _, side_i in matching_tns]

        # Calculate amplicon length
        amplicon_length = self.genome.calc_length_between_scaf_points(jc_L.scaf2, jc_L.pos2, jc_R.pos2)

        return RawTnJc2(
            jc_num_L=jc_L.num,
            jc_num_R=jc_R.num,
            scaf=jc_L.scaf2,
            pos_scaf_L=jc_L.pos2,
            pos_scaf_R=jc_R.pos2,
            pos_tn_L=jc_L.pos1,
            pos_tn_R=jc_R.pos1,
            dir_tn_L=jc_L.dir1,
            dir_tn_R=jc_R.dir1,
            tn_ids=tn_ids,
            tn_orientations=tn_orientations,
            span_origin=jc_L.dir2 == Orientation.REVERSE,
            amplicon_length=amplicon_length,
        )

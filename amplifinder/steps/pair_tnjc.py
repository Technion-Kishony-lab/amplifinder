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
        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[RawTnJc2]:
        """Combine junction pairs."""

        # Convert to list for O(n^2) pairing
        junctions = list(self.tnjcs)
        pairs = self._pair_junctions(junctions)

        raw_tnjc2s = RecordTypedDf.from_records(pairs, RawTnJc2)

        return raw_tnjc2s

    def report_output_message(self, output: RecordTypedDf[RawTnJc2], *, from_cache: bool) -> Optional[str]:
        return f"Found {len(output)} junction pairs (RawTnJc2)"

    def _pair_junctions(self, junctions: List[TnJunction]) -> List[RawTnJc2]:
        """Find all valid junction pairs.

        Based on MATLAB combine_ISJC_pairs.m

        junctions: List[TnJunction] are all junctions that are matched to TN elements.
        NOTE: in all junctions, arm 2 is the chromosome arm
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
        i_ref_tn_sides: List[RefTnSide],
        j_ref_tn_sides: List[RefTnSide],
    ) -> List[Tuple[int, Side]]:
        """Find TN elements that match both junctions on different sides.

        Returns list of (tn_id, side_i) tuples where side_i is the TN side that junction i connects to.
        """
        result = []

        for i_tn_side in i_ref_tn_sides:
            for j_tn_side in j_ref_tn_sides:
                # Same TN ID
                if i_tn_side.tn_id != j_tn_side.tn_id:
                    continue

                # Different sides (left vs right)
                if i_tn_side.side == j_tn_side.side:
                    continue

                result.append((i_tn_side.tn_id, i_tn_side.side))

        return result

    def _create_pair(
        self,
        jc_i: TnJunction,
        jc_j: TnJunction,
        matching_tns: List[Tuple[int, Side]],
    ) -> RawTnJc2:
        """Create a junction pair record.

        Assigns start/end based on forward strand direction:
        - Start (S): junction where forward strand starts (dir2 == FORWARD)
        - End (E): junction where forward strand ends (dir2 == REVERSE)
        """
        # Determine start/end based on forward strand direction
        if jc_i.dir2 == Orientation.FORWARD:
            jc_S = jc_i
            jc_E = jc_j
            swapped = False
        else:
            jc_S = jc_j
            jc_E = jc_i
            swapped = True

        # Extract TN IDs and compute orientations
        tn_ids = [m[0] for m in matching_tns]
        tn_orientations = [Orientation(side_i.value * jc_i.dir2.value * (-1 if swapped else 1))
                           for _, side_i in matching_tns]

        # Create RawTnJc2 with placeholder amplicon_length, then compute it
        pair = RawTnJc2(
            jc_num_S=jc_S.num,
            jc_num_E=jc_E.num,
            scaf=jc_S.scaf2,
            start=jc_S.pos2,
            end=jc_E.pos2,
            pos_tn_S=jc_S.pos1,
            pos_tn_E=jc_E.pos1,
            dir_tn_S=jc_S.dir1,
            dir_tn_E=jc_E.dir1,
            tn_ids=tn_ids,
            tn_orientations=tn_orientations,
        )
        pair.compute_and_store_amplicon_length(self.genome)
        return pair

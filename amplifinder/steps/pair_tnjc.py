"""Step: Pair TN junctions into TnJc2."""

from pathlib import Path
from typing import Optional, List

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.data_types import RecordTypedDf, RefTnJunction, TnJunction, RawTnJc2, Orientation
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

    def _pair_junctions(self, junctions: List[TnJunction]) -> List[RawTnJc2]:
        """Find all valid junction pairs.

        Based on MATLAB combine_ISJC_pairs.m

        junctions: List[TnJunction] are all junctions that are matched to TN elements.
        note: in all junctions, arm 2 is the chromosome arm
        """

        pairs = []
        n = len(junctions)

        for i in range(n - 1):
            for j in range(i + 1, n):
                jc_i = junctions[i]
                jc_j = junctions[j]

                # (a) Same scaffold for chromosome arm (arm 2)
                if jc_i.scaf2 != jc_j.scaf2:
                    continue

                # (a) Opposing chromosome directions
                if jc_i.dir2 == jc_j.dir2:
                    continue

                # (c) Find matching TN: same ID, different sides
                matching_tn_sides = RawTnJc2.find_matching_tn_sides(jc_i.ref_tn_sides, jc_j.ref_tn_sides)
                if not matching_tn_sides:
                    continue

                # (d) A reference TN junction must appear in the matching TNs
                matching_tn_ids = [tn_side_i.tn_id for tn_side_i, tn_side_j in matching_tn_sides]
                if jc_i.is_ref_tn_junction() and jc_i.ref_tn_side.tn_id not in matching_tn_ids:
                    continue
                if jc_j.is_ref_tn_junction() and jc_j.ref_tn_side.tn_id not in matching_tn_ids:
                    continue

                # Create pair record
                pair = self._create_pair(jc_i, jc_j)
                pairs.append(pair)

        return pairs

    def _create_pair(
        self,
        jc_i: TnJunction,
        jc_j: TnJunction,
    ) -> RawTnJc2:
        """Create a junction pair record.

        Assigns start/end based on forward strand direction:
        - Start (S): junction where forward strand starts (dir2 == FORWARD)
        - End (E): junction where forward strand ends (dir2 == REVERSE)
        """
        # Determine start/end based on forward strand direction
        if jc_i.dir2 == Orientation.FORWARD:
            tnjc_left, tnjc_right = jc_i, jc_j
        else:
            tnjc_left, tnjc_right = jc_j, jc_i

        # Get scaffold object
        scaffold = self.genome.get_scaffold(jc_S.scaf2)

        # Create RawTnJc2 with the two junction objects and scaffold
        pair = RawTnJc2(tnjc_left=tnjc_left, tnjc_right=tnjc_right, scaffold=scaffold)
        return pair

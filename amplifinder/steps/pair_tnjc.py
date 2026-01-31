"""Step: Pair TN junctions into TnJc2."""

from pathlib import Path
from typing import Optional, List

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.data_types.genome import Genome
from amplifinder.data_types.scaffold import SeqScaffold
from amplifinder.data_types import RecordTypedDf, TnJunction, RawTnJc2, Orientation, BaseEvent


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
        transposition_threshold: int,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tnjcs: TN-associated junctions from CreateTnJcStep
            genome: Reference genome (for circularity and length per scaffold)
            output_dir: Directory to write output
            transposition_threshold: Max distance for transposition classification
            force: Force re-run
        """
        self.tnjcs = tnjcs
        self.genome = genome
        self.transposition_threshold = transposition_threshold
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
        pair_num = 1

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
                matching_tn_ids = [tn_side_i.tn_id for tn_side_i, tn_side_j in matching_tn_sides]

                # (d) A reference TN junction must appear in the matching TNs
                if jc_i.is_ref_tn_junction() and jc_i.ref_tn_side.tn_id not in matching_tn_ids:
                    continue
                if jc_j.is_ref_tn_junction() and jc_j.ref_tn_side.tn_id not in matching_tn_ids:
                    continue

                # Create pair record
                pair = self._create_pair(jc_i, jc_j, pair_id=pair_num)
                pairs.append(pair)
                pair_num += 1

        return pairs

    def _create_pair(
        self,
        jc_i: TnJunction,
        jc_j: TnJunction,
        pair_id: int,
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
        scaffold = self.genome.get_scaffold(tnjc_left.scaf2)

        # Compute base_event before creating RawTnJc2
        base_event = self._compute_base_event(tnjc_left, tnjc_right)
        
        # Create RawTnJc2 with the two junction objects, scaffold, and base_event
        pair = RawTnJc2(
            pair_id=pair_id, 
            tnjc_left=tnjc_left, 
            tnjc_right=tnjc_right, 
            scaffold=scaffold,
            base_event=base_event,
        )
        return pair
    
    def _compute_base_event(
        self, tnjc_left: TnJunction, tnjc_right: TnJunction
    ) -> BaseEvent:
        """Classify base event type for junction pair."""
        # Check for reference TN
        if tnjc_left.is_ref_tn_junction() and \
                tnjc_right.is_ref_tn_junction() and \
                tnjc_left.ref_tn_side.tn_id == tnjc_right.ref_tn_side.tn_id:
            return BaseEvent.REFERENCE_TN
        
        # Check for transposition
        left_pos = tnjc_left.pos2
        right_pos = tnjc_right.pos2
        if abs(left_pos - right_pos) < self.transposition_threshold:
            return BaseEvent.TRANSPOSITION
        
        return BaseEvent.LOCUS_JOINING

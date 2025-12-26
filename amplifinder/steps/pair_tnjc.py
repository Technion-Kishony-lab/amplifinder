"""Step: Pair TN junctions into TnJc2."""

from pathlib import Path
from typing import Optional, List, Tuple

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.logger import info
from amplifinder.data_types import RecordTypedDf, TnJunction, RawTnJc2, RefTnSide, Side, Orientation
from amplifinder.data_types.genome import Genome
from amplifinder.utils.tools import ensure_dir


class PairTnJcToRawTnJc2Step(RecordTypedDfStep[RawTnJc2]):
    """Combine TN junctions into pairs (candidate amplicons).

    For each pair of junctions, checks:
    - (a) Same scaffold, opposite chromosome directions
    - (b) Same TN element, different sides (left/right)

    Based on MATLAB combine_ISJC_pairs.m and calculate_amplicon_length.m
    """

    def __init__(
        self,
        tnjc: RecordTypedDf[TnJunction],
        genome: Genome,
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            tnjc: TN-associated junctions from CreateTnJcStep
            genome: Reference genome (for circularity and length per scaffold)
            output_dir: Directory to write output
            force: Force re-run
        """
        self.tnjc = tnjc
        self.genome = genome
        
        # Cache scaffold properties for multi-scaffold support
        self._scaf_lengths = {rec.name: len(rec.seq) for rec in genome.records}
        self._scaf_circular = {
            rec.name: rec.annotations.get("topology", "linear").lower() == "circular"
            for rec in genome.records
        }
        
        super().__init__(
            output_dir=output_dir,
            input_files=[],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[RawTnJc2]:
        """Combine junction pairs."""
        ensure_dir(self.output_dir)

        # Convert to list for O(n^2) pairing
        junctions = list(self.tnjc)
        pairs = self._pair_junctions(junctions)

        tnjc2 = RecordTypedDf.from_records(pairs, RawTnJc2)

        info(f"Found {len(tnjc2)} junction pairs (RawTnJc2)")
        return tnjc2

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

                # (a1) Same scaffold for chromosome side (side 2)
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
        switched = jc_i.pos2 > jc_j.pos2
        if switched:
            jc_L, jc_R = jc_j, jc_i
        else:
            jc_L, jc_R = jc_i, jc_j

        # Extract TN IDs and compute orientations
        # MATLAB: orientation = side_i * dir2(i), where i is L
        tn_ids = [m[0] for m in matching_tns]
        # Use dir2 of the left junction (jc_L) for orientation calculation
        tn_orientations = [Orientation(side_i.value * jc_L.dir2.value * (-1 if switched else 1)) for _, side_i in matching_tns]

        # Span origin: if left junction points left (REVERSE), amplicon spans origin
        span_origin = jc_L.dir2 == Orientation.REVERSE

        # Calculate amplicon length
        amplicon_length, complementary_length = self._calculate_amplicon_length(
            jc_L.pos2, jc_R.pos2, jc_L.scaf2, span_origin
        )

        return RawTnJc2(
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
            tn_orientations=tn_orientations,
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

        # Get scaffold-specific properties (supports multi-scaffold genomes)
        scaf_length = self._scaf_lengths[scaf]
        is_circular = self._scaf_circular.get(scaf, False)

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

        amp = int(amplicon_length) if amplicon_length != float("inf") else -1
        comp = int(complementary_length) if complementary_length != float("inf") else -1
        return amp, comp

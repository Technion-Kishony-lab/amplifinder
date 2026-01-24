"""Step 8: Classify junction pair structures based on TN relationships."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import BaseRawEvent, RawEvent
from amplifinder.data_types.genome import Genome
from amplifinder.data_types import RecordTypedDf, CoveredTnJc2, RefTn, Side, SingleLocusLinkedTnJc2, TnJunction

from .base import RecordTypedDfStep


class LinkTnJc2ToSingleLocusPairsStep(RecordTypedDfStep[SingleLocusLinkedTnJc2]):
    """Classify junction pair structures based on TN relationships.

    This step analyzes how junction pairs relate to reference TN elements
    and other single-locus pairs to classify them as:
    - reference: TN in its reference location
    - transposition: short de novo insertion
    - flanked: both sides share IS with single-locus pairs
    - hemi-flanked: one side shares IS
    - unflanked: neither side shares IS
    """

    def __init__(
        self,
        covered_tnjc2s: RecordTypedDf[CoveredTnJc2],
        genome: Genome,
        tn_locs: RecordTypedDf[RefTn],
        output_dir: Path,
        transposition_threshold: int = 30,
        force: Optional[bool] = None,
    ):
        self.covered_tnjc2s = covered_tnjc2s
        self.genome = genome
        self.tn_locs = tn_locs
        self.transposition_threshold = transposition_threshold

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[SingleLocusLinkedTnJc2]:
        """Classify junction pairs."""
        tnjc2s = self.covered_tnjc2s.to_records()
        base_raw_events = [self._compute_base_raw_event(tnjc2) for tnjc2 in tnjc2s]

        # First pass: create all SingleLocusLinkedTnJc2 objects with empty matchings
        linked_tnjc2s: list[SingleLocusLinkedTnJc2] = []
        for i, tncj2_i in enumerate(tnjc2s):
            linked_tnjc2s.append(SingleLocusLinkedTnJc2.from_other(
                tncj2_i,
                single_locus_tnjc2_left_matchings=[],
                single_locus_tnjc2_right_matchings=[],
                base_raw_event=base_raw_events[i],
            ))

        # Second pass: populate the matchings with SingleLocusLinkedTnJc2 objects
        for i, linked_i in enumerate(linked_tnjc2s):
            left_matches = self._find_all_single_locus_matching_tnjc2s(
                linked_i.tnjc_left, linked_tnjc2s, base_raw_events, exclude_idx=i)
            right_matches = self._find_all_single_locus_matching_tnjc2s(
                linked_i.tnjc_right, linked_tnjc2s, base_raw_events, exclude_idx=i)
            linked_i.single_locus_tnjc2_left_matchings = left_matches
            linked_i.single_locus_tnjc2_right_matchings = right_matches

        return RecordTypedDf.from_records(linked_tnjc2s, SingleLocusLinkedTnJc2)

    def _compute_base_raw_event(self, tnjc2: CoveredTnJc2) -> BaseRawEvent:
        if tnjc2.tnjc_left.is_ref_tn_junction() and \
                tnjc2.tnjc_right.is_ref_tn_junction() and \
                tnjc2.tnjc_left.ref_tn_side.tn_id == tnjc2.tnjc_right.ref_tn_side.tn_id and \
                not tnjc2.span_origin:
            return BaseRawEvent.REFERENCE
        elif abs(tnjc2.left - tnjc2.right) < self.transposition_threshold:
            return BaseRawEvent.TRANSPOSITION
        else:
            return BaseRawEvent.LOCUS_JOINING

    def _find_all_single_locus_matching_tnjc2s(
        self, tnjc_i: TnJunction, tnjc2s: list[SingleLocusLinkedTnJc2],
        base_raw_events: list[BaseRawEvent], exclude_idx: int
    ) -> list[tuple[SingleLocusLinkedTnJc2, Side]]:
        matches = []
        for j, tnjc2_j in enumerate(tnjc2s):
            if j == exclude_idx:  # Don't match with self
                continue
            if not base_raw_events[j].is_single_locus():
                continue
            if tnjc_i == tnjc2_j.tnjc_left:
                matches.append((tnjc2_j, Side.LEFT))
            elif tnjc_i == tnjc2_j.tnjc_right:
                matches.append((tnjc2_j, Side.RIGHT))
        return matches

    def report_output_message(self, output: RecordTypedDf[SingleLocusLinkedTnJc2]) -> Optional[str]:
        records: list[SingleLocusLinkedTnJc2] = output.to_records()
        counts: dict[RawEvent, int] = {name: 0 for name in RawEvent}
        for record in records:
            counts[record.raw_event] += 1
        lines = [f"{event.value:22s}: {count:4d}" for event, count in counts.items()]
        return "\n" + "\n".join(lines)

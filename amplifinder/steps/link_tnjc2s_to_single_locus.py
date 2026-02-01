"""Step 8: Classify junction pair structures based on TN relationships."""

from pathlib import Path
from typing import Optional

from amplifinder.logger import logger
from amplifinder.data_types import (
    RecordTypedDf, RawTnJc2, RefTn, Side, SingleLocusLinkedTnJc2, TnJunction, Architecture,
)
from amplifinder.data_types.tnjc2s import TnJc2AndSide

from .base import RecordTypedDfStep, count_classifications, format_classification_counts

REPORT_CATEGORIES = (
    Architecture.REFERENCE_TN,
    Architecture.TRANSPOSITION,
    Architecture.FLANKED,
    Architecture.HEMI_FLANKED_LEFT,
    Architecture.HEMI_FLANKED_RIGHT,
    Architecture.UNFLANKED,
)


def _issue_debug_message_if_multi_single_locus_match(tnjc2, matches):
    if len(matches) > 1:
        logger.debug_message(f"tnjc: {tnjc2}\nmacthes: {matches}", "tnjc2_matches_multiple_single_locus", 1)


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
        covered_tnjc2s: RecordTypedDf[RawTnJc2],
        tn_locs: RecordTypedDf[RefTn],
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        self.covered_tnjc2s = covered_tnjc2s
        self.tn_locs = tn_locs

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[SingleLocusLinkedTnJc2]:
        """Classify junction pairs."""
        tnjc2s = self.covered_tnjc2s.to_records()

        # First pass: convert all to SingleLocusLinkedTnJc2 with empty matchings
        linked_tnjc2s: list[SingleLocusLinkedTnJc2] = []
        for tncj2_i in tnjc2s:
            linked = SingleLocusLinkedTnJc2.from_other(
                tncj2_i,
                single_locus_tnjc2_left_matchings=[],
                single_locus_tnjc2_right_matchings=[],
                base_event=tncj2_i.base_event,
            )
            linked_tnjc2s.append(linked)

        # Second pass: find matchings and update
        final_linked_tnjc2s: list[SingleLocusLinkedTnJc2] = []
        for i, linked_i in enumerate(linked_tnjc2s):
            left_matches = self._find_all_single_locus_matching_tnjc2s(
                linked_i.tnjc_left, linked_tnjc2s, exclude_idx=i)
            right_matches = self._find_all_single_locus_matching_tnjc2s(
                linked_i.tnjc_right, linked_tnjc2s, exclude_idx=i)

            final_linked = SingleLocusLinkedTnJc2.from_other(
                linked_i,
                single_locus_tnjc2_left_matchings=left_matches,
                single_locus_tnjc2_right_matchings=right_matches,
                base_event=linked_i.base_event,
            )
            _issue_debug_message_if_multi_single_locus_match(final_linked, left_matches)
            _issue_debug_message_if_multi_single_locus_match(final_linked, right_matches)
            final_linked_tnjc2s.append(final_linked)

        return RecordTypedDf.from_records(final_linked_tnjc2s, SingleLocusLinkedTnJc2)

    def report_output_message(self, output: RecordTypedDf[SingleLocusLinkedTnJc2]) -> Optional[str]:
        counts = count_classifications(
            [record.raw_event for record in output],
            REPORT_CATEGORIES,
        )
        return format_classification_counts(
            REPORT_CATEGORIES,
            counts,
            label_fn=lambda event: event.description,
        )

    def _find_all_single_locus_matching_tnjc2s(
        self, tnjc_i: TnJunction, tnjc2s: list[SingleLocusLinkedTnJc2], exclude_idx: int
    ) -> list[TnJc2AndSide]:
        matches = []
        for j, tnjc2_j in enumerate(tnjc2s):
            if j == exclude_idx:  # Don't match with self
                continue
            if not tnjc2_j.base_event.is_single_locus():
                continue
            if tnjc_i == tnjc2_j.tnjc_left:
                matches.append(TnJc2AndSide(tnjc2_j, Side.LEFT))
            elif tnjc_i == tnjc2_j.tnjc_right:
                matches.append(TnJc2AndSide(tnjc2_j, Side.RIGHT))
        return matches

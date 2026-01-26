"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, SingleLocusLinkedTnJc2,
)
from amplifinder.steps.base import RecordTypedDfStep


class FilterTnJc2CandidatesStep(RecordTypedDfStep[SingleLocusLinkedTnJc2]):
    """Filter candidates by amplicon length and copy number.

    Filters out candidates that are:
    - Too short (< min_amplicon_length)
    - Too long (> max_amplicon_length)
    - Copy number too low (< copy_number_threshold)

    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        linked_tnjc2s: RecordTypedDf[SingleLocusLinkedTnJc2],
        output_dir: Path,
        min_amplicon_length: int,
        max_amplicon_length: int,
        copy_number_threshold: float,
        force: Optional[bool] = None,
    ):
        self.linked_tnjc2s = linked_tnjc2s
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.copy_number_threshold = copy_number_threshold

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[SingleLocusLinkedTnJc2]:
        """Filter candidates by amplicon length and copy number."""
        filtered_tnjc2s: list[SingleLocusLinkedTnJc2] = []

        for tnjc2 in self.linked_tnjc2s:
            # Filter by length
            if not (self.min_amplicon_length <= tnjc2.amplicon_length <= self.max_amplicon_length):
                continue

            # Filter out candidates without a chosen TN
            if tnjc2.chosen_tn_id is None:
                continue

            # Filter by copy number threshold
            if tnjc2.copy_number < self.copy_number_threshold:
                continue

            filtered_tnjc2s.append(tnjc2)

        return RecordTypedDf.from_records(filtered_tnjc2s, SingleLocusLinkedTnJc2)

    def report_output_message(self, output: RecordTypedDf[SingleLocusLinkedTnJc2]) -> Optional[str]:
        return f"Filtered to {len(output)} candidates (from {len(self.linked_tnjc2s)} classified)"

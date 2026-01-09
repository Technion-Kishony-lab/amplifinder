"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, SingleLocusLinkedTnJc2,
)
from amplifinder.steps.base import RecordTypedDfStep


class FilterTnJc2CandidatesStep(RecordTypedDfStep[SingleLocusLinkedTnJc2]):
    """Filter candidates by amplicon length.

    Filters out candidates that are:
    - Too short (< min_amplicon_length)
    - Too long (> max_amplicon_length)

    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        classified_tnjc2s: RecordTypedDf[SingleLocusLinkedTnJc2],
        output_dir: Path,
        min_amplicon_length: int = 30,
        max_amplicon_length: int = 1_000_000,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[SingleLocusLinkedTnJc2]:
        """Filter candidates by amplicon length."""
        filtered_tnjc2s: list[SingleLocusLinkedTnJc2] = []

        for tnjc2 in self.classified_tnjc2s:
            # Filter by length
            if not (self.min_amplicon_length <= tnjc2.amplicon_length <= self.max_amplicon_length):
                continue

            # Filter out candidates without a chosen TN
            if tnjc2.chosen_tn_id is None:
                continue

            filtered_tnjc2s.append(tnjc2)

        return RecordTypedDf.from_records(filtered_tnjc2s, SingleLocusLinkedTnJc2)

    def report_output_message(self, output: RecordTypedDf[SingleLocusLinkedTnJc2], *, from_cache: bool) -> Optional[str]:
        return f"Filtered to {len(output)} candidates (from {len(self.classified_tnjc2s)} classified)"

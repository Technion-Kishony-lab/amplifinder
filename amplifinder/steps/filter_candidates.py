"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, SingleLocusLinkedTnJc2, Architecture,
)
from amplifinder.steps.base import RecordTypedDfStep


def _count_raw_events(records: list[SingleLocusLinkedTnJc2]) -> dict[Architecture, int]:
    """Count occurrences of each RawEvent type in a list of records."""
    counts: dict[Architecture, int] = {name: 0 for name in Architecture}
    for record in records:
        counts[record.raw_event] += 1
    return counts


class FilterTnJc2CandidatesStep(RecordTypedDfStep[SingleLocusLinkedTnJc2]):
    """Filter candidates by amplicon length and copy number.

    Filters out candidates that are:
    - Too short (< min_amplicon_length)
    - Too long (> max_amplicon_length)
    - Copy number out of range (< del_copy_number_threshold and >= copy_number_threshold)

    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        linked_tnjc2s: RecordTypedDf[SingleLocusLinkedTnJc2],
        output_dir: Path,
        min_amplicon_length: int,
        max_amplicon_length: int,
        replication_copy_number_threshold: float,
        delition_copy_number_threshold: float,
        force: Optional[bool] = None,
    ):
        self.linked_tnjc2s = linked_tnjc2s
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.replication_copy_number_threshold = replication_copy_number_threshold
        self.delition_copy_number_threshold = delition_copy_number_threshold

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

            # Filter by copy number: keep amplifications and deletions
            if not (tnjc2.copy_number >= self.replication_copy_number_threshold or 
                    tnjc2.copy_number < self.delition_copy_number_threshold):
                continue

            filtered_tnjc2s.append(tnjc2)

        return RecordTypedDf.from_records(filtered_tnjc2s, SingleLocusLinkedTnJc2)

    def report_output_message(self, output: RecordTypedDf[SingleLocusLinkedTnJc2]) -> Optional[str]:
        before_counts = _count_raw_events(self.linked_tnjc2s.to_records())
        after_counts = _count_raw_events(output.to_records())
        
        lines = [f"{event.value:22s}: {before_counts[event]:4d} -> {after_counts[event]:4d}" 
                 for event in Architecture]
        return "\n" + "\n".join(lines)

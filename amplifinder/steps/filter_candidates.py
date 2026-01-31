"""Step 9: Filter candidates by copy number."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, Architecture,
)
from amplifinder.steps.base import RecordTypedDfStep


def _count_raw_events(records: list[CoveredTnJc2]) -> dict[Architecture, int]:
    """Count occurrences of each RawEvent type in a list of records."""
    counts: dict[Architecture, int] = {name: 0 for name in Architecture}
    for record in records:
        counts[record.raw_event] += 1
    return counts


class FilterTnJc2CandidatesStep(RecordTypedDfStep[CoveredTnJc2]):
    """Filter candidates by copy number.

    Filters out candidates that are:
    - Copy number out of range (< del_copy_number_threshold and >= copy_number_threshold)
    - Missing a chosen TN

    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        linked_tnjc2s: RecordTypedDf[CoveredTnJc2],
        output_dir: Path,
        replication_copy_number_threshold: float,
        deletion_copy_number_threshold: float,
        force: Optional[bool] = None,
    ):
        self.linked_tnjc2s = linked_tnjc2s
        self.replication_copy_number_threshold = replication_copy_number_threshold
        self.deletion_copy_number_threshold = deletion_copy_number_threshold

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[CoveredTnJc2]:
        """Filter candidates by copy number and chosen TN."""
        filtered_tnjc2s: list[CoveredTnJc2] = []

        for tnjc2 in self.linked_tnjc2s:
            # Filter out candidates without a chosen TN
            if tnjc2.chosen_tn_id is None:
                continue

            # Filter by copy number: keep amplifications and deletions
            if self.deletion_copy_number_threshold <= tnjc2.copy_number < self.replication_copy_number_threshold:
                continue

            filtered_tnjc2s.append(tnjc2)

        return RecordTypedDf.from_records(filtered_tnjc2s, CoveredTnJc2)

    def report_output_message(self, output: RecordTypedDf[CoveredTnJc2]) -> Optional[str]:
        before_counts = _count_raw_events(self.linked_tnjc2s.to_records())
        after_counts = _count_raw_events(output.to_records())

        lines = [(f"{event.description:30s}: {before_counts[event]:4d} -> "
                  f"{after_counts[event]:4d}")
                 for event in Architecture]
        return "\n" + "\n".join(lines)

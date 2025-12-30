"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, FilteredTnJc2,
)
from amplifinder.steps.base import RecordTypedDfStep


class FilterTnJc2CandidatesStep(RecordTypedDfStep[FilteredTnJc2]):
    """Filter candidates by amplicon length.

    Filters out candidates that are:
    - Too short (< min_amplicon_length)
    - Too long (> max_amplicon_length)

    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        output_dir: Path,
        min_amplicon_length: int = 30,
        max_amplicon_length: int = 1_000_000,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length

        super().__init__(
            output_dir=output_dir,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[FilteredTnJc2]:
        """Filter candidates by amplicon length."""
        candidate_records = []

        for i, classified in enumerate(self.classified_tnjc2s):
            # Filter by length
            if (classified.amplicon_length <= self.min_amplicon_length or
                    classified.amplicon_length >= self.max_amplicon_length):
                continue

            # Generate analysis directory name
            # Format: jc_{start}_{end}_{tn_id:03d}_L{read_len}
            # Note: read_len defaults to 150 if not provided in config
            start = classified.start
            end = classified.end
            tn_id = classified.chosen_tn_id or (i + 1)
            # Default read length (can be overridden by config in pipeline)
            read_len = 150

            analysis_dir = f"jc_{start}_{end}_{tn_id:03d}_L{read_len}"

            # Create candidate record
            candidate = FilteredTnJc2.from_other(
                classified,
                analysis_dir=analysis_dir,
            )
            candidate_records.append(candidate)

        result = RecordTypedDf.from_records(candidate_records, FilteredTnJc2)

        return result

    def report_output_message(self, output: RecordTypedDf[FilteredTnJc2], *, from_cache: bool) -> Optional[str]:
        return f"Filtered to {len(output)} candidates (from {len(self.classified_tnjc2s)} classified)"

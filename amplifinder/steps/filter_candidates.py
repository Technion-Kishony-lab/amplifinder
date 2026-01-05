"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, FilteredTnJc2, Side,
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
        read_length: int = 150,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.read_length = read_length

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[FilteredTnJc2]:
        """Filter candidates by amplicon length."""
        filtered_tnjc2s = []

        for i, tnjc2 in enumerate(self.classified_tnjc2s):
            # Filter by length
            if not (self.min_amplicon_length <= tnjc2.amplicon_length <= self.max_amplicon_length):
                continue
            tn_id = tnjc2.chosen_tn_id
            if tn_id is None or side_S is None:
                continue

            # Generate analysis directory name
            # Format: jc_{start}_{end}_{tn_id:03d}_{side}_L{read_len}
            side_S = tnjc2.chosen_tn_side_left
            assert side_S is not None
            side_str = "S" if side_S == Side.START else "E"
            analysis_dir = f"jc_{tnjc2.left}-{tnjc2.right}_{tn_id:03d}{side_str}_{self.read_length}bp"

            # Create candidate record
            filtered_tnjc2s.append(FilteredTnJc2.from_other(
                tnjc2,
                analysis_dir=analysis_dir,
            ))

        return RecordTypedDf.from_records(filtered_tnjc2s, FilteredTnJc2)

    def report_output_message(self, output: RecordTypedDf[FilteredTnJc2], *, from_cache: bool) -> Optional[str]:
        return f"Filtered to {len(output)} candidates (from {len(self.classified_tnjc2s)} classified)"

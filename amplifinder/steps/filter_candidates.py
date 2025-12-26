"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, CandidateTnJc2,
)
from amplifinder.steps.base import Step
from amplifinder.logger import info


class FilterCandidatesStep(Step[RecordTypedDf[CandidateTnJc2]]):
    """Filter candidates by amplicon length.
    
    Filters out candidates that are:
    - Too short (< min_amplicon_length)
    - Too long (> max_amplicon_length)
    
    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        classified_tnjc2: RecordTypedDf[ClassifiedTnJc2],
        output_dir: Path,
        min_amplicon_length: int = 30,
        max_amplicon_length: int = 1_000_000,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2 = classified_tnjc2
        self.output_dir = Path(output_dir)
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        
        self.output_file = output_dir / "tn_jc2_candidates.csv"
        
        super().__init__(
            input_files=[],
            output_files=[self.output_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[CandidateTnJc2]:
        """Filter candidates by amplicon length."""
        candidate_records = []
        
        for i, classified in enumerate(self.classified_tnjc2):
            # Filter by length
            if (classified.amplicon_length <= self.min_amplicon_length or
                classified.amplicon_length >= self.max_amplicon_length):
                continue
            
            # Generate analysis directory name
            # Format: jc_{start}_{end}_{tn_id:03d}_L{read_len}
            # Note: read_len defaults to 150 if not provided in config
            start = classified.pos_chr_L
            end = classified.pos_chr_R
            tn_id = classified.chosen_tn_id or (i + 1)
            # Default read length (can be overridden by config in pipeline)
            read_len = 150
            
            analysis_dir = f"jc_{start}_{end}_{tn_id:03d}_L{read_len}"
            
            # Create candidate record
            candidate = CandidateTnJc2.from_other(
                classified,
                analysis_dir=analysis_dir,
            )
            candidate_records.append(candidate)
        
        result = RecordTypedDf.from_records(candidate_records, CandidateTnJc2)
        info(f"Filtered to {len(result)} candidates (from {len(self.classified_tnjc2)} classified)")
        
        return result

    def _save_output(self, output: RecordTypedDf[CandidateTnJc2]) -> None:
        """Save candidates to CSV."""
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDf[CandidateTnJc2]:
        """Load candidates from CSV."""
        return RecordTypedDf.from_csv(self.output_file, CandidateTnJc2)

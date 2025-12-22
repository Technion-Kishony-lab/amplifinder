"""Step 9: Filter candidates by amplicon length."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDF, ClassifiedTnJc2, CandidateTnJc2,
)
from amplifinder.steps.base import Step
from amplifinder.logger import info


class FilterCandidatesStep(Step[RecordTypedDF[CandidateTnJc2]]):
    """Filter candidates by amplicon length.
    
    Filters out candidates that are:
    - Too short (< min_amplicon_length)
    - Too long (> max_amplicon_length)
    
    Also assigns analysis directory names to each candidate.
    """

    def __init__(
        self,
        classified_tnjc2: RecordTypedDF[ClassifiedTnJc2],
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

    def _calculate_output(self) -> RecordTypedDF[CandidateTnJc2]:
        """Filter candidates by amplicon length."""
        candidate_records = []
        
        for i, classified in enumerate(self.classified_tnjc2):
            # Filter by length
            if (classified.amplicon_length <= self.min_amplicon_length or
                classified.amplicon_length >= self.max_amplicon_length):
                continue
            
            # Generate analysis directory name
            # Format: jc_{start}_{end}_{tn_id:03d}_L{read_len}
            # For now, use index as tn_id and assume read_len=150
            # In practice, read_len should come from config
            start = classified.pos_chr_L
            end = classified.pos_chr_R
            tn_id = classified.chosen_tn_id or (i + 1)
            read_len = 150  # TODO: get from config
            
            analysis_dir = f"jc_{start}_{end}_{tn_id:03d}_L{read_len}"
            
            # Create candidate record
            candidate = CandidateTnJc2.from_other(
                classified,
                analysis_dir=analysis_dir,
            )
            candidate_records.append(candidate)
        
        result = RecordTypedDF.from_records(candidate_records, CandidateTnJc2)
        info(f"Filtered to {len(result)} candidates (from {len(self.classified_tnjc2)} classified)")
        
        return result

    def _save_output(self, output: RecordTypedDF[CandidateTnJc2]) -> None:
        """Save candidates to CSV."""
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[CandidateTnJc2]:
        """Load candidates from CSV."""
        return RecordTypedDF.from_csv(self.output_file, CandidateTnJc2)

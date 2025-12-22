"""Step 14: Export results to CSV."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDF, AnalyzedTnJc2, ISJC2Export
from amplifinder.steps.base import Step
from amplifinder.logger import info


class ExportStep(Step[RecordTypedDF[ISJC2Export]]):
    """Export analyzed candidates to CSV files.
    
    Creates two CSV files:
    1. ISJC2.csv - All analyzed candidates
    2. candidate_amplifications.csv - Filtered candidates based on copy number thresholds
    """

    def __init__(
        self,
        analyzed_candidates: RecordTypedDF[AnalyzedTnJc2],
        output_dir: Path,
        copy_number_threshold: float = 1.5,
        del_copy_number_threshold: float = 0.3,
        filter_amplicon_length: int = 100,
        force: Optional[bool] = None,
    ):
        self.analyzed_candidates = analyzed_candidates
        self.output_dir = Path(output_dir)
        self.copy_number_threshold = copy_number_threshold
        self.del_copy_number_threshold = del_copy_number_threshold
        self.filter_amplicon_length = filter_amplicon_length
        
        self.isjc2_file = output_dir / "ISJC2.csv"
        self.candidates_file = output_dir / "candidate_amplifications.csv"
        
        super().__init__(
            input_files=[],
            output_files=[self.isjc2_file, self.candidates_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDF[ISJC2Export]:
        """Export candidates to CSV."""
        if len(self.analyzed_candidates) == 0:
            info("No candidates to export, creating empty ISJC2.csv with headers")
            export_df = RecordTypedDF.empty(ISJC2Export)
            export_df.to_csv(self.isjc2_file)
            export_df.to_csv(self.candidates_file)
            info(f"Created empty ISJC2.csv and candidate_amplifications.csv with headers")
            return export_df
        
        # Build export records by iterating over typed records
        export_records = []
        for candidate in self.analyzed_candidates:
            export_records.append(ISJC2Export(
                isolate=candidate.iso_name,
                Reference=candidate.ref_name,
                Ancestor=candidate.anc_name,
                Positions_in_chromosome=f"{candidate.pos_chr_L}-{candidate.pos_chr_R}",
                Direction_in_chromosome=f"{candidate.dir_chr_L}/{candidate.dir_chr_R}",
                amplicon_length=candidate.amplicon_length,
                IS_element=','.join(map(str, candidate.tn_ids)) if candidate.tn_ids else None,
                median_copy_number=candidate.amplicon_coverage,
                mode_copy_number=candidate.amplicon_coverage_mode,
                event=candidate.event,
                isolate_architecture=str(candidate.isolate_architecture),
            ))
        
        export_df = RecordTypedDF.from_records(export_records, ISJC2Export)
        
        # Sort by mode_copy_number descending
        export_df = export_df.pipe(
            lambda df: df.sort_values('mode_copy_number', ascending=False)
        )
        
        # Export ISJC2.csv (all candidates)
        export_df.to_csv(self.isjc2_file)
        info(f"Exported {len(export_df)} candidates to {self.isjc2_file}")
        
        # Filter candidates for candidate_amplifications.csv
        filtered = export_df.pipe(
            lambda df: df[
                ((df['mode_copy_number'] > self.copy_number_threshold) |
                 (df['mode_copy_number'] < self.del_copy_number_threshold)) &
                (df['amplicon_length'] > self.filter_amplicon_length)
            ]
        )
        
        filtered.to_csv(self.candidates_file)
        info(f"Exported {len(filtered)} filtered candidates to {self.candidates_file}")
        
        return export_df

    def _save_output(self, output: RecordTypedDF[ISJC2Export]) -> None:
        """Output already saved in _calculate_output."""
        pass

    def load_outputs(self) -> RecordTypedDF[ISJC2Export]:
        """Load exported ISJC2 data."""
        return RecordTypedDF.from_csv(self.isjc2_file, ISJC2Export)

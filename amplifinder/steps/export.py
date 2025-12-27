"""Step 14: Export results to CSV."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, AnalyzedTnJc2, ExportedTnJc2
from amplifinder.steps.base import Step
from amplifinder.logger import info


class ExportTnJc2Step(Step[RecordTypedDf[ExportedTnJc2]]):
    """Export analyzed candidates to CSV files.

    Creates two CSV files:
    1. tnjc2_exported.csv - All analyzed candidates
    2. candidate_amplifications.csv - Filtered candidates based on copy number thresholds
    """

    def __init__(
        self,
        analyzed_tnjc2s: RecordTypedDf[AnalyzedTnJc2],
        output_dir: Path,
        ref_name: str,
        iso_name: str,
        anc_name: Optional[str] = None,
        copy_number_threshold: float = 1.5,
        del_copy_number_threshold: float = 0.3,
        filter_amplicon_length: int = 100,
        force: Optional[bool] = None,
    ):
        self.analyzed_tnjc2s = analyzed_tnjc2s
        self.output_dir = Path(output_dir)
        self.ref_name = ref_name
        self.iso_name = iso_name
        self.anc_name = anc_name
        self.copy_number_threshold = copy_number_threshold
        self.del_copy_number_threshold = del_copy_number_threshold
        self.filter_amplicon_length = filter_amplicon_length

        self.isjc2_file = output_dir / "tnjc2_exported.csv"
        self.candidates_file = output_dir / "candidate_amplifications.csv"

        super().__init__(
            input_files=[],
            output_files=[self.isjc2_file, self.candidates_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[ExportedTnJc2]:
        """Export candidates to CSV."""
        # Build export records by iterating over typed records
        export_records = []
        for analyzed_tnjc2 in self.analyzed_tnjc2s:
            export_records.append(ExportedTnJc2(
                isolate=self.iso_name,
                Reference=self.ref_name,
                Ancestor=self.anc_name,
                Positions_in_chromosome=f"{analyzed_tnjc2.pos_scaf_L}-{analyzed_tnjc2.pos_scaf_R}",
                Direction_in_chromosome=f"{analyzed_tnjc2.dir_scaf_L}/{analyzed_tnjc2.dir_scaf_R}",
                amplicon_length=analyzed_tnjc2.amplicon_length,
                IS_element=','.join(map(str, analyzed_tnjc2.tn_ids)) if analyzed_tnjc2.tn_ids else None,
                median_copy_number=analyzed_tnjc2.amplicon_coverage,
                mode_copy_number=analyzed_tnjc2.amplicon_coverage_mode,
                event=analyzed_tnjc2.event,
                isolate_architecture=str(analyzed_tnjc2.isolate_architecture),
            ))

        export_df = RecordTypedDf.from_records(export_records, ExportedTnJc2)

        # Sort by mode_copy_number descending
        export_df = export_df.pipe(
            lambda df: df.sort_values('mode_copy_number', ascending=False)
        )

        # Export tnjc2_exported.csv (all candidates) - to_csv handles empty DataFrames automatically
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

        # to_csv handles empty DataFrames automatically
        filtered.to_csv(self.candidates_file)
        info(f"Exported {len(filtered)} filtered candidates to {self.candidates_file}")

        return export_df

    def _save_output(self, output: RecordTypedDf[ExportedTnJc2]) -> None:
        """Output already saved in _calculate_output."""
        pass

    def load_outputs(self) -> RecordTypedDf[ExportedTnJc2]:
        """Load exported TnJc2 data."""
        return RecordTypedDf.from_csv(self.isjc2_file, ExportedTnJc2)

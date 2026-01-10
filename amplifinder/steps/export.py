"""Step 14: Export results to CSV."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, ClassifiedTnJc2, ExportedTnJc2, Orientation, Genome
from amplifinder.steps.base import OutputStep


class ExportTnJc2Step(OutputStep[RecordTypedDf[ExportedTnJc2]]):
    """Export analyzed candidates to CSV files.

    Creates two CSV files:
    1. tnjc2_exported.csv - All analyzed candidates
    2. candidate_amplifications.csv - Filtered candidates based on copy number thresholds
    """

    def __init__(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        genome: Genome,
        output_dir: Path,
        ref_name: str,
        iso_name: str,
        anc_name: Optional[str] = None,
        copy_number_threshold: float = 1.5,
        del_copy_number_threshold: float = 0.3,
        filter_amplicon_length: int = 100,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self.genome = genome
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
            output_files=[self.isjc2_file, self.candidates_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[ExportedTnJc2]:
        """Build export dataframe (writing handled in _save_output)."""
        # Build export records by iterating over typed records
        export_records = []
        for tnjc2 in self.classified_tnjc2s:
            # Get span_origin from segment scaffold
            span_origin = tnjc2.span_origin
            export_records.append(ExportedTnJc2(
                isolate=self.iso_name,
                Reference=self.ref_name,
                Ancestor=self.anc_name,
                Positions_in_chromosome=f"{tnjc2.left}-{tnjc2.right}",
                Direction_in_chromosome=(
                    f"{Orientation.REVERSE if span_origin else Orientation.FORWARD}/"
                    f"{Orientation.FORWARD if span_origin else Orientation.REVERSE}"),
                amplicon_length=tnjc2.amplicon_length,
                IS_element=','.join(map(str, tnjc2.tn_ids)) if tnjc2.tn_ids else None,
                median_copy_number=tnjc2.copy_number,
                mode_copy_number=tnjc2.copy_number,
                event=tnjc2.event,
                isolate_architecture=str(tnjc2.isolate_architecture),
            ))

        export_df = RecordTypedDf.from_records(export_records, ExportedTnJc2)

        # Sort by mode_copy_number descending
        export_df = export_df.pipe(
            lambda df: df.sort_values('mode_copy_number', ascending=False)
        )

        return export_df

    def _save_output(self, output: RecordTypedDf[ExportedTnJc2]) -> None:
        """Write export CSVs."""
        # Export tnjc2_exported.csv (all candidates)
        output.to_csv(self.isjc2_file, index=False)

        # Filter candidates for candidate_amplifications.csv
        filtered = self._filter_export_df(output)
        filtered.to_csv(self.candidates_file, index=False)

    def report_output_message(self, output: RecordTypedDf[ExportedTnJc2], *, from_cache: bool) -> Optional[str]:
        filtered_len = len(self._filter_export_df(output))
        return (
            f"Exported {len(output)} candidates to {self.isjc2_file}; "
            f"{filtered_len} filtered candidates to {self.candidates_file}"
        )

    def _filter_export_df(self, df: RecordTypedDf[ExportedTnJc2]) -> RecordTypedDf[ExportedTnJc2]:
        return df.pipe(
            lambda d: d[
                ((d['mode_copy_number'] > self.copy_number_threshold) |
                 (d['mode_copy_number'] < self.del_copy_number_threshold)) &
                (d['amplicon_length'] > self.filter_amplicon_length)
            ]
        )

"""Step: Run breseq alignment pipeline."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import OutputStep
from amplifinder.tools.breseq import run_breseq, parse_breseq_output
from amplifinder.data_types.record_types import BreseqJunction
from amplifinder.data_types.typed_df import RecordTypedDf


class BreseqStep(OutputStep[RecordTypedDf[BreseqJunction]]):
    """Run breseq alignment pipeline.

    Output is breseq run directory (output.gd), not a saved CSV.
    Returns RecordTypedDf[BreseqJunction] from parsed breseq output.
    """

    STEP_LOCK_TIMEOUT = 7200  # breseq can run for hours

    def __init__(
        self,
        breseq_path: Path,
        fastq_path: Optional[Path] = None,
        ref_file: Optional[Path] = None,
        docker: bool = False,
        threads: int = 1,
        force: Optional[bool] = None,
    ):
        self.fastq_path = Path(fastq_path) if fastq_path else None
        self.ref_file = Path(ref_file) if ref_file else None
        self.docker = docker
        self.threads = threads

        # Breseq output location
        self.breseq_path = Path(breseq_path)
        output_gd = self.breseq_path / "output" / "output.gd"

        # If output doesn't exist, fastq_path and ref_file are required
        if not output_gd.exists() and (fastq_path is None or ref_file is None):
            raise ValueError("fastq_path and ref_file are required when output doesn't exist")

        input_files = []
        if fastq_path:
            input_files.append(fastq_path)
        if ref_file:
            input_files.append(ref_file)

        super().__init__(
            input_files=input_files if input_files else None,
            artifact_files=[output_gd],
            force=force,
        )

    def _generate_artifacts(self) -> None:
        """Run breseq to produce output.gd and related files."""
        run_breseq(
            ref_paths=[self.ref_file],
            fastq_path=self.fastq_path,
            output_path=self.breseq_path,
            docker=self.docker,
            threads=self.threads,
        )

    def _calculate_output(self) -> RecordTypedDf[BreseqJunction]:
        """Parse breseq junction output into RecordTypedDf.

        Returns:
            RecordTypedDf[BreseqJunction] with column renaming applied
        """
        outputs = parse_breseq_output(self.breseq_path)
        jc_df = outputs.get("JC")

        if jc_df is None or jc_df.empty:
            # Return empty RecordTypedDf
            return RecordTypedDf.from_records([], BreseqJunction)

        # Rename breseq columns to our naming convention
        jc_df = jc_df.rename(columns={
            'flanking_left': 'flanking1',
            'flanking_right': 'flanking2'
        })

        return RecordTypedDf(jc_df, BreseqJunction)

    def report_output_message(self, output: RecordTypedDf[BreseqJunction], *, from_cache: bool) -> Optional[str]:
        """Report junction count."""
        prefix = 'Artifacts cached,' if from_cache else 'Artifacts created,'
        return f"{prefix} {len(output)} breseq junctions."


class AncBreseqStep(BreseqStep):
    """Run breseq alignment pipeline for ancestor."""
    pass

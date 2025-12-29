"""Step: Run breseq alignment pipeline."""

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.tools.breseq import run_breseq, parse_breseq_output


class BreseqStep(Step[Dict[str, pd.DataFrame]]):
    """Run breseq alignment pipeline."""

    STEP_LOCK_TIMEOUT = 7200  # breseq can run for hours

    def __init__(
        self,
        output_path: Path,
        fastq_path: Optional[Path] = None,
        ref_file: Optional[Path] = None,
        docker: bool = False,
        threads: int = 1,
        force: Optional[bool] = None,
    ):
        self.output_path = Path(output_path)
        self.fastq_path = Path(fastq_path) if fastq_path else None
        self.ref_file = Path(ref_file) if ref_file else None
        self.docker = docker
        self.threads = threads

        # Check if output exists
        output_file = output_path / "output" / "output.gd"
        output_exists = output_file.exists()

        # If output doesn't exist, fastq_path and ref_file are required
        if not output_exists and (fastq_path is None or ref_file is None):
            raise ValueError("fastq_path and ref_file are required when output doesn't exist")

        input_files = []
        if fastq_path:
            input_files.append(fastq_path)
        if ref_file:
            input_files.append(ref_file)

        super().__init__(
            input_files=input_files if input_files else None,
            output_files=[output_file],
            force=force,
        )

    def _output_labels(self) -> list[str]:
        """Human-readable labels for outputs."""
        return [str(self.output_path.name)]

    def _calculate_output(self) -> Dict[str, pd.DataFrame]:
        """Run breseq and return parsed output."""
        run_breseq(
            ref_paths=[self.ref_file],
            fastq_path=self.fastq_path,
            output_path=self.output_path,
            docker=self.docker,
            threads=self.threads,
        )
        return self.load_outputs()

    def load_outputs(self) -> Dict[str, pd.DataFrame]:
        """Parse breseq output into DataFrames.

        Returns:
            Dict with keys: JC, SNP, MOB, DEL, UN
        """
        return parse_breseq_output(self.output_path)


class AncBreseqStep(BreseqStep):
    """Run breseq alignment pipeline for ancestor."""
    pass

"""Step: Run breseq alignment pipeline."""

from pathlib import Path
from typing import Dict, Optional

import pandas as pd

from amplifinder.steps.base import Step
from amplifinder.tools.breseq import run_breseq, parse_breseq_output


class BreseqStep(Step):
    """Run breseq alignment pipeline."""

    def __init__(
        self,
        fastq_path: Path,
        ref_file: Path,
        output_path: Path,
        docker: bool,
        threads: int,
        force: Optional[bool] = None,
    ):
        self.fastq_path = Path(fastq_path)
        self.ref_file = Path(ref_file)
        self.output_path = Path(output_path)
        self.docker = docker
        self.threads = threads

        super().__init__(
            inputs=[fastq_path, ref_file],
            outputs=[output_path / "output" / "output.gd"],
            force=force,
        )

    def _run(self) -> None:
        """Run breseq."""
        run_breseq(
            ref_paths=[self.ref_file],
            fastq_path=self.fastq_path,
            output_path=self.output_path,
            docker=self.docker,
            threads=self.threads,
        )

    def read_outputs(self) -> Dict[str, pd.DataFrame]:
        """Parse breseq output into DataFrames.

        Returns:
            Dict with keys: JC, SNP, MOB, DEL, UN
        """
        return parse_breseq_output(self.output_path)

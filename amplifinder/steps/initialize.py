"""Step 0: Initialize output directories."""

from pathlib import Path

from amplifinder.steps.base import Step
from amplifinder.logger import info


class InitializingStep(Step[Path]):
    """Create output directories for the pipeline."""

    def __init__(self, output_dir: Path, iso_name: str, force: bool = None):
        self.output_dir = Path(output_dir)
        self.iso_name = iso_name
        self.iso_output = self.output_dir / iso_name
        super().__init__(output_files=[self.iso_output], force=force)

    def _calculate_output(self) -> Path:
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.iso_output.mkdir(parents=True, exist_ok=True)
        info(f"Created output directory: {self.iso_output}")
        return self.iso_output

    def load_outputs(self) -> Path:
        return self.iso_output

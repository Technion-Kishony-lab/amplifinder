"""Step 0: Initialize output directories."""

from pathlib import Path
from typing import Optional, Tuple

from amplifinder.steps.base import OutputStep
from amplifinder.config import Config
from amplifinder.utils.file_utils import ensure_dir


class InitializingStep(OutputStep[Tuple[Path, Optional[Path]]]):
    """Create output directories for the pipeline.

    Folder structure: {output_dir}/{ref_name}/{anc_name}/{iso_name}/
    When iso_name == anc_name, creates {output_dir}/{ref_name}/{name}/{name}/
    """

    def __init__(
        self,
        config: Config,
        force: Optional[bool] = None,
    ):
        self.config = config
        self.iso_run_dir = config.iso_run_dir

        # Build output_files list - include ancestor run dir if needed
        output_files = [self.iso_run_dir]
        if config.has_ancestor:
            self.anc_run_dir = config.anc_run_dir
            output_files.append(self.anc_run_dir)
        else:
            self.anc_run_dir = None

        super().__init__(artifact_files=output_files, force=force)

    def _generate_artifacts(self) -> None:
        """Create directories and save config."""
        ensure_dir(self.iso_run_dir)

        # Also create ancestor run directory if ancestor exists (needed for breseq output)
        if self.anc_run_dir is not None:
            ensure_dir(self.anc_run_dir)

        # Save config to run directory
        self.config.save(self.iso_run_dir)

    def _calculate_output(self) -> Tuple[Path, Optional[Path]]:
        """Return run directories."""
        return self.iso_run_dir, self.anc_run_dir

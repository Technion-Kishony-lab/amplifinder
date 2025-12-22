"""Step 0: Initialize output directories."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import Step
from amplifinder.config import Config, save_config, get_run_dir
from amplifinder.logger import info


class InitializingStep(Step[Path]):
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
        self.run_dir = get_run_dir(config)
        super().__init__(output_files=[self.run_dir], force=force)

    def _calculate_output(self) -> Path:
        """Create directories and save config."""
        self.run_dir.mkdir(parents=True, exist_ok=True)
        info(f"Created output directory: {self.run_dir}")
        
        # Save config to run directory
        save_config(self.config, self.run_dir)
        
        return self.run_dir

    def load_outputs(self) -> Path:
        return self.run_dir

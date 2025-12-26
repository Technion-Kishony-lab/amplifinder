"""Step 0: Initialize output directories."""

from pathlib import Path
from typing import Optional, Tuple

from amplifinder.steps.base import Step
from amplifinder.config import Config, save_config, get_iso_run_dir, get_anc_run_dir
from amplifinder.logger import info
from amplifinder.utils.tools import ensure_dir


class InitializingStep(Step[Tuple[Path, Optional[Path]]]):
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
        self.iso_run_dir = get_iso_run_dir(config)
        
        # Build output_files list - include ancestor run dir if needed
        output_files = [self.iso_run_dir]
        if config.has_ancestor:
            self.anc_run_dir = get_anc_run_dir(config)
            output_files.append(self.anc_run_dir)
        else:
            self.anc_run_dir = None
        
        super().__init__(output_files=output_files, force=force)

    def _calculate_output(self) -> Tuple[Path, Optional[Path]]:
        """Create directories and save config."""
        ensure_dir(self.iso_run_dir)
        
        # Also create ancestor run directory if ancestor exists (needed for breseq output)
        if self.anc_run_dir is not None:
            ensure_dir(self.anc_run_dir)
        
        # Save config to run directory
        save_config(self.config, self.iso_run_dir)
        
        return self.iso_run_dir, self.anc_run_dir

    def load_outputs(self) -> Tuple[Path, Optional[Path]]:
        return self.iso_run_dir, self.anc_run_dir

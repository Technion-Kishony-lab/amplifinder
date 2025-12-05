"""Step: Fetch/load reference genome."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import Step
from amplifinder.data_types.genome import Genome, get_genome


class GetReferenceStep(Step):
    """Fetch/load reference genome from NCBI or local cache."""
    
    def __init__(
        self,
        ref_name: str,
        ref_path: Path,
        ncbi: bool = True,
        force: Optional[bool] = None,
    ):
        self.ref_name = ref_name
        self.ref_path = Path(ref_path)
        self.ncbi = ncbi
        
        # Output: mapping file
        self.mapping_file = self.ref_path / f"{ref_name}.json"
        
        super().__init__(inputs=[], outputs=[self.mapping_file], force=force)

    def has_output(self) -> bool:
        """Check if output exists."""
        if not super().has_output():
            return False
        try:
            self.read_outputs()
            return True
        except FileNotFoundError:
            return False
    
    def _run(self) -> None:
        """Fetch genome from NCBI or load from cache."""
        get_genome(self.ref_name, self.ref_path, self.ncbi)
    
    def read_outputs(self) -> Genome:
        """Load genome from cached files."""
        return get_genome(self.ref_name, self.ref_path, ncbi=False)
    
    def run_and_read_outputs(self) -> Genome:
        """Run step and return genome."""
        return super().run_and_read_outputs()

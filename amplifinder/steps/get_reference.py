"""Step: Fetch/load reference genome."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import Step
from amplifinder.data_types.genome import Genome, get_genome


class GetRefGenomeStep(Step[Genome]):
    """Fetch/load reference genome from NCBI or local cache."""

    def __init__(
        self,
        ref_name: str,
        ref_path: Path,
        ncbi: bool,
        force: Optional[bool] = None,
    ):
        self.ref_name = ref_name
        self.ref_path = Path(ref_path)
        self.ncbi = ncbi

        # Output: mapping file
        self.mapping_file = self.ref_path / f"{ref_name}.json"

        super().__init__(input_files=[], output_files=[self.mapping_file], force=force)

    def has_output_files(self) -> bool:
        """Check if output exists and is valid."""
        if not super().has_output_files():
            return False
        try:
            self.read_outputs()
            return True
        except FileNotFoundError:
            return False

    def _calculate_output(self) -> Genome:
        """Fetch genome from NCBI or load from cache."""
        return get_genome(self.ref_name, self.ref_path, self.ncbi)

    def load_outputs(self) -> Genome:
        """Load genome from cached files."""
        return get_genome(self.ref_name, self.ref_path, ncbi=False)

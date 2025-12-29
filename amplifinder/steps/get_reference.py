"""Step: Fetch/load reference genome."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import Step
from amplifinder.data_types.genome import Genome, get_genome, exists_genome


class GetRefGenomeStep(Step[Genome]):
    """Fetch/load reference genome from NCBI or local cache.

    Uses file locking to prevent race conditions when multiple
    parallel runs try to download the same reference genome.
    """

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

        super().__init__(output_files=[self.mapping_file], force=force)

    def has_output_files(self) -> bool:
        """Check if output exists and is valid."""
        if not super().has_output_files():
            return False
        return exists_genome(self.ref_name, self.ref_path)

    def _calculate_output(self) -> Genome:
        """Fetch genome from NCBI or load from cache.

        Locking is handled by Step.run() via _get_lock_target.
        """
        self.log(f"Fetching reference genome: {self.ref_name}")
        return get_genome(self.ref_name, self.ref_path, self.ncbi)

    def _get_lock_target(self) -> Path:
        """Lock on shared genome directory to serialize cache writes."""
        return self.ref_path / ".lock_genome"

    def load_outputs(self) -> Genome:
        """Load genome from cached files."""
        return get_genome(self.ref_name, self.ref_path, ncbi=False)

    def read_outputs(self) -> Genome:
        """Alias for load_outputs (used by has_output_files validation)."""
        return self.load_outputs()

"""Step: Fetch/load reference genome."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import Step
from amplifinder.data_types.genome import Genome, get_genome, exists_genome
from amplifinder.utils.file_lock import locked_resource


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

        Uses locking to prevent parallel downloads of the same genome.
        The lock is acquired around the entire fetch operation.
        """
        # Lock around genome fetch to prevent parallel downloads
        # The base Step.run() already handles locking for the output file,
        # but we add an additional lock here for the genome directory
        # in case multiple runs with different output dirs share ref_path
        with locked_resource(self.ref_path, f"genome_{self.ref_name}", timeout=600):
            # Re-check if genome exists (another process may have downloaded it)
            if self.has_output_files():
                self.log(f"{self.name}: genome already cached (verified under lock)")
                return self.load_outputs()

            self.log(f"Fetching reference genome: {self.ref_name}")
            return get_genome(self.ref_name, self.ref_path, self.ncbi)

    def load_outputs(self) -> Genome:
        """Load genome from cached files."""
        return get_genome(self.ref_name, self.ref_path, ncbi=False)

    def read_outputs(self) -> Genome:
        """Alias for load_outputs (used by has_output_files validation)."""
        return self.load_outputs()

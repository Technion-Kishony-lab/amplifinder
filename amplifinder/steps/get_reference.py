"""Step: Fetch/load reference genome."""

from pathlib import Path
from typing import Optional

from amplifinder.steps.base import OutputStep
from amplifinder.data_types.genome import Genome, get_genome, exists_genome


class GetRefGenomeStep(OutputStep[Genome]):
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

        super().__init__(artifact_files=[self.mapping_file], force=force)

    def _get_lock_target(self) -> Path:
        """Lock on shared genome directory to serialize cache writes."""
        return self.ref_path

    def has_artifact_files(self) -> bool:
        """Check if genome mapping file exists and is valid."""
        return exists_genome(self.ref_name, self.ref_path)

    def _generate_artifacts(self) -> None:
        """Fetch genome from NCBI or load from cache (creates mapping file)."""
        get_genome(self.ref_name, self.ref_path, self.ncbi)

    def _calculate_output(self) -> Genome:
        """Load genome from cached files."""
        return get_genome(self.ref_name, self.ref_path, ncbi=False)

"""Global environment configuration.

Loads server-specific paths from config.yaml in project root.
"""

from pathlib import Path
from typing import Optional
import yaml

from amplifinder.utils.flag_utils import MutableFlag


# Search paths for config.yaml
CONFIG_SEARCH_PATHS = [
    Path(__file__).parent.parent / "config.yaml",  # Project root
    Path.home() / ".amplifinder" / "config.yaml",  # User home
    Path("/etc/amplifinder/config.yaml"),          # System-wide
]


def _find_config() -> Optional[Path]:
    """Find config.yaml in standard locations."""
    for path in CONFIG_SEARCH_PATHS:
        if path.exists():
            return path
    return None


def _load_config() -> dict:
    """Load global config from config.yaml."""
    config_path = _find_config()
    if config_path is None:
        return {}

    with open(config_path) as f:
        return yaml.safe_load(f) or {}


# Load config at import time
_CONFIG = _load_config()


# Export paths as module-level variables
def _get_path(key: str) -> Optional[Path]:
    """Get path from config."""
    return Path(_CONFIG[key]) if _CONFIG.get(key) is not None else None


BLAST_PATH: Optional[Path] = _get_path("blast_path")
SAMTOOLS_PATH: Optional[Path] = _get_path("samtools_path")
BOWTIE2_PATH: Optional[Path] = _get_path("bowtie2_path")
BRESEQ_DOCKER: bool = _CONFIG.get("breseq_docker", False)  # TODO: This is not used yet

# Global parameters for read type classification
IGNORE_DUPLICATES: bool = True  # If False, keep all duplicate read_ids as a list
LINK_PAIRED_END: bool = True  # If False, don't merge LEFT and RIGHT reads into PAIRED

# Maximum distance from junction to allow indels 
# (None = no limit, -1 = no indels, 0 = only indels precisely at junction, >0 = maximum distance)
ALLOW_INDELS_AT_JUNCTION_DISTANCE: Optional[int] = 4

DEBUG = MutableFlag(False)

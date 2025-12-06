"""Global environment configuration.

Loads server-specific paths from config.yaml in project root.
"""

from pathlib import Path
from typing import Optional
import yaml


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
ISDB_PATH: Optional[Path] = _get_path("isdb_path")
BRESEQ_DOCKER: bool = _CONFIG.get("breseq_docker", False)

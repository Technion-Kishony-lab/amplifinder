"""ISfinder database utilities."""

from pathlib import Path


def get_builtin_isfinder_db_path() -> Path:
    """Get path to bundled ISfinder database (IS.fna)."""
    return Path(__file__).parent / "ISfinderDB" / "IS.fna"

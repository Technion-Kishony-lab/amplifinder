"""ISfinder database utilities."""

from pathlib import Path
import os
import shutil

from amplifinder.utils.file_utils import is_writable_dir


def get_builtin_isfinder_db_path() -> Path:
    """Get path to bundled ISfinder database (IS.fna).

    If the bundled directory is not writable, copy it to a user-writable cache
    location and return the cached copy instead.
    """
    bundled = Path(__file__).parent / "ISfinderDB"
    is_fna = bundled / "IS.fna"
    if is_writable_dir(bundled):
        return is_fna

    cache_root = Path(
        os.environ.get(
            "AMPLIFINDER_ISFINDER_CACHE",
            str(Path.home() / ".amplifinder" / "ISfinderDB"),
        )
    )
    cache_root.mkdir(parents=True, exist_ok=True)
    cached_is_fna = cache_root / "IS.fna"
    if not cached_is_fna.exists():
        shutil.copytree(bundled, cache_root, dirs_exist_ok=True)
    return cached_is_fna

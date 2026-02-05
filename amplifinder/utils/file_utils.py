"""File and directory utilities."""

import os
import shutil
from pathlib import Path
from typing import Union


def fmt_count(num: int, total: int | None = None) -> str:
    """Format count with optional percentage.

    Args:
        num: Number to format
        total: Optional total for percentage calculation

    Returns:
        Formatted string: 'num:6' or 'num:6 (pct:4.1f%)' if total provided

    Examples:
        >>> fmt_count(100)
        '   100'
        >>> fmt_count(75, 100)
        '    75 (  75.0%)'
    """
    if total is None or total == 0:
        return f"{num:6}"
    pct = num / total * 100
    return f"{num:6} ({pct:4.1f}%)"


def ensure_dir(path: Union[str, Path], cleanup: bool = False) -> Path:
    """Create directory and all parent directories if they don't exist.

    Args:
        path: Directory path to create (str or Path)
        cleanup: If True, remove existing directory/file before creating

    Returns:
        Path object of the created directory
    """
    path = Path(path)
    if cleanup and path.exists():
        remove_file_or_dir(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def ensure_parent_dir(path: Union[str, Path]) -> Path:
    """Create parent directory of a file/directory path if it doesn't exist.

    Args:
        path: File or directory path whose parent should be created (str or Path)

    Returns:
        Path object of the input path (not the parent)
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def remove_file_or_dir(path: Union[str, Path]) -> None:
    """Remove file or directory if it exists.

    Args:
        path: File or directory path to remove (str or Path)
    """
    path = Path(path)
    if path.exists():
        if path.is_dir():
            shutil.rmtree(path)
        else:
            path.unlink()


def is_writable_dir(path: Path) -> bool:
    """Check if a directory is writable by the current user."""
    return path.exists() and os.access(path, os.W_OK)

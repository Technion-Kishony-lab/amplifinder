"""File-based locking for parallel run safety.

Prevents race conditions when multiple AmpliFinder processes access
shared resources. Uses filelock for cross-platform file locks.
"""

from contextlib import contextmanager
from pathlib import Path
from typing import Optional
import hashlib

from filelock import FileLock, Timeout

from amplifinder.logger import logger
from amplifinder.utils.file_utils import is_writable_dir


DEFAULT_LOCK_TIMEOUT = 3600  # 1 hour


class LockTimeout(Exception):
    """Lock acquisition timed out."""
    pass


@contextmanager
def _locked_operation(
    lock_filepath: Path,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
    description: str = "operation",
):
    """File-based lock using filelock."""
    lock_filepath = Path(lock_filepath)
    lock_filepath.parent.mkdir(parents=True, exist_ok=True)
    lock = FileLock(lock_filepath)
    try:
        lock.acquire(timeout=timeout)
        logger.debug_message(
            f"Lock acquired for {description}: {lock_filepath}",
            category="file_lock",
            max_prints=None
        )
    except Timeout as exc:
        raise LockTimeout(
            f"Could not acquire lock for {description} within {timeout}s. "
            f"Lock file: {lock_filepath}"
        ) from exc

    try:
        yield
    finally:
        if lock.is_locked:
            lock.release()
            logger.debug_message(
                f"Lock released for {description}",
                category="file_lock",
                max_prints=None
            )


@contextmanager
def locked_resource(
    resource_path: Optional[Path],
    resource_type: str,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
):
    """Lock a shared resource for exclusive access.

    If resource_path is None, this becomes a no-op context manager (no locking).
    This allows conditional locking without code duplication.
    """
    if resource_path is None:
        # No lock needed - just yield
        yield
        return
    resource_path = Path(resource_path)
    if resource_path.is_dir():
        lock_filepath = resource_path / f".{resource_type}.lock"
    else:
        lock_filepath = resource_path.parent / f".{resource_path.stem}.{resource_type}.lock"

    if not is_writable_dir(lock_filepath.parent):
        lock_filepath = _get_fallback_lock_path(resource_path, resource_type)
        logger.warning(
            "Lock path is not writable; using fallback lock file: "
            f"{lock_filepath}"
        )
    with _locked_operation(lock_filepath, timeout, f"{resource_type} at {resource_path}"):
        yield


def _get_fallback_lock_path(resource_path: Path, resource_type: str) -> Path:
    """Fallback lock path in user space for read-only resources."""
    lock_root = Path.home() / ".amplifinder" / "locks"
    lock_root.mkdir(parents=True, exist_ok=True)
    resource_id = hashlib.sha256(str(resource_path).encode()).hexdigest()[:12]
    name = f"{resource_path.name}.{resource_type}.{resource_id}.lock"
    return lock_root / name

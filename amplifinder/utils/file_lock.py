"""File-based locking for parallel run safety.

Prevents race conditions when multiple AmpliFinder processes access
shared resources. Uses filelock for cross-platform file locks.
"""

from contextlib import contextmanager
from pathlib import Path
from amplifinder.env import DEBUG

from filelock import FileLock, Timeout

from amplifinder.logger import logger


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
        if DEBUG:
            logger.debug(f"Lock acquired for {description}: {lock_filepath}")
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
            if DEBUG:
                logger.debug(f"Lock released for {description}")


@contextmanager
def locked_resource(
    resource_path: Path,
    resource_type: str,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
):
    """Lock a shared resource for exclusive access."""
    if resource_path.is_dir():
        lock_filepath = resource_path / f".{resource_type}.lock"
    else:
        lock_filepath = resource_path.parent / f".{resource_path.stem}.{resource_type}.lock"
    with _locked_operation(lock_filepath, timeout, f"{resource_type} at {resource_path}"):
        yield

"""Folder-based locking for parallel run safety.

Prevents race conditions when multiple AmpliFinder processes access 
shared resources. Uses atomic mkdir() - no external dependencies.
"""

from contextlib import contextmanager
from pathlib import Path
from typing import Generator
import time

from amplifinder.logger import debug


DEFAULT_LOCK_TIMEOUT = 3600  # 1 hour
POLL_INTERVAL = 1.0


class LockTimeout(Exception):
    """Lock acquisition timed out."""
    pass


@contextmanager
def locked_operation(
    lock_path: Path,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
    description: str = "operation",
) -> Generator[None, None, None]:
    """Folder-based lock using atomic mkdir()."""
    lock_path = Path(lock_path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    
    deadline = time.time() + timeout
    while True:
        try:
            lock_path.mkdir(exist_ok=False)
            debug(f"Lock acquired for {description}: {lock_path}")
            break
        except FileExistsError:
            if time.time() >= deadline:
                raise LockTimeout(
                    f"Could not acquire lock for {description} within {timeout}s. "
                    f"Lock folder: {lock_path}"
                )
            time.sleep(POLL_INTERVAL)
    
    try:
        yield
    finally:
        lock_path.rmdir()
        debug(f"Lock released for {description}")


def get_step_lock_path(output_file: Path, step_name: str) -> Path:
    """Get lock file path for a step based on its output file.
    
    Lock files are placed in the same directory as the output,
    with a hidden name based on the step name.
    
    Args:
        output_file: Primary output file of the step
        step_name: Name of the step (used in lock filename)
        
    Returns:
        Path to the lock file
    """
    return output_file.parent / f".{step_name}.lock"


def get_resource_lock_path(resource_path: Path, resource_type: str) -> Path:
    """Get lock file path for a shared resource.
    
    Args:
        resource_path: Path to the resource (file or directory)
        resource_type: Type of resource (e.g., 'blast_db', 'genome')
        
    Returns:
        Path to the lock file
    """
    if resource_path.is_dir():
        return resource_path / f".{resource_type}.lock"
    else:
        return resource_path.parent / f".{resource_path.stem}.{resource_type}.lock"


@contextmanager
def locked_resource(
    resource_path: Path,
    resource_type: str,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
) -> Generator[None, None, None]:
    """Lock a shared resource for exclusive access."""
    lock_path = get_resource_lock_path(resource_path, resource_type)
    with locked_operation(lock_path, timeout, f"{resource_type} at {resource_path}"):
        yield

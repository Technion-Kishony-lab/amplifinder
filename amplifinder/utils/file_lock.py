"""File locking utilities for parallel run safety.

Provides mechanisms to prevent race conditions when multiple AmpliFinder
processes access shared resources (reference genomes, BLAST DBs, ancestor runs).
"""

from contextlib import contextmanager
from pathlib import Path
from typing import Generator, Optional

from filelock import FileLock, Timeout

from amplifinder.logger import debug, info


# Default timeout for locks (seconds)
DEFAULT_LOCK_TIMEOUT = 3600  # 1 hour (breseq can take a while)


@contextmanager
def locked_operation(
    lock_path: Path,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
    description: str = "operation",
) -> Generator[None, None, None]:
    """Context manager for file-locked operations.
    
    Acquires an exclusive lock before yielding, releases on exit.
    Other processes attempting to acquire the same lock will block.
    
    Args:
        lock_path: Path to the lock file (will be created if doesn't exist)
        timeout: Maximum seconds to wait for lock (-1 for infinite)
        description: Description for logging
        
    Yields:
        None (use as context manager)
        
    Raises:
        Timeout: If lock cannot be acquired within timeout
        
    Example:
        >>> with locked_operation(Path("/tmp/myfile.lock"), description="download"):
        ...     download_file()
    """
    lock_path = Path(lock_path)
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    
    lock = FileLock(str(lock_path), timeout=timeout)
    
    debug(f"Acquiring lock for {description}: {lock_path}")
    try:
        with lock:
            debug(f"Lock acquired for {description}")
            yield
    except Timeout:
        raise Timeout(
            f"Could not acquire lock for {description} within {timeout}s. "
            f"Another process may be running. Lock file: {lock_path}"
        )
    finally:
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
    """Lock a shared resource for exclusive access.
    
    Args:
        resource_path: Path to the resource
        resource_type: Type of resource for lock naming
        timeout: Lock timeout in seconds
        
    Yields:
        None (use as context manager)
    """
    lock_path = get_resource_lock_path(resource_path, resource_type)
    with locked_operation(lock_path, timeout, f"{resource_type} at {resource_path}"):
        yield


class DoneMarker:
    """Marker file to indicate a multi-file operation completed successfully.
    
    Used to safely detect completion of operations that produce multiple outputs.
    The marker is only written after all outputs are verified.
    """
    
    def __init__(self, directory: Path, operation_name: str = "run"):
        """Initialize done marker.
        
        Args:
            directory: Directory where marker file will be placed
            operation_name: Name used in marker filename
        """
        self.directory = Path(directory)
        self.marker_path = self.directory / f".{operation_name}.done"
    
    def exists(self) -> bool:
        """Check if operation completed successfully."""
        return self.marker_path.exists()
    
    def mark_done(self) -> None:
        """Mark operation as complete."""
        self.directory.mkdir(parents=True, exist_ok=True)
        self.marker_path.touch()
        debug(f"Marked done: {self.marker_path}")
    
    def clear(self) -> None:
        """Remove done marker (for re-runs)."""
        if self.marker_path.exists():
            self.marker_path.unlink()
            debug(f"Cleared done marker: {self.marker_path}")


@contextmanager
def locked_step_execution(
    output_path: Optional[Path],
    step_name: str,
    timeout: int = DEFAULT_LOCK_TIMEOUT,
) -> Generator[bool, None, None]:
    """Lock-protected step execution with re-check pattern.
    
    Implements the double-check locking pattern:
    1. Acquire lock
    2. Re-check if output exists (another process may have created it)
    3. Yield whether execution is needed
    4. Release lock
    
    Args:
        output_path: Primary output path (used for lock location)
        step_name: Name of the step
        timeout: Lock timeout
        
    Yields:
        bool: True if execution is needed, False if output already exists
        
    Example:
        >>> with locked_step_execution(output_csv, "MyStep") as needs_run:
        ...     if needs_run:
        ...         result = compute_expensive_result()
        ...         save_result(result)
    """
    if output_path is None:
        # No file output, no locking needed
        yield True
        return
    
    lock_path = get_step_lock_path(output_path, step_name)
    
    with locked_operation(lock_path, timeout, f"step {step_name}"):
        # Re-check under lock: another process may have created output
        if output_path.exists():
            info(f"{step_name}: output exists (checked under lock), skipping")
            yield False
        else:
            yield True


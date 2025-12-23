"""Pipeline step base class with caching logic."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, List, Optional, TypeVar
import shutil

from amplifinder.logger import info, debug
from amplifinder.utils.file_lock import locked_operation, get_step_lock_path

T = TypeVar("T")

# Default lock timeout for steps (seconds)
STEP_LOCK_TIMEOUT = 7200  # 2 hours (breseq can be very slow)


def _delete_file_or_dir_if_exists(p: Path) -> None:
    """Delete file or directory if it exists."""
    if p.exists():
        if p.is_dir():
            shutil.rmtree(p)
        else:
            p.unlink()


class Step(ABC, Generic[T]):
    """Base class for pipeline steps with input/output file tracking.

    Handles:
    - Skip if all outputs exist (unless force=True)
    - Clean partial outputs before re-run
    - Global and step-specific force control
    """

    # Global force flag (applies to all steps)
    global_force: bool = False
    # Global verbose flag (applies to all steps)
    global_verbose: bool = False

    def __init__(
        self,
        input_files: Optional[List[Path]] = None,
        output_files: Optional[List[Path]] = None,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            input_files: Required input files/dirs (must exist)
                    None - no file inputs.

            output_files: Output files/dirs (produced by this step)
                     None - no file outputs (only calculate result in memory)
            force: Step-specific force flag (None = use global)
        """
        self.input_files: List[Path] = [Path(p) for p in input_files] if input_files else []
        self.output_files: Optional[List[Path]] = [Path(p) for p in output_files] if output_files else None
        self._force = force
        self.run_count = 0

    @property
    def name(self) -> str:
        """Step name (class name by default)."""
        return self.__class__.__name__

    @property
    def force(self) -> bool:
        """Effective force flag (step-specific overrides global)."""
        if self._force is not None:
            return self._force
        return Step.global_force

    @abstractmethod
    def _calculate_output(self) -> T:
        """Execute the step logic. Override in subclass."""
        pass

    def _save_output(self, output: T) -> None:
        """Save output to files. Override in subclass to save to files."""
        pass

    def _save_output_and_verify(self, output: T) -> None:
        """Save output to files and verify that it was created."""
        self._save_output(output)
        if missing_out := self.missing_output_files():
            raise FileNotFoundError(f"{self.name}: expected outputs not created: {missing_out}")

    def load_outputs(self) -> T:
        """Load outputs from files. Override in subclass to return typed data."""
        raise NotImplementedError

    def missing_input_files(self) -> list[Path]:
        """Check if all inputs exist."""
        return [p for p in self.input_files if not p.exists()]

    def missing_output_files(self) -> Optional[list[Path]]:
        """Check if all output files exist. None if no output files."""
        if self.output_files is None:
            return None
        return [p for p in self.output_files if not p.exists()]

    def has_output_files(self) -> bool:
        """True if output_files defined and all exist."""
        missing = self.missing_output_files()
        return missing is not None and len(missing) == 0

    def run(self) -> T:
        """Execute step with caching logic and parallel-safe locking.
        
        Uses double-check locking pattern to prevent race conditions:
        1. Quick check without lock (fast path for cached results)
        2. Acquire lock
        3. Re-check under lock (another process may have created output)
        4. Execute if still needed
        5. Release lock
        """
        # Verbose reporting
        if self.global_verbose:
            if self.output_files:
                output_str = ", ".join(str(p) for p in self.output_files)
                info(f"{self.name}: running (outputs: {output_str})")
            else:
                info(f"{self.name}: running (no file outputs)")
        
        # Check inputs exist
        if missing_input := self.missing_input_files():
            raise FileNotFoundError(f"{self.name}: missing inputs: {missing_input}")

        # Fast path: check if can skip without lock (common case)
        if not self.force and self.has_output_files():
            info(f"{self.name}: skipped (outputs exist)")
            return self.load_outputs()

        # Steps without file outputs don't need locking
        if self.output_files is None:
            return self._run_unlocked()

        # Acquire lock and re-check (TOCTOU fix)
        lock_path = self._get_lock_path()
        with locked_operation(lock_path, STEP_LOCK_TIMEOUT, f"step {self.name}"):
            # Re-check under lock: another process may have created output
            if not self.force and self.has_output_files():
                info(f"{self.name}: skipped (outputs exist, verified under lock)")
                return self.load_outputs()
            
            return self._run_unlocked()

    def _run_unlocked(self) -> T:
        """Execute step logic (assumes lock is held or not needed)."""
        # Clean partial outputs
        self._clean_outputs()

        # Run
        info(f"{self.name}: running...")
        self.run_count += 1
        output = self._calculate_output()

        # Save output (only if output_files defined)
        if self.output_files is not None:
            self._save_output_and_verify(output)

        return output

    def _get_lock_path(self) -> Path:
        """Get lock file path for this step."""
        if self.output_files:
            return get_step_lock_path(self.output_files[0], self.name)
        # Fallback for steps without output files (shouldn't happen in locked path)
        raise ValueError(f"{self.name}: cannot get lock path without output files")

    def _clean_outputs(self) -> None:
        """Remove existing outputs before re-run."""
        if self.output_files is None:
            return
        for p in self.output_files:
            _delete_file_or_dir_if_exists(p)

    @classmethod
    def set_global_force(cls, force: bool) -> None:
        """Set global force flag for all steps."""
        cls.global_force = force

    @classmethod
    def set_global_verbose(cls, verbose: bool) -> None:
        """Set global verbose flag for all steps."""
        cls.global_verbose = verbose

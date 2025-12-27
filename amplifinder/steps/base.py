"""Pipeline step base class with caching logic."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, List, Optional, TypeVar, get_args, get_origin, Type

from amplifinder.logger import info
from amplifinder.utils.file_lock import locked_operation, get_step_lock_path
from amplifinder.utils.file_utils import remove_file_or_dir, ensure_dir
from amplifinder.data_types.typed_df import RecordTypedDf
from amplifinder.data_types.records import Record
from amplifinder.steps.io_naming import default_path

T = TypeVar("T")

# Default lock timeout for steps (seconds)
STEP_LOCK_TIMEOUT = 7200  # 2 hours (breseq can be very slow)


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
    # Should save output flag (False = do not save output)
    should_save: bool = True

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

        # Save output (only if output_files defined and should_save is True)
        if self.output_files is not None and self.should_save:
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
            remove_file_or_dir(p)

    @classmethod
    def set_global_force(cls, force: bool) -> None:
        """Set global force flag for all steps."""
        cls.global_force = force

    @classmethod
    def set_global_verbose(cls, verbose: bool) -> None:
        """Set global verbose flag for all steps."""
        cls.global_verbose = verbose


R = TypeVar("R", bound=Record)


class RecordTypedDfStep(Step[RecordTypedDf[R]], Generic[R]):
    """Base class for steps that output RecordTypedDf to CSV.

    Automatically handles:
    - Output file path from io_naming.default_path()
    - CSV save/load using RecordTypedDf

    Subclasses should:
    - Set class var `record_cls` (or it will be auto-deduced from typing)
    - Override `_calculate_output()` to return RecordTypedDf[R]
    """

    record_cls: Optional[Type[R]] = None  # Can be set explicitly or auto-deduced

    def __init__(
        self,
        output_dir: Optional[Path] = None,
        output_file: Optional[Path] = None,
        input_files: Optional[List[Path]] = None,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            output_dir: Directory for output CSV file (uses default filename from io_naming)
            output_file: Full path to output CSV file (overrides output_dir)
            input_files: Required input files/dirs
            force: Step-specific force flag
        """
        if output_dir is not None and output_file is not None:
            raise ValueError("Cannot specify both output_dir and output_file")
        if output_dir is None and output_file is None:
            raise ValueError("Must specify either output_dir or output_file")

        if output_file is not None:
            # Use provided output_file directly
            self.output_file = Path(output_file)
            self.output_dir = self.output_file.parent
        else:
            # Use output_dir with default filename
            self.output_dir = Path(output_dir)
            record_type = self._get_record_cls()
            self.output_file = default_path(self.output_dir, record_type)

        super().__init__(
            input_files=input_files,
            output_files=[self.output_file],
            force=force,
        )

    @classmethod
    def _get_record_cls(cls) -> Type[R]:
        """Get record class from class var or auto-deduce from typing."""
        # Check class var first
        if cls.record_cls is not None:
            return cls.record_cls

        # Auto-deduce from RecordTypedDfStep[R] typing
        # Look for RecordTypedDfStep[...] in __orig_bases__
        if hasattr(cls, '__orig_bases__'):
            for base in cls.__orig_bases__:
                origin = get_origin(base)
                # Check if it's RecordTypedDfStep[...]
                if origin is RecordTypedDfStep or (
                    hasattr(
                        RecordTypedDfStep,
                        '__origin__') and origin == RecordTypedDfStep.__origin__):
                    args = get_args(base)
                    if args:
                        return args[0]  # R

        raise ValueError(
            f"{cls.__name__}: cannot deduce record_cls. "
            "Set record_cls class var or use RecordTypedDfStep[RecordType] typing."
        )

    def _save_output(self, output: RecordTypedDf[R]) -> None:
        """Save RecordTypedDf to CSV."""
        ensure_dir(self.output_file.parent)
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDf[R]:
        """Load RecordTypedDf from CSV."""
        record_type = self._get_record_cls()
        return RecordTypedDf.from_csv(self.output_file, record_type)

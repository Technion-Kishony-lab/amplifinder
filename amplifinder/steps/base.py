"""Pipeline step base class with caching logic."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, List, Optional, TypeVar, get_args, get_origin, Type

from amplifinder.logger import info
from amplifinder.utils.file_lock import locked_resource
from amplifinder.utils.file_utils import remove_file_or_dir, ensure_dir
from amplifinder.data_types.typed_df import RecordTypedDf
from amplifinder.data_types.records import Record
from amplifinder.steps.io_naming import default_path

T = TypeVar("T")


class Step(ABC, Generic[T]):
    """Base class for pipeline steps with input/output file tracking.

    Handles:
    - Skip if all outputs exist (unless force=True)
    - Clean partial outputs before re-run
    - Global and step-specific force control
    """

    # Default lock timeout for steps (seconds); override per step if needed
    STEP_LOCK_TIMEOUT = 1800

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

    def _get_header(self) -> dict[str, str]:
        """Build logger extras with optional header color."""
        header = self.name
        if self.global_verbose and self.output_files is not None:
            labels = self._output_labels()
            output_str = ", ".join(labels) if labels else "none"
            header = f"{self.name} (outputs: {output_str})"
        header_color = f"\033[36m{header}\033[0m"
        return {"header": f"{header} ", "header_color": f"{header_color} "}

    def log(self, msg: str, *, verbose_only: bool = True, extra: Optional[dict[str, str]] = None) -> None:
        """Step-aware info log that honors global verbosity."""
        if verbose_only and not self.global_verbose:
            return
        info(msg, extra=extra)

    def _output_labels(self) -> list[str]:
        """Human-readable labels for outputs (override for custom logging)."""
        return [str(p.name) for p in self.output_files]

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

    @property
    def is_saving(self) -> bool:
        """True if output_files defined and should_save is True."""
        return self.output_files is not None and len(self.output_files) > 0 and self.should_save

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
        header_extra = self._get_header()

        # Check inputs exist
        if missing_input := self.missing_input_files():
            raise FileNotFoundError(f"{self.name}: missing inputs: {missing_input}")

        # Fast path: check if can skip without lock (common case)
        if not self.force and self.has_output_files():
            self.log("skipped (outputs exist)", extra=header_extra)
            return self.load_outputs()

        # Steps without file outputs or not saving don't need locking
        if not self.is_saving:
            self.log("running...", extra=header_extra)
            return self._run_unlocked()

        # Acquire lock and re-check (TOCTOU fix)
        lock_target = self._get_lock_target()
        with locked_resource(lock_target, self.name, timeout=self.STEP_LOCK_TIMEOUT):
            # Re-check under lock: another process may have created output
            if not self.force and self.has_output_files():
                self.log("skipped (outputs exist, verified under lock)", extra=header_extra)
                return self.load_outputs()
            self.log("running (under lock)", extra=header_extra)
            return self._run_unlocked()

    def _get_lock_target(self) -> Path:
        """Path used for step lock; override to customize lock scope."""
        return self.output_files[0]

    def _run_unlocked(self) -> T:
        """Execute step logic (assumes lock is held or not needed)."""
        # Clean partial outputs
        self._clean_outputs()

        # Run
        self.run_count += 1
        output = self._calculate_output()

        # Save output (only if output_files defined and should_save is True)
        if self.is_saving:
            self._save_output_and_verify(output)

        return output

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

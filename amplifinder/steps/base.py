"""Pipeline step base class with caching logic."""
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, List, Optional, TypeVar, get_args, get_origin, Type

from line_profiler import LineProfiler

from amplifinder.logger import logger
from amplifinder.utils.file_lock import locked_resource
from amplifinder.utils.file_utils import remove_file_or_dir, ensure_dir
from amplifinder.utils.timing import timer
from amplifinder.data_types import RecordTypedDf
from amplifinder.records.base_records import Record
from amplifinder.steps.io_naming import default_filename
from amplifinder.exceptions import PrematureTerminationError


T = TypeVar("T")


class Step(ABC):
    """Pipeline step that creates artifacts via external tools.

    Base class for steps that:
    - Generate artifact files (BAMs, indexes, databases, etc.)
    - Use caching: skip generation if artifacts exist and not force
    - Don't return data to next step (use OutputStep for that)

    For steps that also return data to next step, use OutputStep.
    """

    # Default lock timeout for steps (seconds); override per step if needed
    STEP_LOCK_TIMEOUT = 1800

    # Global force flag (applies to all steps)
    global_force: bool = False

    # Global profile flag (applies to all steps)
    should_profile: bool = False

    def __init__(
        self, *,
        input_files: Optional[List[Path]] = None,
        artifact_files: Optional[List[Path]] = None,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            input_files: Required input files/dirs (must exist)
                    None - no file inputs.

            artifact_files: Files/dirs created by external tools. Generation
                    can be skipped if all exist and force is False.

            force: Step-specific force flag (None = use global)
        """
        self.input_files: List[Path] = [Path(p) for p in input_files] if input_files else []
        self.artifact_files: Optional[List[Path]] = [Path(p) for p in artifact_files] if artifact_files else None
        self._force = force

        # Whether artifacts were created (True if created, False if skipped, None if not run yet)
        self._artifacts_generated: Optional[bool] = None

    """ Logging """

    @property
    def name(self) -> str:
        """Step name (class name by default)."""
        return self.__class__.__name__

    """ Artifact files """

    @classmethod
    def set_global_force(cls, force: bool) -> None:
        """Set global force flag for all steps."""
        cls.global_force = force

    @property
    def force(self) -> bool:
        """Effective force flag (step-specific overrides global)."""
        if self._force is not None:
            return self._force
        return Step.global_force

    def _artifact_labels(self) -> list[str]:
        """Human-readable labels for artifacts (override for custom logging)."""
        return [str(p.name) for p in self.artifact_files] if self.artifact_files else []

    def _generate_artifacts(self) -> None:
        """Create artifact files (external tools). Subclasses should override."""
        pass

    def missing_artifact_files(self) -> Optional[list[Path]]:
        """Check if all artifact files exist. None if no artifacts tracked."""
        if self.artifact_files is None:
            return None
        return [p for p in self.artifact_files if not p.exists()]

    def has_artifact_files(self) -> bool:
        """True if artifact_files defined and all exist."""
        missing = self.missing_artifact_files()
        return missing is not None and len(missing) == 0

    """ Running """

    def missing_input_files(self) -> list[Path]:
        """Check if all inputs exist."""
        return [p for p in self.input_files if not p.exists()]

    def _terminate_run(self, reason: str, details: Optional[dict] = None) -> None:
        """Raise PrematureTerminationError with step name.

        Args:
            reason: Human-readable explanation of why termination occurred
            details: Optional dict with additional context (metrics, thresholds, etc.)

        Raises:
            PrematureTerminationError: Always raises with step name prepended
        """
        raise PrematureTerminationError(
            reason=reason,
            step=self.name,
            details=details or {}
        )

    def _should_skip_artifacts(self) -> bool:
        """Check if we should skip artifact generation."""
        return not self.force and self.has_artifact_files()

    def _get_status_message(self) -> str:
        """Determine what work will be done and return status message."""
        if self.artifact_files:
            return "skip artifacts" if self._should_skip_artifacts() else "create artifacts"
        else:
            return "no artifacts to generate"

    def _do_work(self):
        """Execute the work. Returns None for base Step."""
        self._generate_artifacts_if_needed()
        return None

    def _generate_artifacts_if_needed(self) -> None:
        """Generate artifacts if needed."""
        self._artifacts_generated = False

        if self.artifact_files:
            # Fast path: skip artifact generation if already present
            if self._should_skip_artifacts():
                return
            
            # Acquire lock if needed (None = no-op) and re-check (TOCTOU)
            lock_target = self._get_lock_target()
            with locked_resource(lock_target, self.name, timeout=self.STEP_LOCK_TIMEOUT):
                if self._should_skip_artifacts():
                    return
                else:
                    if self.force:
                        self._clean_artifacts()
                    self._run_generate_artifacts()

    def _run_generate_artifacts(self):
        """Run the creation of artifacts."""
        if self.should_profile:
            lp = self._create_profiler()
            lp.add_function(self._generate_artifacts)
            lp.runcall(self._generate_artifacts)
            self._save_profiler_stats(lp, suffix="artifacts")
        else:
            self._generate_artifacts()
        self._artifacts_generated = True

    def run(self):
        """Execute step: generate artifacts with caching and optional locking."""
        # Check inputs exist
        if missing_input := self.missing_input_files():
            raise FileNotFoundError(f"{self.name}: missing inputs: {missing_input}")

        if self._artifacts_generated is not None:
            raise RuntimeError(f"{self.name}: already ran; cannot run again")

        # Determine and log what we're about to do
        status_msg = self._get_status_message()
        self._log_run_status(status_msg)

        # Do the work (timed)
        with timer(log=False) as t:
            output = self._do_work()

        logger.info("-" * 107 + f" [cyan]{t.elapsed:.1f} sec[/cyan] --------\n")
        return output

    def _log_run_status(self, log_msg: str) -> None:
        """Standardized run status logging."""
        step_name = self.name
        step_name_color = f"[cyan]{step_name}[/cyan]"
        log_msg_color = f"[cyan]{log_msg}[/cyan]"

        labels = self._artifact_labels()
        output_str = ", ".join(labels) if labels else "none"

        remaining = 100 - len(step_name) - len(log_msg) - 2
        formatted = f"{step_name_color} " + "=" * remaining + f" {log_msg_color} ========"
        logger.info(formatted, timestamp=True)
        logger.info(f"Output: {output_str}")

    def _get_lock_target(self) -> Optional[Path]:
        """Path used for step lock; override to customize lock scope.
        
        Returns:
            Path to lock file, or None if no locking is needed (isolate-specific files).
            
        Override this method to:
        - Return a shared directory path for steps writing to shared resources
        - Return None for steps that only write to isolate-specific directories
        """
        # By default, no lock needed (isolate-specific files)
        return None

    def _clean_artifacts(self) -> None:
        """Remove existing artifacts before regeneration."""
        if self.artifact_files is None:
            return
        for p in self.artifact_files:
            remove_file_or_dir(p)

    """ Profiling """

    @classmethod
    def set_global_profile(cls, profile: bool) -> None:
        """Set global profile flag for all steps."""
        cls.should_profile = profile

    def _create_profiler(self) -> LineProfiler:
        """Create and configure LineProfiler with custom functions."""
        try:
            from line_profiler import LineProfiler
        except ImportError:
            raise ImportError(
                "line_profiler is required for profiling. "
                "Install it with: pip install 'amplifinder[dev]' or pip install line_profiler"
            )

        lp = LineProfiler()
        for func in self._get_profiler_functions():
            lp.add_function(func)
        return lp

    def _get_profiler_functions(self) -> list:
        """Override in subclass to specify functions to profile.

        Returns:
            List of functions to add to the line profiler.
        """
        return []

    def _get_profiler_stats_path(self, suffix: str = "") -> Optional[Path]:
        """Get path for saving profiler stats. Override in subclass to customize."""
        if self.artifact_files:
            suffix_str = f"_{suffix}" if suffix else ""
            return self.artifact_files[0].parent / f"line_profiler_stats{suffix_str}.lprof"
        return None

    def _save_profiler_stats(self, lp, suffix: str = "") -> None:
        """Save and optionally print line profiler statistics.

        Args:
            lp: LineProfiler instance with collected stats.
            suffix: Optional suffix for stats filename (e.g., "artifacts", "calculation").
        """
        from io import StringIO

        # Save stats to file
        stats_file = self._get_profiler_stats_path(suffix)
        if stats_file:
            lp.dump_stats(str(stats_file))
            logger.info(f"Line profiler stats saved to {stats_file}")

        # Print stats if verbose
        if logger.verbose:
            output_stream = StringIO()
            lp.print_stats(stream=output_stream, stripzeros=True)
            stats_text = output_stream.getvalue()
            if stats_text.strip():
                logger.info(f"Line profiler stats ({suffix or 'all'}):")
                for line in stats_text.strip().split('\n'):
                    logger.info(f"  {line}")


class OutputStep(Step, Generic[T]):
    """Step that creates artifacts AND returns output data to next step.

    Output files (CSVs, etc.) are treated as artifacts - included in artifact_files.
    Overrides _generate_artifacts_if_needed to also compute and return output.

    Extends Step with:
    - _calculate_output() -> T: Compute output from artifacts
    - _generate_artifacts_if_needed() -> T: Returns calculated output (not None)
    - _save_output(output): Called after artifacts are generated to save output files

    Use this when step needs to:
    - Generate artifacts (BAMs, CSVs, etc.)
    - Return data to next pipeline step

    Note: run() is inherited from Step and automatically returns the result of
    _generate_artifacts_if_needed(), which is T for OutputStep.
    """

    def __init__(
        self, *,
        input_files: Optional[List[Path]] = None,
        artifact_files: Optional[List[Path]] = None,
        output_files: Optional[List[Path]] = None,
        force: Optional[bool] = None,
    ):
        """Initialize output step.

        Args:
            input_files: Required input files/dirs (must exist)
            artifact_files: Files/dirs created by external tools (BAMs, etc.)
            output_files: Output files to create (CSVs, etc.) - added to artifacts
            force: Step-specific force flag (None = use global)
        """
        # Merge output_files into artifact_files
        if output_files is not None:
            artifact_files = (artifact_files or []) + output_files

        super().__init__(
            input_files=input_files,
            artifact_files=artifact_files,
            force=force,
        )

        self.output_files = output_files
        self._output_calculated: Optional[bool] = None

    @abstractmethod
    def _calculate_output(self) -> T:
        """Calculate output using existing artifacts and step attributes."""
        pass

    def report_output_message(self, output: T) -> Optional[str]:
        """Message to report after producing/loading output. Override per step."""
        return None

    def _save_output(self, output: T) -> None:
        """Save output to artifact files (e.g., CSV). Override in subclass."""
        pass

    def _get_status_message(self) -> str:
        """Determine what work will be done and return status message."""
        artifact_status = super()._get_status_message()
        return f"{artifact_status}, calc output"

    def _do_work(self) -> T:
        """Execute the work and return output."""
        # Generate artifacts
        super()._generate_artifacts_if_needed()

        # Calculate output
        if self.should_profile:
            lp = self._create_profiler()
            lp.add_function(self._calculate_output)
            output = lp.runcall(self._calculate_output)
            self._save_profiler_stats(lp, suffix="calculation")
        else:
            output = self._calculate_output()
        self._output_calculated = True

        # Save output files if artifacts were just generated (not cached)
        if self._artifacts_generated:
            self._save_output(output)

        # Report
        msg = self.report_output_message(output)
        if msg:
            logger.log_always(msg)

        return output


R = TypeVar("R", bound=Record)


class RecordTypedDfStep(OutputStep[RecordTypedDf[R]], Generic[R]):
    """Base class for steps that output RecordTypedDf to CSV (view only).

    Automatically handles:
    - Output file path from io_naming.default_path()
    - CSV save via RecordTypedDf

    Subclasses should:
    - Set class var `record_cls` (or it will be auto-deduced from typing)
    - Override `_calculate_output()` to return RecordTypedDf[R]
    """

    record_cls: Optional[Type[R]] = None  # Can be set explicitly or auto-deduced

    def __init__(
        self,
        *,
        output_dir: Optional[Path | str] = None,
        output_file: Optional[Path | str] = None,
        input_files: Optional[List[Path]] = None,
        artifact_files: Optional[List[Path]] = None,
        force: Optional[bool] = None,
    ):
        """Initialize step.

        Args:
            output_dir: Directory for output CSV file (str or Path; uses default filename from io_naming)
            output_file: Full path to output CSV file (str or Path; overrides output_dir)
            input_files: Required input files/dirs
            force: Step-specific force flag
        """
        if output_dir is not None and output_file is not None:
            raise ValueError("Cannot specify both output_dir and output_file")

        if output_file is None and output_dir is not None:
            # Use output_dir with default filename
            record_type = self._get_record_cls()
            output_file = Path(output_dir) / default_filename(record_type)

        super().__init__(
            input_files=input_files,
            artifact_files=artifact_files,
            output_files=[output_file] if output_file is not None else None,
            force=force,
        )

    @property
    def output_file(self) -> Optional[Path]:
        """Output file."""
        return self.output_files[0] if self.output_files else None

    @property
    def output_dir(self) -> Optional[Path]:
        """Output directory."""
        return self.output_files[0].parent if self.output_files else None

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

    def report_output_message(self, output: RecordTypedDf[R]) -> Optional[str]:
        """Uniform record count logging for RecordTypedDf steps."""
        record_cls = self._get_record_cls()
        return f"Created {len(output)} {record_cls.NAME}."

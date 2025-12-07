"""Pipeline step base class with caching logic."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import Generic, List, Optional, TypeVar
import shutil

from amplifinder.logger import info

T = TypeVar("T")


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
        """Execute step with caching logic."""
        # Check inputs exist
        if missing_input := self.missing_input_files():
            raise FileNotFoundError(f"{self.name}: missing inputs: {missing_input}")

        # Check if can skip (only if saving to files)
        if not self.force and self.has_output_files():
            info(f"{self.name}: skipped (outputs exist)")
            return self.load_outputs()

        # Clean partial outputs
        self._clean_outputs()

        # Run
        info(f"{self.name}: running...")
        output = self._calculate_output()

        # Save output (only if output_files defined)
        if self.output_files is not None:
            self._save_output_and_verify(output)

        return output

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

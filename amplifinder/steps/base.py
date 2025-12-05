"""Pipeline step base class with caching logic."""

from abc import ABC, abstractmethod
from pathlib import Path
from typing import List, Optional
import shutil

from amplifinder.logger import info, warning


class Step(ABC):
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
        inputs: Optional[List[Path]] = None,
        outputs: Optional[List[Path]] = None,
        force: Optional[bool] = None,
    ):
        """Initialize step.
        
        Args:
            inputs: Required input files/dirs (must exist)
            outputs: Output files/dirs (produced by this step)
            force: Step-specific force flag (None = use global)
        """
        self.inputs = [Path(p) for p in (inputs or [])]
        self.outputs = [Path(p) for p in (outputs or [])]
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
    def _run(self) -> None:
        """Execute the step logic. Override in subclass."""
        pass
    
    def read_outputs(self):
        """Load outputs from files. Override in subclass to return typed data."""
        return None

    def run_and_read_outputs(self):
        """Run step and return typed outputs."""
        self.run()
        return self.read_outputs()
    
    def has_input(self) -> bool:
        """Check if all inputs exist."""
        return all(p.exists() for p in self.inputs)
    
    def has_output(self) -> bool:
        """Check if all outputs exist."""
        return all(p.exists() for p in self.outputs)
    
    def run(self) -> bool:
        """Execute step with caching logic.
        
        Returns:
            True if step ran, False if skipped
        """
        # Check inputs exist
        if not self.has_input():
            missing = [p for p in self.inputs if not p.exists()]
            raise FileNotFoundError(f"{self.name}: missing inputs: {missing}")
        
        # Check if can skip
        if not self.force and self.has_output():
            info(f"{self.name}: skipped (outputs exist)")
            return False
        
        # Clean partial outputs
        self._clean_outputs()
        
        # Run
        info(f"{self.name}: running...")
        self._run()
        
        # Verify outputs created
        missing_out = [p for p in self.outputs if not p.exists()]
        if missing_out:
            warning(f"{self.name}: expected outputs not created: {missing_out}")
        
        return True
    
    def _clean_outputs(self) -> None:
        """Remove existing outputs before re-run."""
        for p in self.outputs:
            if p.exists():
                if p.is_dir():
                    shutil.rmtree(p)
                else:
                    p.unlink()
    
    @classmethod
    def set_global_force(cls, force: bool) -> None:
        """Set global force flag for all steps."""
        cls.global_force = force


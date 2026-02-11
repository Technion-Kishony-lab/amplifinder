"""Exception types for AmpliFinder."""


CONFIG_HELP_MSG = (
    "\n\nTo configure paths and settings, copy the bundled amplifinder.yaml to "
    "~/.amplifinder/amplifinder.yaml and edit it.\n"
    "See docs/configuration.md for details."
)


class PrematureTerminationError(Exception):
    """Raised when pipeline should terminate early due to data quality issues.

    This exception signals intentional early termination (e.g., breseq output
    suggests wrong reference) rather than an actual error. The pipeline catches
    this and handles it according to configuration.

    Attributes:
        reason: Human-readable explanation of why termination occurred
        step: Optional name of the step that triggered termination
        details: Optional dict with additional context (e.g., metrics, thresholds)
    """

    def __init__(self, reason: str, step: str = None, details: dict = None):
        """Initialize premature termination error.

        Args:
            reason: Human-readable explanation
            step: Name of the step that triggered termination
            details: Additional context (metrics, thresholds, etc.)
        """
        self.reason = reason
        self.step = step
        self.details = details or {}

        super().__init__(self._format_message())

    def _format_message(self) -> str:
        """Format message for exception display.

        Returns:
            Single-line message with step prefix if present
        """
        if self.step:
            return f"[{self.step}] {self.reason}"
        return self.reason

    def get_detailed_message(self) -> str:
        """Get full message with step and details for logging/status files.

        Returns:
            Multi-line string with reason, step, and all details
        """
        msg = self.reason
        if self.step:
            msg += f"\nstep: {self.step}"
        if self.details:
            for key, val in self.details.items():
                msg += f"\n{key}: {val}"
        return msg


class ToolNotFoundError(FileNotFoundError):
    """Tool executable not found with helpful message about amplifinder.yaml."""
    
    def __init__(self, message: str, include_help: bool = True):
        if include_help:
            message = f"{message}{CONFIG_HELP_MSG}"
        super().__init__(message)

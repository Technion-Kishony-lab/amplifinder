"""Simple custom logger for AmpliFinder."""

from collections import defaultdict
import sys
import re

from datetime import datetime
from pathlib import Path
from typing import Optional

from amplifinder.env import DEBUG
from amplifinder.utils.file_utils import ensure_parent_dir
from rich.console import Console


# ANSI color codes (fallback)
class Colors:
    """ANSI color codes for terminal output."""
    RESET = "\033[0m"

    # Level colors
    DEBUG = "\033[36m"      # Cyan
    INFO = "\033[32m"       # Green
    WARNING = "\033[33m"    # Yellow
    ERROR = "\033[31m"      # Red

    # Component colors
    TIMESTAMP = "\033[90m"  # Gray

    # Highlight colors
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"
    WHITE = "\033[37m"
    BOLD = "\033[1m"


class SimpleLogger:
    """Simple logger that writes to screen and/or file."""

    LEVEL_COLORS = {
        "DEBUG": "cyan",
        "INFO": "green",
        "WARNING": "yellow",
        "ERROR": "red",
    }

    CATEGORIES_TO_COUNTS: dict[str, int] = defaultdict(int)
    
    # Track console-printed messages for deduplication (batch mode)
    _console_printed: set[str] = set()

    def __init__(self, log_file: Optional[Path] = None, warnings_file: Optional[Path] = None, debug_file: Optional[Path] = None, use_colors: bool = True, verbose: bool = False):
        """Initialize logger.

        Args:
            log_file: Path to log file (None for screen only)
            warnings_file: Path to warnings-only file (None to skip)
            debug_file: Path to debug-only file (None to skip)
            use_colors: Use rich colors and auto-highlighting for screen output
            verbose: If True, show INFO/DEBUG messages; if False, only WARNING/ERROR
        """
        self.log_file = Path(log_file) if log_file else None
        self.warnings_file = Path(warnings_file) if warnings_file else None
        self.debug_file = Path(debug_file) if debug_file else None
        self.use_colors = use_colors
        self.verbose = verbose
        self.console = Console(file=sys.stdout, highlight=True) if use_colors else None

    def log(
        self,
        msg: str,
        *,
        level: str = "INFO",
        color: Optional[str] = None,
        timestamp: bool = False,
        to_screen: bool = True,
        to_file: bool = True,
        end: str = "\n",
        force_screen: bool = False,
        console_once: Optional[str] = None,
    ) -> None:
        """Log message to screen and/or file.

        Args:
            msg: Message to log. Can contain:
                 - Plain text (auto-highlighted by rich if use_colors=True)
                 - Rich markup like [cyan]text[/cyan] for inline colors
                 - Use colorize() helper for inline coloring
            level: Log level (DEBUG/INFO/WARNING/ERROR)
            color: Override entire message color (e.g., 'cyan', 'magenta', 'bold cyan')
                   If None, uses level-based color for WARNING/ERROR, auto-highlight for others
            timestamp: Include timestamp
            to_screen: Print to screen (subject to verbose filtering unless force_screen=True)
            to_file: Write to log file
            end: Line ending
            force_screen: If True, always show on screen regardless of verbose mode
            console_once: Dedupe key - if provided, prints to console only once per key.
                          Always writes to file regardless. Useful for shared resource warnings
                          in batch mode (e.g., console_once="genome_U00096")

        Examples:
            # Auto-highlighting (default)
            logger.info("Found 10 records")  # '10' auto-highlighted

            # Inline colors with colorize()
            logger.info(f"Status: {colorize('OK', 'green')}")

            # Inline colors with rich markup
            logger.info("Status: [green]OK[/green]")

            # Override entire message color
            logger.info("Important message", color="bold magenta")
            
            # Deduplicated console output (batch mode)
            logger.warning("Genome issue", console_once="genome_U00096")
        """
        level = level.upper()

        # Filter INFO/DEBUG messages if not verbose (always show WARNING/ERROR)
        # unless force_screen is True
        if not force_screen and not self.verbose and level in ("INFO", "DEBUG"):
            to_screen = False
        
        # Console deduplication: skip console if already printed
        if console_once and to_screen:
            if console_once in self._console_printed:
                to_screen = False
            else:
                self._console_printed.add(console_once)

        # Format for screen (with rich colors and auto-highlighting)
        if to_screen:
            if self.console:
                screen_msg = self._format_rich_message(msg, level, timestamp, color)
                self.console.print(screen_msg, end=end)
            else:
                screen_msg = self._format_message(msg, level, timestamp, use_colors=False)
                print(screen_msg, end=end, flush=True)

        # Format for file (no colors) - shared formatting for all files
        is_warning = level in ("WARNING", "ERROR") and self.warnings_file
        is_debug = level == "DEBUG" and self.debug_file
        if to_file and (self.log_file or is_warning or is_debug):
            file_msg = self._format_message(msg, level, timestamp, use_colors=False) + end
            
            # Write to main log file
            if self.log_file:
                with open(self.log_file, 'a') as f:
                    f.write(file_msg)
                    f.flush()
            
            # Write warnings/errors to separate warnings file
            if is_warning:
                with open(self.warnings_file, 'a') as f:
                    f.write(file_msg)
                    f.flush()
            
            # Write debug messages to separate debug file
            if is_debug:
                with open(self.debug_file, 'a') as f:
                    f.write(file_msg)
                    f.flush()

    def _format_rich_message(self, msg: str, level: str, timestamp: bool, color: Optional[str]) -> str:
        """Format message with rich colors and auto-highlighting.

        Args:
            msg: Message text (may contain rich markup like [cyan]text[/cyan])
            level: Log level
            timestamp: Include timestamp
            color: Override message color (None = use default behavior)

        Returns:
            String with rich markup (console.print will apply auto-highlighting)
        """
        parts = []

        if timestamp:
            time_str = datetime.now().strftime("%H:%M:%S")
            parts.append(f"[dim]{time_str}[/dim]")

            # Add level with color
            level_color = self.LEVEL_COLORS.get(level, "white")
            parts.append(f"[{level_color}]{level}[/{level_color}]:")

        # Priority: explicit color > level-based color > auto-highlighting
        if color:
            # User specified explicit color - wrap entire message
            parts.append(f"[{color}]{msg}[/{color}]")
        elif level in ("WARNING", "ERROR"):
            # Warnings/errors get level color unless they already have markup
            if '[' in msg and ']' in msg:
                # Message has markup, keep it
                parts.append(msg)
            else:
                # Plain message, apply level color
                level_color = self.LEVEL_COLORS.get(level, "white")
                parts.append(f"[{level_color}]{msg}[/{level_color}]")
        else:
            # Info/debug: pass message as-is for auto-highlighting + markup parsing
            parts.append(msg)

        return " ".join(parts)

    def _format_message(self, msg: str, level: str, timestamp: bool, use_colors: bool) -> str:
        """Format message with optional timestamp and colors (fallback)."""

        # Strip rich markup tags for plain text output
        msg = re.sub(r'\[/?[a-z\s]+\]', '', msg)

        parts = []

        if timestamp:
            time_str = datetime.now().strftime("%H:%M:%S")
            if use_colors and self.use_colors:
                parts.append(f"{Colors.TIMESTAMP}{time_str}{Colors.RESET}")
            else:
                parts.append(time_str)

            # Add level
            if use_colors and self.use_colors:
                level_color = self.LEVEL_COLORS.get(level, "")
                parts.append(f"{level_color}{level}{Colors.RESET}")
            else:
                parts.append(level)

            parts.append(msg)
            return " ".join(parts) + ":"
        else:
            return msg

    def debug(self, msg: str, **kwargs) -> None:
        """Log debug message."""
        self.log(msg, level="DEBUG", **kwargs)

    def info(self, msg: str, **kwargs) -> None:
        """Log info message."""
        self.log(msg, level="INFO", **kwargs)

    def warning(self, msg: str, **kwargs) -> None:
        """Log warning message."""
        self.log(msg, level="WARNING", **kwargs)

    def error(self, msg: str, **kwargs) -> None:
        """Log error message."""
        self.log(msg, level="ERROR", **kwargs)

    def set_verbose(self, verbose: bool) -> None:
        """Set verbose mode (True = show INFO/DEBUG, False = only WARNING/ERROR)."""
        self.verbose = verbose

    def reconfigure(self, log_file: Optional[Path] = None, warnings_file: Optional[Path] = None, debug_file: Optional[Path] = None, use_colors: bool = True, verbose: bool = False) -> None:
        """Reconfigure this logger instance (updates in place)."""
        if log_file:
            log_file = ensure_parent_dir(log_file)
        
        if warnings_file:
            warnings_file = ensure_parent_dir(warnings_file)
        
        if debug_file:
            debug_file = ensure_parent_dir(debug_file)

        self.log_file = Path(log_file) if log_file else None
        self.warnings_file = Path(warnings_file) if warnings_file else None
        self.debug_file = Path(debug_file) if debug_file else None
        self.use_colors = use_colors
        self.verbose = verbose
        self.console = Console(file=sys.stdout, highlight=True) if use_colors else None

    def log_always(self, msg: str, **kwargs) -> None:
        """Log message that always shows (ignores verbose mode)."""
        self.log(msg, force_screen=True, **kwargs)

    def print_progress(self, msg: str, end: str = "\n") -> None:
        """Print without timestamp (respects verbose mode, screen only)."""
        self.log(msg, timestamp=False, to_file=False, force_screen=False, end=end)

    def debug_message(self, message,
                      category: Optional[str] = None,
                      folder: Optional[Path] = None,
                      max_prints: Optional[int] = 1):
        """Debug message with category-based console deduplication.
        
        Args:
            message: Debug message to log
            category: Category for grouping/deduplication
            folder: Deprecated - category files now written to debug.txt
            max_prints: Max console prints per category (None = unlimited). 
                       Always writes to debug.txt regardless of limit.
        """
        if not DEBUG:
            return
        
        # Track counts per category for console limiting
        if category is None:
            count = 0
        else:
            count = self.CATEGORIES_TO_COUNTS[category]
            self.CATEGORIES_TO_COUNTS[category] += 1
        
        # Determine if should print to console based on max_prints
        should_print_console = (max_prints is None or count < max_prints)
        
        # Format message with category header
        full_msg = f"DEBUG[{category}]\n{message}" if category else message
        
        # Use standard debug() - always writes to debug.txt, conditional console
        self.debug(full_msg, to_screen=should_print_console)
        
        # Legacy folder support (deprecated, debug.txt is preferred)
        if folder and category:
            filepath = folder / (category + ".txt")
            with open(filepath, 'a') as f:
                f.write(message + "\n")


# Global logger instance
_logger: Optional[SimpleLogger] = None


def setup_logger(log_path: Optional[Path] = None, warnings_path: Optional[Path] = None, debug_path: Optional[Path] = None, use_colors: bool = True, verbose: bool = False) -> SimpleLogger:
    """Set up global logger (reconfigures existing instance).

    Args:
        log_path: Path to log file (None for screen only)
        warnings_path: Path to warnings-only file (None to skip)
        debug_path: Path to debug-only file (None to skip)
        use_colors: Use ANSI colors for screen output
        verbose: If True, show INFO/DEBUG messages; if False, only WARNING/ERROR

    Returns:
        Configured logger instance
    """
    logger = get_logger()
    logger.reconfigure(log_file=log_path, warnings_file=warnings_path, debug_file=debug_path, use_colors=use_colors, verbose=verbose)
    return logger


def get_logger() -> SimpleLogger:
    """Get the global logger, creating default if needed.

    Returns:
        Logger instance
    """
    global _logger
    if _logger is None:
        _logger = SimpleLogger()
    return _logger


# =============================================================================
# Color helpers for manual highlighting
# =============================================================================

def colorize(text: str, color: str) -> str:
    """Wrap text in rich markup for inline colored output.

    Use this for coloring specific parts of a message while keeping auto-highlighting
    for the rest. For coloring the entire message, use the `color` parameter instead.

    Args:
        text: Text to colorize
        color: Rich color/style name (cyan, magenta, yellow, red, green, blue,
               bold, dim, italic, or combinations like "bold cyan")

    Returns:
        Rich markup string (e.g., "[cyan]text[/cyan]")

    Examples:
        # Inline colors (specific text)
        logger.info(f"Found {colorize('10', 'cyan')} records")
        logger.info(f"Status: {colorize('OK', 'green')}")
        logger.info(f"File: {colorize('data.csv', 'bold magenta')}")

        # Entire message color (use color parameter instead)
        logger.info("Important message", color="bold cyan")
    """
    return f"[{color}]{text}[/{color}]"


def c(text: str, color: str) -> str:
    """Shorthand for colorize() - wrap text in rich markup.

    Args:
        text: Text to colorize
        color: Rich color/style name

    Returns:
        Rich markup string

    Examples:
        logger.info(f"Found {c('10', 'cyan')} records in {c('data.csv', 'magenta')}")
    """
    return colorize(text, color)


# =============================================================================
# Logger singleton instance (preferred API)
# =============================================================================

# Export the global logger singleton for convenient access
# Usage:
#   from amplifinder import logger
#   from amplifinder.logger import colorize, c  # for inline colors
#
#   # Auto-highlighting (default)
#   logger.info("Found 10 records")  # '10' auto-highlighted
#
#   # Inline colors
#   logger.info(f"Status: {c('OK', 'green')}")
#
#   # Entire message color
#   logger.info("Important", color="bold cyan")
logger = get_logger()

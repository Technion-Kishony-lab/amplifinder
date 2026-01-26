"""Simple custom logger for AmpliFinder."""

import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

from rich.console import Console
from rich.text import Text


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
    
    def __init__(self, log_file: Optional[Path] = None, use_colors: bool = True):
        """Initialize logger.
        
        Args:
            log_file: Path to log file (None for screen only)
            use_colors: Use rich colors and auto-highlighting for screen output
        """
        self.log_file = Path(log_file) if log_file else None
        self.use_colors = use_colors
        self.console = Console(file=sys.stdout, highlight=True) if use_colors else None
    
    def log(
        self,
        msg: str,
        *,
        level: str = "INFO",
        timestamp: bool = True,
        to_screen: bool = True,
        to_file: bool = True,
        end: str = "\n",
    ) -> None:
        """Log message to screen and/or file.
        
        Args:
            msg: Message to log (numbers/strings auto-highlighted if use_colors=True)
            level: Log level (DEBUG/INFO/WARNING/ERROR)
            timestamp: Include timestamp
            to_screen: Print to screen
            to_file: Write to log file
            end: Line ending
        """
        level = level.upper()
        
        # Format for screen (with rich colors and auto-highlighting)
        if to_screen:
            if self.console:
                screen_msg = self._format_rich_message(msg, level, timestamp)
                self.console.print(screen_msg, end=end)
            else:
                screen_msg = self._format_message(msg, level, timestamp, use_colors=False)
                print(screen_msg, end=end, flush=True)
        
        # Format for file (no colors)
        if to_file and self.log_file:
            file_msg = self._format_message(msg, level, timestamp, use_colors=False)
            with open(self.log_file, 'a') as f:
                f.write(file_msg + end)
    
    def _format_rich_message(self, msg: str, level: str, timestamp: bool) -> Text:
        """Format message with rich colors and auto-highlighting."""
        text = Text()
        
        if timestamp:
            time_str = datetime.now().strftime("%H:%M:%S")
            text.append(time_str, style="dim")
            text.append(" ")
            
            # Add level with color
            level_color = self.LEVEL_COLORS.get(level, "white")
            text.append(level, style=level_color)
            text.append(": ")
        
        # For warnings/errors, color entire message
        # For info/debug, apply auto-highlighting
        if level in ("WARNING", "ERROR"):
            text.append(msg, style=self.LEVEL_COLORS.get(level, "white"))
        else:
            # Add message (rich will auto-highlight numbers, strings, etc.)
            text.append(msg)
        
        return text
    
    def _format_message(self, msg: str, level: str, timestamp: bool, use_colors: bool) -> str:
        """Format message with optional timestamp and colors (fallback)."""
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


# Global logger instance
_logger: Optional[SimpleLogger] = None


def setup_logger(log_path: Optional[Path] = None, use_colors: bool = True) -> SimpleLogger:
    """Set up global logger.
    
    Args:
        log_path: Path to log file (None for screen only)
        use_colors: Use ANSI colors for screen output
    
    Returns:
        Configured logger instance
    """
    global _logger
    
    if log_path:
        from amplifinder.utils.file_utils import ensure_parent_dir
        log_path = ensure_parent_dir(log_path)
    
    _logger = SimpleLogger(log_path, use_colors)
    return _logger


def get_logger() -> SimpleLogger:
    """Get the global logger, creating default if needed.
    
    Returns:
        Logger instance
    """
    global _logger
    if _logger is None:
        _logger = SimpleLogger()
    return _logger


# Convenience functions
def log(msg: str, level: str = "info", **kwargs) -> None:
    """Log a message at the specified level."""
    logger = get_logger()
    getattr(logger, level.lower(), logger.info)(msg, **kwargs)


def info(msg: str, **kwargs) -> None:
    """Log info message."""
    get_logger().info(msg, **kwargs)


def warning(msg: str, **kwargs) -> None:
    """Log warning message."""
    get_logger().warning(msg, **kwargs)


def error(msg: str, **kwargs) -> None:
    """Log error message."""
    get_logger().error(msg, **kwargs)


def debug(msg: str, **kwargs) -> None:
    """Log debug message."""
    get_logger().debug(msg, **kwargs)


# Color helpers for manual highlighting
def colorize(text: str, color: str) -> str:
    """Wrap text in rich markup for colored output.
    
    Args:
        text: Text to colorize
        color: Color name (blue, cyan, magenta, yellow, red, green, bold)
    
    Returns:
        Rich markup string (e.g., "[cyan]text[/cyan]")
    
    Examples:
        >>> self.log(f"Found {colorize('10', 'cyan')} records")
        >>> self.log(f"Status: {colorize('OK', 'green')}")
    """
    return f"[{color}]{text}[/{color}]"

"""Logging utilities for AmpliFinder."""

import logging
import sys
from pathlib import Path
from typing import Optional

# Module-level logger
_logger: Optional[logging.Logger] = None


class HeaderFormatter(logging.Formatter):
    """Ensures optional header fields exist."""

    def format(self, record: logging.LogRecord) -> str:
        if not hasattr(record, "header"):
            record.header = ""
        if not hasattr(record, "header_color"):
            record.header_color = ""
        return super().format(record)


class ColorFormatter(HeaderFormatter):
    """ANSI color formatter for console logs."""

    COLORS = {
        logging.DEBUG: "\033[36m",  # cyan
        logging.INFO: "\033[32m",  # green
        logging.WARNING: "\033[33m",  # yellow
        logging.ERROR: "\033[31m",  # red
        logging.CRITICAL: "\033[35m",  # magenta
    }
    RESET = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        msg = super().format(record)
        color = self.COLORS.get(record.levelno)
        return f"{color}{msg}{self.RESET}" if color else msg


def setup_logger(
    log_path: Optional[Path] = None,
    level: int = logging.INFO,
    name: str = "amplifinder",
) -> logging.Logger:
    """Set up logger with console and optional file output.

    Args:
        log_path: Path to log file (None for console only)
        level: Logging level (default: INFO)
        name: Logger name

    Returns:
        Configured logger instance
    """
    global _logger

    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Clear existing handlers
    logger.handlers.clear()

    fmt = "%(asctime)s %(levelname)s: %(header_color)s%(message)s"
    plain_fmt = "%(asctime)s %(levelname)s: %(header)s%(message)s"
    date_fmt = "%H:%M:%S"
    color_formatter = ColorFormatter(fmt, datefmt=date_fmt)
    plain_formatter = HeaderFormatter(plain_fmt, datefmt=date_fmt)

    # Console handler (stdout)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(color_formatter)
    logger.addHandler(console_handler)

    # File handler (if path provided)
    if log_path is not None:
        from amplifinder.utils.file_utils import ensure_parent_dir
        log_path = ensure_parent_dir(log_path)

        file_handler = logging.FileHandler(log_path, mode="a")
        file_handler.setLevel(level)
        file_handler.setFormatter(plain_formatter)
        logger.addHandler(file_handler)

    _logger = logger
    return logger


def get_logger() -> logging.Logger:
    """Get the module logger, creating a default one if needed.

    Returns:
        Logger instance
    """
    global _logger
    if _logger is None:
        _logger = setup_logger()
    return _logger


def log(msg: str, level: str = "info", **kwargs) -> None:
    """Log a message at the specified level.

    Args:
        msg: Message to log
        level: Log level (debug, info, warning, error)
    """
    logger = get_logger()
    level_map = {
        "debug": logger.debug,
        "info": logger.info,
        "warning": logger.warning,
        "error": logger.error,
    }
    log_func = level_map.get(level.lower(), logger.info)
    log_func(msg, **kwargs)


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

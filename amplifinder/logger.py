"""Logging utilities for AmpliFinder."""

import logging
import sys
from pathlib import Path
from typing import Optional

# Module-level logger
_logger: Optional[logging.Logger] = None


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
    
    # Format: TYPE: message | timestamp
    formatter = logging.Formatter(
        "%(levelname)s: %(message)s | %(asctime)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    
    # Console handler (stdout)
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    # File handler (if path provided)
    if log_path is not None:
        log_path = Path(log_path)
        log_path.parent.mkdir(parents=True, exist_ok=True)
        
        file_handler = logging.FileHandler(log_path, mode="a")
        file_handler.setLevel(level)
        file_handler.setFormatter(formatter)
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


def log(msg: str, level: str = "info") -> None:
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
    log_func(msg)


def info(msg: str) -> None:
    """Log info message."""
    get_logger().info(msg)


def warning(msg: str) -> None:
    """Log warning message."""
    get_logger().warning(msg)


def error(msg: str) -> None:
    """Log error message."""
    get_logger().error(msg)


def debug(msg: str) -> None:
    """Log debug message."""
    get_logger().debug(msg)


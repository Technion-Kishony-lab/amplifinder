"""Lightweight timing helpers."""

from contextlib import contextmanager
from time import perf_counter
from typing import Dict, Optional

from amplifinder.logger import logger

_timers: Dict[str, float] = {}


class TimerResult:
    """Result object returned by timer context manager."""
    def __init__(self):
        self.elapsed: float = 0.0


@contextmanager
def timer(msg: Optional[str] = None, log: bool = True, extra: Optional[dict[str, str]] = None):
    """Context manager for timing code blocks."""
    result = TimerResult()
    start = perf_counter()
    try:
        yield result
    finally:
        elapsed = perf_counter() - start
        result.elapsed = elapsed
        if log and msg:
            logger.info(f"{msg} completed in {elapsed:.3f}s", extra=extra)


@contextmanager
def print_timer(start_msg: str, end_msg: Optional[str] = None, time_format: str = "{:.1f} sec",
                should_log: bool = True, seperate_prints: bool = True, use_log: bool = False, to_file: bool = False):
    """Context manager that prints start message, runs code, then prints time.

    Args:
        start_msg: Message to print before running code
        end_msg: Message to print after (default: time only)
        time_format: Format string for time (default: "{:.1f} sec")
        should_log: If False, don't print anything (default: True)
        seperated_prints: If True, print start message and end message separately (default: True)
        use_log: If True, use logger (respects verbose) instead of print (default: True)
        to_file: If False, don't write to log file (default: False)
    Example:
        with print_timer("Building index ... ", end_msg="\n"):
            build_index()
        # Output: "Building index ... \n12.3 sec\n"
    """
    end_msg = end_msg or ""
    if should_log and seperate_prints:
        if use_log:
            logger.info(start_msg, to_file=to_file)
        else:
            print(start_msg, end="", flush=True)
    result = TimerResult()
    start = perf_counter()
    try:
        yield result
    finally:
        elapsed = perf_counter() - start
        result.elapsed = elapsed
        if should_log:
            prefix = "" if seperate_prints else start_msg
            time_str = time_format.format(elapsed)
            msg = f"{prefix}{time_str}{end_msg}"
            if use_log:
                logger.info(msg, to_file=to_file)
            else:
                print(msg, flush=True, end="")

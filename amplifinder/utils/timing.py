"""Lightweight timing helpers."""

from contextlib import contextmanager
from time import perf_counter
from typing import Dict, Optional

from amplifinder.logger import info

_timers: Dict[str, float] = {}


class TimerResult:
    """Result object returned by timer context manager."""
    def __init__(self):
        self.elapsed: float = 0.0


def start_timer(timer: str = "default") -> None:
    _timers[timer] = perf_counter()


def end_timer(msg: Optional[str], timer: str = "default", log: bool = True, extra: Optional[dict[str, str]] = None) -> float:
    start = _timers.pop(timer, None)
    if start is None:
        raise ValueError(f"Timer {timer} not started")
    elapsed = perf_counter() - start
    if log and msg:
        info(f"{msg}: {elapsed:.3f}s", extra=extra)
    return elapsed


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
            info(f"{msg} completed in {elapsed:.3f}s", extra=extra)


@contextmanager
def print_timer(start_msg: str, end_msg: Optional[str] = None, time_format: str = "{:.1f} sec", should_log: bool = True):
    """Context manager that prints start message, runs code, then prints time.
    
    Args:
        start_msg: Message to print before running code (printed with end="", flush=True)
        end_msg: Message to print after (default: time only)
        time_format: Format string for time (default: "{:.1f} sec")
        should_log: If False, don't print anything (default: True)
    
    Example:
        with print_timer("Building index ... ", end_msg=", "):
            build_index()
        # Output: "Building index ... (12.3 sec), "
    """
    if should_log:
        print(start_msg, end="", flush=True)
    result = TimerResult()
    start = perf_counter()
    try:
        yield result
    finally:
        elapsed = perf_counter() - start
        result.elapsed = elapsed
        if should_log:
            time_str = time_format.format(elapsed)
            if end_msg is not None:
                print(f"({time_str}){end_msg}", end="", flush=True)
            else:
                print(f"({time_str})", end="", flush=True)

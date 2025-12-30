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
def print_timer(start_msg: str, end_msg: Optional[str] = None, time_format: str = "{:.1f} sec",
                should_log: bool = True, seperate_prints: bool = False, use_log: bool = False):
    """Context manager that prints start message, runs code, then prints time.
    
    Args:
        start_msg: Message to print before running code (printed with end="", flush=True)
        end_msg: Message to print after (default: time only)
        time_format: Format string for time (default: "{:.1f} sec")
        should_log: If False, don't print anything (default: True)
        seperated_prints: If True, print start message and end message separately
        use_log: If True, use logger.info instead of print
    Example:
        with print_timer("Building index ... ", end_msg=", "):
            build_index()
        # Output: "Building index ... (12.3 sec), "
    """
    end_msg = end_msg or ""
    if should_log and seperate_prints:
        if use_log:
            info(start_msg)
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
                info(msg)
            else:
                print(msg, flush=True)

"""Lightweight timing helpers."""

from contextlib import contextmanager
from time import perf_counter
from typing import Dict, Optional

from amplifinder.logger import info

_timers: Dict[str, float] = {}


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
    start = perf_counter()
    try:
        yield
    finally:
        elapsed = perf_counter() - start
        if log and msg:
            info(f"{msg} completed in {elapsed:.3f}s", extra=extra)

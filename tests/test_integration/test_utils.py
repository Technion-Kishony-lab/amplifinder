"""Utility functions for integration tests."""

from contextlib import contextmanager


def print_color(txt):
    """Print in magenta color for test output visibility."""
    MAGENTA = '\033[95m'
    RESET = '\033[0m'
    print(f"{MAGENTA}{txt}{RESET}", flush=True)


@contextmanager
def force_step():
    """Context manager to temporarily set Step.global_force = True."""
    from amplifinder.steps.base import Step
    Step.global_force = True
    try:
        yield
    finally:
        Step.global_force = False

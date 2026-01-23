"""Flag utilities."""

from contextlib import contextmanager


class MutableFlag:
    """A flag that can be temporarily set using a context manager."""

    def __init__(self, initial_value: bool = False):
        self.value = initial_value

    @contextmanager
    def temp_set(self, value: bool):
        """Temporarily set the flag value."""
        old_value = self.value
        self.value = value
        try:
            yield
        finally:
            self.value = old_value

    def __bool__(self):
        return self.value

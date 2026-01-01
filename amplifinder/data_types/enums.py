"""Enum definitions for AmpliFinder."""

from enum import Enum


class ReversibleIntEnum(int, Enum):
    """Base class for int enums with opposite() method."""

    def opposite(self):
        """Return the enum member with negated value."""
        return type(self)(-self.value)


class Side(ReversibleIntEnum):
    """Side of a TN element (left or right)."""
    LEFT = -1
    RIGHT = 1


class Orientation(ReversibleIntEnum):
    """Orientation relative to reference (forward, reverse, or both/mixed)."""
    FORWARD = 1
    REVERSE = -1
    BOTH = 0


class AverageMethod(str, Enum):
    """Method for calculating average coverage statistics."""
    MEDIAN = "median"
    MODE = "mode"
    MEAN = "mean"

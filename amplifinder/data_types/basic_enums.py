"""Enum definitions for AmpliFinder."""
from __future__ import annotations
import numpy as np
import operator

from dataclasses import dataclass
from enum import Enum
from typing import ClassVar, Optional, TypeAlias


class ReversibleIntEnum(int, Enum):
    """Base class for int enums with opposite() method."""

    def opposite(self):
        """Return the enum member with negated value."""
        return type(self)(-self.value)


class Terminal(ReversibleIntEnum):
    """Ends of a genetic element (start or end)."""
    START = -1
    END = 1


class Side(ReversibleIntEnum):
    """Side of a junction (left, middle, or right)."""
    LEFT = -1
    MIDDLE = 0
    RIGHT = 1


class Orientation(ReversibleIntEnum):
    """Orientation relative to reference (forward, reverse, or both/mixed)."""
    REVERSE = -1
    BOTH = 0   # TODO: This is not used yet
    FORWARD = 1


class AverageMethod(str, Enum):
    """Method for calculating average coverage statistics."""
    MEDIAN = "median"
    MODE = "mode"
    MEAN = "mean"

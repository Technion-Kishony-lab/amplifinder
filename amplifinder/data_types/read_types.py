"""Read types and junction read counts for AmpliFinder."""
from __future__ import annotations
import operator

from dataclasses import dataclass
from enum import Enum
from typing import Optional

from amplifinder.data_types.basic_enums import Side


class ReadType(str, Enum):
    """Type of a read at a junction."""
    LEFT_FAR = "left_far"
    LEFT = "left"
    LEFT_MARGINAL = "left_marginal"
    SPANNING = "spanning"
    PAIRED = "paired"
    RIGHT_MARGINAL = "right_marginal"
    RIGHT = "right"
    RIGHT_FAR = "right_far"

    def is_marginal(self) -> bool:
        """Return True if this is a marginal read type."""
        return self in [self.LEFT_MARGINAL, self.RIGHT_MARGINAL]

    def is_far(self) -> bool:
        """Return True if this is a far read type."""
        return self in [self.LEFT_FAR, self.RIGHT_FAR]

    def get_side(self) -> Side:
        """Return the side of the read type."""
        if self in [self.LEFT_FAR, self.LEFT, self.LEFT_MARGINAL]:
            return Side.LEFT
        elif self in [self.RIGHT_FAR, self.RIGHT, self.RIGHT_MARGINAL]:
            return Side.RIGHT
        return Side.MIDDLE


@dataclass
class JunctionReadCounts:
    """Read counts at a junction.
    Order: LEFT_FAR, LEFT, LEFT_MARGINAL, SPANNING, PAIRED, RIGHT_MARGINAL, RIGHT, RIGHT_FAR.
    """
    left_far: int = 0        # reads far on left side of junction
    left: int = 0            # reads on left side of junction
    left_marginal: int = 0   # reads partially overlapping from left
    spanning: int = 0        # reads spanning the junction
    paired: int = 0          # paired-end reads with different classifications
    right_marginal: int = 0  # reads partially overlapping from right
    right: int = 0           # reads on right side of junction
    right_far: int = 0       # reads far on right side of junction

    @classmethod
    def get_field(cls, read_type: ReadType) -> str:
        """Get the field name for a read type."""
        return read_type.value

    def __getitem__(self, read_type: ReadType) -> int:
        """Get count for a read type."""
        return getattr(self, self.get_field(read_type))

    def __setitem__(self, read_type: ReadType, value: int) -> None:
        """Set count for a read type."""
        setattr(self, self.get_field(read_type), value)

    @property
    def counts(self) -> dict[ReadType, int]:
        """Return counts as a dict mapping ReadType to count."""
        return {k: self[k] for k in ReadType}

    @property
    def total(self) -> int:
        """Return total number of reads."""
        return sum(self.counts.values())

    def increment(self, read_type: Optional[ReadType]) -> None:
        """Increment the count for the given read type."""
        self[read_type] += 1

    @classmethod
    def from_scalar(cls, scalar: int | float) -> JunctionReadCounts:
        """Create a JunctionReadCounts object from a scalar."""
        return JunctionReadCounts(**{cls.get_field(rt): scalar for rt in ReadType})

    def _apply_op(self, other: JunctionReadCounts | int | float, op) -> JunctionReadCounts:
        """Apply a binary operation to all fields."""

        if not isinstance(other, JunctionReadCounts):
            other = self.from_scalar(other)

        return JunctionReadCounts(**{
            self.get_field(rt): op(self[rt], other[rt]) for rt in ReadType
        })

    def max(self, other: JunctionReadCounts | int | float) -> JunctionReadCounts:
        """Return element-wise maximum of self and other."""
        return self._apply_op(other, max)

    # Type stubs for dynamically added operators
    def __add__(self, other: JunctionReadCounts | int | float) -> JunctionReadCounts: ...
    def __sub__(self, other: JunctionReadCounts | int | float) -> JunctionReadCounts: ...
    def __mul__(self, other: JunctionReadCounts | int | float) -> JunctionReadCounts: ...
    def __truediv__(self, other: JunctionReadCounts | int | float) -> JunctionReadCounts: ...
    def __floordiv__(self, other: JunctionReadCounts | int | float) -> JunctionReadCounts: ...


# Dynamically add arithmetic operators to JunctionReadCounts
for _op_name, _op_func in [
    ('__add__', operator.add),
    ('__sub__', operator.sub),
    ('__mul__', operator.mul),
    ('__truediv__', operator.truediv),
    ('__floordiv__', operator.floordiv),
]:
    setattr(
        JunctionReadCounts,
        _op_name,
        lambda self, other, op=_op_func: self._apply_op(other, op)
    )

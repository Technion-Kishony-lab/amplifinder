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


class BaseRawEvent(str, Enum):
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"
    LOCUS_JOINING = "locus-joining"

    def is_single_locus(self) -> bool:
        """Return True if this is a single locus event."""
        return self in [self.REFERENCE, self.TRANSPOSITION]


class RawEvent(str, Enum):
    """Structural classification based on junction pair relationships (Step 8)."""
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"

    # Locus-joining:
    UNFLANKED = "unflanked"
    HEMI_FLANKED_LEFT = "hemi-flanked left"
    HEMI_FLANKED_RIGHT = "hemi-flanked right"
    FLANKED = "flanked"
    MULTIPLE_SINGLE_LOCUS = "multiple single locus"  # TODO: This is not used yet
    UNRESOLVED = "unresolved"


class Element(str, Enum):
    """Elements (chromosome, cassette, IS)."""
    CHR = "chr"  # chromosome
    AMP = "amp"  # cassette
    TN = "tn"    # IS


class JunctionType(Enum):
    """The 7 synthetic junction types.

    Amplicon structure:

    ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~

    legend:
    ~~~    chromosome
    >>>    IS
    ====== cassette (amplicon)

    """
    CHR_TO_AMP_LEFT       = (1, '~~-==', Element.CHR, Element.AMP, -1)  # noqa: E221, E203
    CHR_TO_TN_LEFT        = (2, '~~->>', Element.CHR, Element.TN , -2)  # noqa: E221, E203
    AMP_RIGHT_TO_TN_LEFT  = (3, '==->>', Element.AMP, Element.TN , +3)  # noqa: E221, E203
    AMP_RIGHT_TO_AMP_LEFT = (4, '==-==', Element.AMP, Element.AMP,  0)  # noqa: E221, E203
    TN_RIGHT_TO_AMP_LEFT  = (5, '>>-==', Element.TN,  Element.AMP, -3)  # noqa: E221, E203
    TN_RIGHT_TO_CHR       = (6, '>>-~~', Element.TN,  Element.CHR, +2)  # noqa: E221, E203
    AMP_RIGHT_TO_CHR      = (7, '==-~~', Element.AMP, Element.CHR, +1)  # noqa: E221, E203

    @classmethod
    def sorted(cls) -> list[JunctionType]:
        """Return the JunctionType enum members sorted by value."""
        return sorted(cls, key=lambda x: x.num)

    @property
    def num(self) -> int:
        """Return the number of the junction. (1-7, matlab legacy)"""
        return self.value[0]

    @property
    def symbol(self) -> str:
        """Return the symbol of the junction."""
        return self.value[1]

    @property
    def elements(self) -> tuple[Element, Element]:
        """Return the elements of the junction."""
        return self.value[2], self.value[3]

    @property
    def order(self) -> int:
        """Return the number of the junction."""
        return self.value[4]

    @property
    def side(self) -> Side:
        """Return the side of the junction."""
        return Side(np.sign(self.order))


JcCall: TypeAlias = bool | None
"""Junction call: True if the junction is called, False if not, None if not enough evidence."""


class EventModifier(str, Enum):
    """Modifiers for classified events (Step 13)."""
    ANCESTRAL = "ancestral"
    DENOVO_LEFT = "de novo left"
    DENOVO_RIGHT = "de novo right"
    LOW_COVERAGE = "low coverage near junction"  # TODO: This is not used yet


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
        if read_type is None:
            return
        self[read_type] += 1

    @classmethod
    def from_scalar(cls, scalar: int | float) -> JunctionReadCounts:
        """Create a JunctionReadCounts object from a scalar."""
        return JunctionReadCounts(
            left=scalar, left_marginal=scalar, spanning=scalar, right_marginal=scalar, right=scalar
        )

    def _apply_op(self, other: JunctionReadCounts | int | float, op) -> JunctionReadCounts:
        """Apply a binary operation to all fields."""

        if not isinstance(other, JunctionReadCounts):
            other = self.from_scalar(other)

        return JunctionReadCounts(
            left=op(self.left, other.left),
            left_marginal=op(self.left_marginal, other.left_marginal),
            spanning=op(self.spanning, other.spanning),
            right_marginal=op(self.right_marginal, other.right_marginal),
            right=op(self.right, other.right)
        )

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

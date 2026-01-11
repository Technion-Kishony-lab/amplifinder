"""Enum definitions for AmpliFinder."""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


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


class JunctionType(int, Enum):
    """The 7 synthetic junction types.

    Amplicon structure:

    ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~

    ~~-==  ~~->>  ==->>  ==-==  >>-==  >>-~~  ==-~~
      1      2      3      4      5      6      7

    legend:
    ~~~    chromosome
    >>>    IS
    ====== cassette (amplicon)
    """
    CHR_TO_AMP_LEFT = 1          # (1) ~~-==  left reference (chromosome-cassette)
    CHR_TO_TN_LEFT = 2           # (2) ~~->>  left IS transposition (chromosome-IS)
    AMP_RIGHT_TO_TN_LEFT = 3     # (3) ==->>  left of mid IS (cassette-IS)
    AMP_RIGHT_TO_AMP_LEFT = 4    # (4) ==-==  lost IS (cassette-cassette, no IS)
    TN_RIGHT_TO_AMP_LEFT = 5     # (5) >>-==  right of mid IS (IS-cassette)
    TN_RIGHT_TO_CHR = 6          # (6) >>-~~  right IS transposition (IS-chromosome)
    AMP_RIGHT_TO_CHR = 7         # (7) ==-~~  right reference (cassette-chromosome)

    @classmethod
    def sorted(cls) -> list[JunctionType]:
        """Return the JunctionType enum members sorted by value."""
        return sorted(cls, key=lambda x: x.value)


class EventModifier(str, Enum):
    """Modifiers for classified events (Step 13)."""
    ANCESTRAL = "ancestral"
    DENOVO_LEFT = "de novo left"
    DENOVO_RIGHT = "de novo right"
    LOW_COVERAGE = "low coverage near junction"  # TODO: This is not used yet


@dataclass
class JunctionReadCounts:
    """Read counts at a junction."""
    left: int = 0      # reads on left side of junction
    right: int = 0     # reads on right side of junction
    spanning: int = 0  # reads spanning the junction
    undetermined: int = 0     # reads partially overlapping the junction

    def add_read(self, start: int, end: int, junction_point: int, min_overlap_len: int) -> str:
        # Determine read type and update counts
        if end <= junction_point:
            read_type = 'left'
            self.left += 1
        elif start > junction_point:
            read_type = 'right'
            self.right += 1
        elif start <= junction_point - min_overlap_len and end >= junction_point + min_overlap_len:
            read_type = 'spanning'
            self.spanning += 1
        else:
            read_type = 'undetermined'
            self.undetermined += 1
        return read_type

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


class BaseRawEvent(str, Enum):
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"
    LOCUS_JOINING = "locus-joining"


class RawEvent(str, Enum):
    """Structural classification based on junction pair relationships (Step 8)."""
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"

    # Locus-joining:
    UNFLANKED = "unflanked"
    HEMI_FLANKED_LEFT = "hemi-flanked left"
    HEMI_FLANKED_RIGHT = "hemi-flanked right"
    FLANKED = "flanked"
    MULTIPLE_SINGLE_LOCUS = "multiple single locus"
    UNRESOLVED = "unresolved"


class JunctionType(int, Enum):
    """The 7 synthetic junction types.

    Amplicon structure: ~~~>>>======>>>======>>>~~~
    (1) ~~==  left reference (chromosome-cassette)
    (2) ~~>>  left IS transposition (chromosome-IS)
    (3) ==>>  left of mid IS (cassette-IS)
    (4) ====  lost IS (cassette-cassette, no IS)
    (5) >>==  right of mid IS (IS-cassette)
    (6) >>~~  right IS transposition (IS-chromosome)
    (7) ==~~  right reference (cassette-chromosome)
    """
    LEFT_REF = 1
    LEFT_IS_TRANS = 2
    LEFT_MID_IS = 3
    LOST_IS = 4
    RIGHT_MID_IS = 5
    RIGHT_IS_TRANS = 6
    RIGHT_REF = 7


class EventModifier(str, Enum):
    """Modifiers for classified events (Step 13)."""
    ANCESTRAL = "ancestral"
    DE_NOVO = "de novo"
    LOW_COVERAGE = "low coverage near junction"


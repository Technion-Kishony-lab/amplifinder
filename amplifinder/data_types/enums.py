"""Enum definitions for AmpliFinder."""

from enum import Enum
from typing import NamedTuple


class ReversibleIntEnum(int, Enum):
    """Base class for int enums with opposite() method."""

    def opposite(self):
        """Return the enum member with negated value."""
        return type(self)(-self.value)


class Side(ReversibleIntEnum):
    """Side of a TN element (start or end)."""
    START = -1
    END = 1


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

    Amplicon structure: 
    
    ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~
    
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


class EventModifier(str, Enum):
    """Modifiers for classified events (Step 13)."""
    ANCESTRAL = "ancestral"
    DE_NOVO = "de novo"
    LOW_COVERAGE = "low coverage near junction"


class JunctionCoverage(NamedTuple):
    """Read coverage at a synthetic junction."""
    spanning: int  # reads crossing junction
    left: int      # reads ending at junction
    right: int     # reads starting at junction

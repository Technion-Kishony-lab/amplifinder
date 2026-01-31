"""Event classification enums for AmpliFinder."""
from __future__ import annotations

from enum import Enum
from typing import Optional


class BaseEvent(str, Enum):
    REFERENCE_TN = "reference-tn"
    TRANSPOSITION = "transposition"
    LOCUS_JOINING = "locus-joining"

    def is_single_locus(self) -> bool:
        """Return True if this is a single locus event."""
        return self in [self.REFERENCE_TN, self.TRANSPOSITION]


class Architecture(Enum):
    """Structural classification based on junction pair relationships (Step 8)."""
    # NOTE: For each event, the tuple contains:
    # - description
    # - (left_flanked, right_flanked)
    # - is_singleton

    # Single-locus events:
    REFERENCE_TN = "reference-tn", (None, None), None
    TRANSPOSITION = "transposition", (None, None), None
    TRANSPOSITION_SINGLETON = "transposition-singleton", (None, None), None

    # Locus-joining (multi-copy):
    UNFLANKED = "unflanked", (False, False), False
    HEMI_FLANKED_LEFT = "hemi-flanked left", (False, True), False
    HEMI_FLANKED_RIGHT = "hemi-flanked right", (True, False), False
    FLANKED = "flanked", (True, True), False
    
    # Locus-joining (singleton):
    HEMI_FLANKED_LEFT_SINGLETON = "hemi-flanked left singleton", (False, True), True
    HEMI_FLANKED_RIGHT_SINGLETON = "hemi-flanked right singleton", (True, False), True
    AMPLICON_ONLY = "amplicon-only", (False, False), True
    
    UNRESOLVED = "unresolved", (None, None), None

    @property
    def description(self) -> str:
        """Return the description of the event."""
        return self.value[0]

    @property
    def is_flanking(self) -> tuple[Optional[bool], Optional[bool]]:
        """Return True if this is a flanking event."""
        return self.value[1]

    @property
    def is_singleton(self) -> Optional[bool]:
        """Return True if this is a single locus event."""
        return self.value[2]


class EventDescriptor(str, Enum):
    """Descriptors for classified events."""
    ANCESTRAL = "ancestral"
    DENOVO_LEFT = "de novo left"
    DENOVO_RIGHT = "de novo right"
    LOST_LEFT = "lost left"
    LOST_RIGHT = "lost right"
    TO_SINGLETON = "to singleton"
    FROM_SINGLETON = "from singleton"
    AMP_AMP = "observed amplicon-amplicon (IS loss)"
    LOW_COVERAGE = "low coverage near junction"  # TODO: This is not used yet

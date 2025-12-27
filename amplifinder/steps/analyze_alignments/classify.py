"""Junction architecture classification."""

from typing import Optional, List

from amplifinder.data_types import RawEvent, EventModifier
from amplifinder.steps.analyze_alignments.parse_bam import JunctionReadCounts


def classify_architecture(jc_cov: List[JunctionReadCounts], min_jct_cov: int = 5) -> RawEvent:
    """Classify junction architecture from read coverage patterns.
    
    Based on MATLAB classify_candidates.m
    
    Junction patterns (1=covered, 0=not):
        Pattern         Name
        [0,1,1,0,1,1,0] flanked
        [0,1,1,0,1,0,1] hemi-flanked left
        [1,0,1,0,1,1,0] hemi-flanked right
        [1,0,1,0,1,0,1] unflanked
        [0,1,0,0,1,0,1] hemi-flanked left single
        [1,0,1,0,0,1,0] hemi-flanked right single
        [1,0,0,0,0,0,1] no IS (reference)
        [0,1,0,0,0,1,0] deletion
    
    Args:
        jc_cov: List of 7 JunctionReadCounts (one per junction type)
        min_jct_cov: Minimum spanning reads to consider junction "covered"
    
    Returns:
        RawEvent classification
    """
    # Build binary pattern from spanning coverage
    pattern = [1 if jc.spanning >= min_jct_cov else 0 for jc in jc_cov]
    
    # Pattern matching
    patterns = {
        (0, 1, 1, 0, 1, 1, 0): RawEvent.FLANKED,
        (0, 1, 1, 0, 1, 0, 1): RawEvent.HEMI_FLANKED_LEFT,
        (1, 0, 1, 0, 1, 1, 0): RawEvent.HEMI_FLANKED_RIGHT,
        (1, 0, 1, 0, 1, 0, 1): RawEvent.UNFLANKED,
        (0, 1, 0, 0, 1, 0, 1): RawEvent.HEMI_FLANKED_LEFT,  # single variant
        (1, 0, 1, 0, 0, 1, 0): RawEvent.HEMI_FLANKED_RIGHT,  # single variant
        (1, 0, 0, 0, 0, 0, 1): RawEvent.REFERENCE,
        (0, 1, 0, 0, 0, 1, 0): RawEvent.TRANSPOSITION,
    }
    
    pattern_tuple = tuple(pattern)
    return patterns.get(pattern_tuple, RawEvent.UNRESOLVED)


def classify_event(
    iso_arch: RawEvent,
    anc_arch: Optional[RawEvent],
    min_jct_cov: int = 5,
) -> tuple[str, List[EventModifier]]:
    """Classify final event based on isolate and ancestor architectures.
    
    Args:
        iso_arch: Isolate architecture
        anc_arch: Ancestor architecture (None if no ancestor)
        min_jct_cov: Minimum coverage threshold
    
    Returns:
        (event_description, list_of_modifiers)
    """
    modifiers = []
    event = iso_arch.value
    
    if anc_arch is None:
        # No ancestor comparison
        return event, modifiers
    
    if iso_arch == anc_arch:
        # Same pattern - ancestral
        modifiers.append(EventModifier.ANCESTRAL)
        event = f"{iso_arch.value} (ancestral)"
    else:
        # Different patterns - check for de novo
        # Simplified logic - full implementation would check specific transitions
        modifiers.append(EventModifier.DE_NOVO)
        event = f"{iso_arch.value} (de novo)"
    
    return event, modifiers

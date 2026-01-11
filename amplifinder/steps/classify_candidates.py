"""Step 13: Final classification of candidates based on junction analysis."""

import numpy as np

from pathlib import Path
from typing import Optional, List, Dict

from amplifinder.data_types import JunctionReadCounts, RecordTypedDf, AnalyzedTnJc2, ClassifiedTnJc2, RawEvent, \
    EventModifier, JunctionType
from amplifinder.steps.base import RecordTypedDfStep


def classify_iso_vs_anc(
    iso_arch: RawEvent,
    anc_arch: RawEvent,
) -> List[EventModifier]:
    """Classify event based on isolate vs ancestor architecture comparison.

    Based on MATLAB classify_candidates.m

    Args:
        iso_arch: Isolate architecture from junction analysis
        anc_arch: Ancestor architecture from junction analysis

    Returns:
        list_of_modifiers
    """
    from_iso_anc_events_to_denovo_left_right = {
        # flanked vs single patterns
        (RawEvent.FLANKED, RawEvent.HEMI_FLANKED_LEFT): (False, True),   # de novo right
        (RawEvent.FLANKED, RawEvent.HEMI_FLANKED_RIGHT): (True, False),  # de novo left
        (RawEvent.FLANKED, RawEvent.REFERENCE): (True, True),            # de novo both
        # hemi-flanked transitions
        (RawEvent.HEMI_FLANKED_LEFT, RawEvent.REFERENCE): (True, False),   # de novo left
        (RawEvent.HEMI_FLANKED_RIGHT, RawEvent.REFERENCE): (False, True),  # de novo right
        # unflanked
        (RawEvent.UNFLANKED, RawEvent.REFERENCE): (False, False),
    }

    modifiers = []

    if iso_arch == anc_arch:
        # Same architecture - ancestral
        modifiers.append(EventModifier.ANCESTRAL)
        return modifiers

    # Look up transition
    denovo_left_right = from_iso_anc_events_to_denovo_left_right.get((iso_arch, anc_arch))

    if denovo_left_right is None:
        # Unrecognized transition
        return modifiers

    denovo_left, denovo_right = denovo_left_right

    if denovo_left:
        modifiers.append(EventModifier.DENOVO_LEFT)
    if denovo_right:
        modifiers.append(EventModifier.DENOVO_RIGHT)

    return modifiers


def classify_architecture(jc_calls: Dict[JunctionType, Optional[bool]]) -> RawEvent:
    """Classify junction architecture from read coverage patterns.

    Args:
        jc_calls: Dict mapping JunctionType to True (covered), False (not covered), or None (undetermined)

    Returns:
        RawEvent classification

    Based on MATLAB classify_candidates.m
    """

    # Amplicon structure:
    # 
    # ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~
    # 
    # ~~-==  ~~->>  ==->>  ==-==  >>-==  >>-~~  ==-~~
    #   1      2      3      4      5      6      7
    # 
    # Junction patterns (1=covered, 0=not):

    architecture_patterns: Dict[tuple[int, int, int, int, int, int, int], RawEvent] = {
    #    -- Junction Type --
    #    1  2  3  4  5  6  7   event
        (0, 1, 1, 0, 1, 1, 0): RawEvent.FLANKED,
        (0, 1, 1, 0, 1, 0, 1): RawEvent.HEMI_FLANKED_LEFT,
        (1, 0, 1, 0, 1, 1, 0): RawEvent.HEMI_FLANKED_RIGHT,
        (1, 0, 1, 0, 1, 0, 1): RawEvent.UNFLANKED,
        (0, 1, 0, 0, 1, 0, 1): RawEvent.HEMI_FLANKED_LEFT,   # single variant
        (1, 0, 1, 0, 0, 1, 0): RawEvent.HEMI_FLANKED_RIGHT,  # single variant
        (1, 0, 0, 0, 0, 0, 1): RawEvent.REFERENCE,
        (0, 1, 0, 0, 0, 1, 0): RawEvent.TRANSPOSITION,
    }

    pattern = tuple(jc_calls[jt] for jt in JunctionType.sorted())

    if None in pattern:
        return RawEvent.UNRESOLVED
    return architecture_patterns.get(pattern, RawEvent.UNRESOLVED)


class ClassifyTnJc2CandidatesStep(RecordTypedDfStep[ClassifiedTnJc2]):
    """Final classification of candidates based on iso/anc comparison.

    This step refines the event classification from AnalyzeAlignmentsStep
    by performing detailed iso vs anc pattern comparison.

    Classification depends on run type:
    - has_ancestor=False: limited classification (no de novo detection)
    - has_ancestor=True: full iso vs anc pattern comparison
    """

    def __init__(
        self,
        analyzed_tnjc2s: RecordTypedDf[AnalyzedTnJc2],
        output_dir: Path,
        force: Optional[bool] = None,
    ):
        self.analyzed_tnjc2s = analyzed_tnjc2s
        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[ClassifiedTnJc2]:
        """Reclassify events based on iso/anc comparison."""
        classified_records = []

        for analyzed_tnjc2 in self.analyzed_tnjc2s:
            # Compute architectures from junction coverage
            iso_arch = classify_architecture(analyzed_tnjc2.jc_calls)

            anc_arch, modifiers = None, []
            if analyzed_tnjc2.jc_calls_anc:
                anc_arch = classify_architecture(analyzed_tnjc2.jc_calls_anc)
                modifiers = classify_iso_vs_anc(iso_arch, anc_arch)

            # Create ClassifiedTnJc2 with architecture and classification
            classified = ClassifiedTnJc2.from_other(
                analyzed_tnjc2,
                isolate_architecture=iso_arch,
                ancestor_architecture=anc_arch,
                event_modifiers=modifiers,
            )

            classified_records.append(classified)

        return RecordTypedDf.from_records(classified_records, ClassifiedTnJc2)

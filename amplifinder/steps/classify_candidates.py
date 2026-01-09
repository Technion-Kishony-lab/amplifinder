"""Step 13: Final classification of candidates based on junction analysis."""

from pathlib import Path
from typing import Optional, List

from amplifinder.data_types import JunctionReadCounts, RecordTypedDf, AnalyzedTnJc2, ClassifiedTnJc2, RawEvent, \
    EventModifier, JunctionType
from amplifinder.steps.base import RecordTypedDfStep

# Iso/Anc pattern transition rules (from MATLAB classify_candidates.m)
# Format: (iso_pattern, anc_pattern) -> (de_novo_left, de_novo_right)
ISO_ANC_TRANSITIONS = {
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

# Junction coverage patterns -> RawEvent classification
# Pattern: (jct1, jct2, ..., jct7) where 1=covered, 0=not
ARCHITECTURE_PATTERNS = {
    (0, 1, 1, 0, 1, 1, 0): RawEvent.FLANKED,
    (0, 1, 1, 0, 1, 0, 1): RawEvent.HEMI_FLANKED_LEFT,
    (1, 0, 1, 0, 1, 1, 0): RawEvent.HEMI_FLANKED_RIGHT,
    (1, 0, 1, 0, 1, 0, 1): RawEvent.UNFLANKED,
    (0, 1, 0, 0, 1, 0, 1): RawEvent.HEMI_FLANKED_LEFT,   # single variant
    (1, 0, 1, 0, 0, 1, 0): RawEvent.HEMI_FLANKED_RIGHT,  # single variant
    (1, 0, 0, 0, 0, 0, 1): RawEvent.REFERENCE,
    (0, 1, 0, 0, 0, 1, 0): RawEvent.TRANSPOSITION,
}


def classify_iso_vs_anc(
    iso_arch: RawEvent,
    anc_arch: Optional[RawEvent],
    min_jct_cov: int = 5,
) -> List[EventModifier]:
    """Classify event based on isolate vs ancestor architecture comparison.

    Based on MATLAB classify_candidates.m

    Args:
        iso_arch: Isolate architecture from junction analysis
        anc_arch: Ancestor architecture (None if no ancestor run)
        min_jct_cov: Minimum coverage threshold

    Returns:
        list_of_modifiers
    """
    modifiers = []

    if anc_arch is None:
        # No ancestor - can't determine de novo status
        return modifiers

    if iso_arch == anc_arch:
        # Same architecture - ancestral
        modifiers.append(EventModifier.ANCESTRAL)
        return modifiers

    # Look up transition
    transition = ISO_ANC_TRANSITIONS.get((iso_arch, anc_arch))

    if transition is None:
        # Unrecognized transition
        return modifiers

    de_novo_left, de_novo_right = transition

    if de_novo_left:
        modifiers.append(EventModifier.DENOVO_LEFT)
    if de_novo_right:
        modifiers.append(EventModifier.DENOVO_RIGHT)

    return modifiers


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
    pattern = tuple(1 if jc.spanning >= min_jct_cov else 0 for jc in jc_cov)
    return ARCHITECTURE_PATTERNS.get(pattern, RawEvent.UNRESOLVED)


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
        has_ancestor: bool = False,
        min_jct_cov: int = 5,
        force: Optional[bool] = None,
    ):
        self.analyzed_tnjc2s = analyzed_tnjc2s
        self.has_ancestor = has_ancestor
        self.min_jct_cov = min_jct_cov

        super().__init__(
            output_dir=output_dir,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[ClassifiedTnJc2]:
        """Reclassify events based on iso/anc comparison."""
        classified_records = []

        for analyzed_tnjc2 in self.analyzed_tnjc2s:
            # Compute architectures from junction coverage
            iso_cov_list = [analyzed_tnjc2.jc_cov[jt] for jt in JunctionType.sorted()]
            iso_arch = classify_architecture(iso_cov_list, self.min_jct_cov)

            anc_arch = None
            if self.has_ancestor and analyzed_tnjc2.jc_cov_anc:
                anc_cov_list = [analyzed_tnjc2.jc_cov_anc[jt] for jt in JunctionType.sorted()]
                anc_arch = classify_architecture(anc_cov_list, self.min_jct_cov)

            # Reclassify
            modifiers = classify_iso_vs_anc(iso_arch, anc_arch, self.min_jct_cov)

            # Create ClassifiedTnJc2 with architecture and classification
            classified = ClassifiedTnJc2.from_other(
                analyzed_tnjc2,
                isolate_architecture=iso_arch,
                ancestor_architecture=anc_arch,
                event_modifiers=modifiers,
            )

            classified_records.append(classified)

        return RecordTypedDf.from_records(classified_records, ClassifiedTnJc2)

    def report_output_message(self, output: RecordTypedDf[ClassifiedTnJc2], *, from_cache: bool) -> Optional[str]:
        if self.has_ancestor:
            ancestral = sum(1 for r in output if EventModifier.ANCESTRAL in r.event_modifiers)
            de_novo = sum(1 for r in output if EventModifier.DENOVO_LEFT in r.event_modifiers
                          or EventModifier.DENOVO_RIGHT in r.event_modifiers)
            return f"Classification: {ancestral} ancestral, {de_novo} de novo, {len(output)} total"
        return f"Classification: {len(output)} candidates (no ancestor comparison)"

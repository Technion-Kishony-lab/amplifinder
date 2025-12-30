"""Step 13: Final classification of candidates based on junction analysis."""

from pathlib import Path
from typing import Optional, List, Tuple

from amplifinder.data_types import (
    RecordTypedDf, AnalyzedTnJc2, RawEvent, EventModifier,
)
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


def classify_iso_vs_anc(
    iso_arch: RawEvent,
    anc_arch: Optional[RawEvent],
    min_jct_cov: int = 5,
) -> Tuple[str, List[EventModifier]]:
    """Classify event based on isolate vs ancestor architecture comparison.

    Based on MATLAB classify_candidates.m

    Args:
        iso_arch: Isolate architecture from junction analysis
        anc_arch: Ancestor architecture (None if no ancestor run)
        min_jct_cov: Minimum coverage threshold

    Returns:
        (event_name, list_of_modifiers)
    """
    modifiers = []

    if anc_arch is None:
        # No ancestor - can't determine de novo status
        return iso_arch.value, modifiers

    if iso_arch == anc_arch:
        # Same architecture - ancestral
        modifiers.append(EventModifier.ANCESTRAL)
        return f"{iso_arch.value} (ancestral)", modifiers

    # Look up transition
    transition = ISO_ANC_TRANSITIONS.get((iso_arch, anc_arch))

    if transition is None:
        # Unrecognized transition
        return f"{iso_arch.value} (unresolved iso-anc pair)", modifiers

    de_novo_left, de_novo_right = transition

    event_parts = [iso_arch.value]
    if de_novo_left:
        event_parts.append("de novo left")
        modifiers.append(EventModifier.DE_NOVO)
    if de_novo_right:
        event_parts.append("de novo right")
        if EventModifier.DE_NOVO not in modifiers:
            modifiers.append(EventModifier.DE_NOVO)

    return " ".join(event_parts), modifiers


class ClassifyTnJc2CandidatesStep(RecordTypedDfStep[AnalyzedTnJc2]):
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

    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Reclassify events based on iso/anc comparison."""
        classified_records = []

        for analyzed_tnjc2 in self.analyzed_tnjc2s:
            iso_arch = analyzed_tnjc2.isolate_architecture
            anc_arch = analyzed_tnjc2.ancestor_architecture if self.has_ancestor else None

            # Reclassify
            event, modifiers = classify_iso_vs_anc(iso_arch, anc_arch, self.min_jct_cov)

            # Update record with new classification
            classified = analyzed_tnjc2.model_copy(update={
                "event": event,
                "event_modifiers": modifiers,
            })

            classified_records.append(classified)

        result = RecordTypedDf.from_records(classified_records, AnalyzedTnJc2)

        return result

    def report_output_message(self, output: RecordTypedDf[AnalyzedTnJc2], *, from_cache: bool) -> Optional[str]:
        if self.has_ancestor:
            ancestral = sum(1 for r in output if EventModifier.ANCESTRAL in r.event_modifiers)
            de_novo = sum(1 for r in output if EventModifier.DE_NOVO in r.event_modifiers)
            return f"Classification: {ancestral} ancestral, {de_novo} de novo, {len(output)} total"
        return f"Classification: {len(output)} candidates (no ancestor comparison)"

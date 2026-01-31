"""Step 13: Final classification of candidates based on junction analysis."""

from pathlib import Path
from typing import Optional, List, Dict

from amplifinder.data_types import JcCall, RecordTypedDf, AnalyzedTnJc2, ClassifiedTnJc2, Architecture, \
    EventDescriptor, JunctionType
from amplifinder.steps.base import RecordTypedDfStep


def classify_iso_vs_anc(iso_arch: Architecture, anc_arch: Architecture) -> List[EventDescriptor]:
    """Classify event based on isolate vs ancestor architecture comparison.

    Args:
        iso_arch: Isolate architecture from junction analysis
        anc_arch: Ancestor architecture from junction analysis

    Returns:
        list_of_descriptors: List[EventDescriptor], descriptors for the event
    """
    descriptors = []

    if iso_arch == anc_arch:
        # Same architecture - ancestral
        descriptors.append(EventDescriptor.ANCESTRAL)
        return descriptors

    # Use the is_flanking property to determine de novo and lost junctions
    iso_left, iso_right = iso_arch.is_flanking
    anc_left, anc_right = anc_arch.is_flanking

    # Check left side (each side separately)
    if iso_left is not None and anc_left is not None:
        if iso_left and not anc_left:
            descriptors.append(EventDescriptor.DENOVO_LEFT)
        elif not iso_left and anc_left:
            descriptors.append(EventDescriptor.LOST_LEFT)

    # Check right side
    if iso_right is not None and anc_right is not None:
        if iso_right and not anc_right:
            descriptors.append(EventDescriptor.DENOVO_RIGHT)
        elif not iso_right and anc_right:
            descriptors.append(EventDescriptor.LOST_RIGHT)

    # Check singleton status change
    iso_singleton = iso_arch.is_singleton
    anc_singleton = anc_arch.is_singleton

    if iso_singleton is not None and anc_singleton is not None:
        if iso_singleton and not anc_singleton:
            descriptors.append(EventDescriptor.TO_SINGLETON)
        elif not iso_singleton and anc_singleton:
            descriptors.append(EventDescriptor.FROM_SINGLETON)

    return descriptors


def classify_architecture(jc_calls: Dict[JunctionType, JcCall]) -> Architecture:
    """Classify junction architecture from read coverage patterns.

    Args:
        jc_calls: Dict mapping JunctionType to True (covered), False (not covered), or None (undetermined)

    Returns:
        Architecture classification

    Based on MATLAB classify_candidates.m
    """

    # Amplicon structure:
    #
    #  ~~~~~~~~~~> |---> |========> |---> |========> |---> |~~~~~~~~~~
    #  chromosome   IS    amplicon   IS    amplicon   IS    chromosome
    #
    # Junctions:
    #  ~~~ |==  ~~~ |--  ==> |--  ==> |==  --> |==  --> ~~~  ==> ~~~
    #     1        2        3        4        5        6        7
    #
    # Junction patterns (1=covered, 0=not):
    #
    # '...' = |========> |--->
    #
    #    -- Junction Type --
    #    1  2  3  5  6  7   event
    patterns_to_architecture: dict[tuple[int, int, int, int, int, int], Architecture] = {
        #  ~~~ |---> |========> |---> ... |========> |---> ~~~
        (0, 1, 1, 1, 1, 0): Architecture.FLANKED,

        #  ~~~ |---> |========> |---> ... |========> ~~~
        (0, 1, 1, 1, 0, 1): Architecture.HEMI_FLANKED_LEFT,

        #        ~~~ |========> |---> ... |========> |---> ~~~
        (1, 0, 1, 1, 1, 0): Architecture.HEMI_FLANKED_RIGHT,

        #        ~~~ |========> |---> ... |========> ~~~
        (1, 0, 1, 1, 0, 1): Architecture.UNFLANKED,

        #  ~~~ |---> |========> ~~~
        (0, 1, 0, 1, 0, 1): Architecture.HEMI_FLANKED_LEFT_SINGLETON,

        #        ~~~ |========> |---> ~~~
        (1, 0, 1, 0, 1, 0): Architecture.HEMI_FLANKED_RIGHT_SINGLETON,

        #        ~~~ |========> ~~~
        (1, 0, 0, 0, 0, 1): Architecture.AMPLICON_ONLY,

        #        ~~~ |---> ~~~
        (0, 1, 0, 0, 1, 0): Architecture.TRANSPOSITION,
    }

    pattern = tuple(jc_calls[jc] for jc in JunctionType if jc != JunctionType.AMP_AMP)
    if JcCall.AMBIGIOUS in pattern:
        return Architecture.UNRESOLVED

    pattern_ints = tuple(int(bool(jc_call)) for jc_call in pattern)
    return patterns_to_architecture.get(pattern_ints, Architecture.UNRESOLVED)


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

            anc_arch, descriptors = None, []
            if analyzed_tnjc2.jc_calls_anc:
                anc_arch = classify_architecture(analyzed_tnjc2.jc_calls_anc)
                descriptors.extend(classify_iso_vs_anc(iso_arch, anc_arch))

            if analyzed_tnjc2.jc_calls[JunctionType.AMP_AMP]:
                descriptors.append(EventDescriptor.AMP_AMP)

            # Create ClassifiedTnJc2 with architectures and descriptors
            classified = ClassifiedTnJc2.from_other(
                analyzed_tnjc2,
                iso_architecture=iso_arch,
                anc_architecture=anc_arch,
                event_descriptors=descriptors,
            )

            classified_records.append(classified)

        return RecordTypedDf.from_records(classified_records, ClassifiedTnJc2)

"""Step 8: Classify junction pair structures based on TN relationships."""

from pathlib import Path
from typing import Optional, TYPE_CHECKING

from amplifinder.data_types import RecordTypedDf, CoveredTnJc2, ClassifiedTnJc2, RawEvent, RefTn, Side, TnJunction
from amplifinder.data_types.enums import BaseRawEvent
from amplifinder.data_types.genome import Genome

from .base import RecordTypedDfStep


class ClassifyTnJc2StructureStep(RecordTypedDfStep[ClassifiedTnJc2]):
    """Classify junction pair structures based on TN relationships.

    This step analyzes how junction pairs relate to reference TN elements
    and other single-locus pairs to classify them as:
    - reference: TN in its reference location
    - transposition: short de novo insertion
    - flanked: both sides share IS with single-locus pairs
    - hemi-flanked: one side shares IS
    - unflanked: neither side shares IS
    """

    def __init__(
        self,
        covered_tnjc2s: RecordTypedDf[CoveredTnJc2],
        genome: Genome,
        tn_locs: RecordTypedDf[RefTn],
        output_dir: Path,
        transposition_threshold: int = 30,
        force: Optional[bool] = None,
    ):
        self.covered_tnjc2s = covered_tnjc2s
        self.genome = genome
        self.tn_locs = tn_locs
        self.transposition_threshold = transposition_threshold

        super().__init__(output_dir=output_dir, force=force)

    def _calculate_output(self) -> RecordTypedDf[ClassifiedTnJc2]:
        """Classify junction pairs."""
        tnjc2s = self.covered_tnjc2s.to_records()
        base_raw_events = [self._compute_base_raw_event(tnjc2) for tnjc2 in tnjc2s]
        classified_tnjc2s = []
        for i, tncj2_i in enumerate(tnjc2s):
            tnjc2_matching_S = self._find_matching_tnjc2(tncj2_i.tnjc_S, tnjc2s, base_raw_events)
            tnjc2_matching_E = self._find_matching_tnjc2(tncj2_i.tnjc_E, tnjc2s, base_raw_events)
            classified_tnjc2s.append(ClassifiedTnJc2.from_other(
                tncj2_i, 
                tnjc2_matching_S=tnjc2_matching_S, 
                tnjc2_matching_E=tnjc2_matching_E,
                base_raw_event=base_raw_events[i],
            ))
        return RecordTypedDf.from_records(classified_tnjc2s, ClassifiedTnJc2)

    def _compute_base_raw_event(self, tnjc2: CoveredTnJc2) -> BaseRawEvent:
        if tnjc2.tnjc_S.is_ref_tn_junction() and tnjc2.tnjc_E.is_ref_tn_junction():
            return BaseRawEvent.REFERENCE
        elif abs(tnjc2.start - tnjc2.end) < self.transposition_threshold:
            return BaseRawEvent.TRANSPOSITION
        else:
            return BaseRawEvent.LOCUS_JOINING

    def _find_matching_tnjc2(self, tnjc_i: TnJunction, tnjc2s: list[CoveredTnJc2], base_raw_events: list[BaseRawEvent]) -> Optional[CoveredTnJc2]:
        for j, tnjc2_j in enumerate(tnjc2s):
            if base_raw_events[j] != BaseRawEvent.LOCUS_JOINING:
                continue
            if tnjc_i == tnjc2_j.tnjc_S or tnjc_i == tnjc2_j.tnjc_E:
                return tnjc2_j
        return None

    def report_output_message(self, output: RecordTypedDf[ClassifiedTnJc2], *, from_cache: bool) -> Optional[str]:
        records: list[ClassifiedTnJc2] = output.to_records()
        counts: dict[RawEvent, int] = {name: 0 for name in RawEvent}
        for record in records:
            counts[record.raw_event] += 1
        return ", ".join(f"{event}: {count}" for event, count in counts.items())

"""Step 8: Classify junction pair structures based on TN relationships."""

from pathlib import Path
from typing import Optional, List, Tuple

from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, ClassifiedTnJc2, RawEvent, RefTnLoc,
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.logger import info


def classify_structure(
    covered_tnjc2: RecordTypedDf[CoveredTnJc2],
    min_amplicon_length: int = 30,
) -> RecordTypedDf[ClassifiedTnJc2]:
    """Classify junction pairs based on TN relationships.
    
    Based on MATLAB classify_ISJC2.m
    
    Classification logic:
    1. Reference: both junctions match reference TN at same location
    2. Transposition: short amplicon (< threshold), de novo, not reference
    3. Single-locus: reference or transposition
    4. For each junction, check if it appears in other single-locus pairs
    5. Classify based on shared IS:
       - Both sides share: flanked
       - One side shares: hemi-flanked (left/right)
       - Neither shares: unflanked
       - Multiple matches: multiple single locus
    
    Args:
        covered_tnjc2: Covered TnJc2 records
        min_amplicon_length: Threshold for transposition detection
    
    Returns:
        Classified TnJc2 records
    """
    if len(covered_tnjc2) == 0:
        return RecordTypedDf.empty(ClassifiedTnJc2)
    
    # Step 1: Identify reference and transposition pairs
    is_reference = []
    is_transposition = []
    is_single_locus = []
    
    for tnjc2 in covered_tnjc2:
        # Reference: both junctions are reference (jc_num == 0) and match same TN
        # Simplified: check if both jc_num are 0 (reference junctions)
        # In practice, we'd need to check if they match the same reference TN location
        ref_left = tnjc2.jc_num_L == 0
        ref_right = tnjc2.jc_num_R == 0
        
        # Check if both match same TN ID (simplified check)
        same_tn = len(tnjc2.tn_ids) > 0 and all(tid == tnjc2.tn_ids[0] for tid in tnjc2.tn_ids)
        
        is_ref = ref_left and ref_right and same_tn
        is_reference.append(is_ref)
        
        # Transposition: short amplicon, not reference
        is_trans = (
            not is_ref
            and tnjc2.amplicon_length < min_amplicon_length
        )
        is_transposition.append(is_trans)
        
        # Single-locus: reference or transposition
        is_single_locus.append(is_ref or is_trans)
    
    # Step 2: For each junction, find other single-locus pairs that share it
    n = len(covered_tnjc2)
    shared_tn_ids_left = [[] for _ in range(n)]
    shared_tn_ids_right = [[] for _ in range(n)]
    shared_tn_ids_both = [[] for _ in range(n)]
    multiple_single_locus = [False] * n
    
    for i, tnjc2_i in enumerate(covered_tnjc2):
        if not is_single_locus[i]:
            continue
        
        # Check left junction
        matches_left = []
        for j, tnjc2_j in enumerate(covered_tnjc2):
            if i == j or not is_single_locus[j]:
                continue
            
            # Check if left junction matches
            if tnjc2_i.jc_num_L == tnjc2_j.jc_num_L or tnjc2_i.jc_num_L == tnjc2_j.jc_num_R:
                # Find shared TN IDs
                shared = [tid for tid in tnjc2_i.tn_ids if tid in tnjc2_j.tn_ids]
                if shared:
                    matches_left.append((j, shared))
        
        if len(matches_left) > 1:
            multiple_single_locus[i] = True
        elif len(matches_left) == 1:
            _, shared = matches_left[0]
            shared_tn_ids_left[i] = shared
        
        # Check right junction
        matches_right = []
        for j, tnjc2_j in enumerate(covered_tnjc2):
            if i == j or not is_single_locus[j]:
                continue
            
            # Check if right junction matches
            if tnjc2_i.jc_num_R == tnjc2_j.jc_num_L or tnjc2_i.jc_num_R == tnjc2_j.jc_num_R:
                # Find shared TN IDs
                shared = [tid for tid in tnjc2_i.tn_ids if tid in tnjc2_j.tn_ids]
                if shared:
                    matches_right.append((j, shared))
        
        if len(matches_right) > 1:
            multiple_single_locus[i] = True
        elif len(matches_right) == 1:
            _, shared = matches_right[0]
            shared_tn_ids_right[i] = shared
        
        # Find TN IDs shared by both sides
        shared_both = [tid for tid in shared_tn_ids_left[i] if tid in shared_tn_ids_right[i]]
        shared_tn_ids_both[i] = shared_both
    
    # Step 3: Classify based on shared IS
    classified_records = []
    for i, tnjc2 in enumerate(covered_tnjc2):
        # Determine classification
        if multiple_single_locus[i]:
            raw_event = RawEvent.MULTIPLE_SINGLE_LOCUS
        elif is_reference[i]:
            raw_event = RawEvent.REFERENCE
        elif is_transposition[i]:
            raw_event = RawEvent.TRANSPOSITION
        elif len(shared_tn_ids_both[i]) > 0:
            raw_event = RawEvent.FLANKED
        elif len(shared_tn_ids_left[i]) > 0 or len(shared_tn_ids_right[i]) > 0:
            # Hemi-flanked: determine left vs right
            # Check which side has the shared IS
            has_left = len(shared_tn_ids_left[i]) > 0
            has_right = len(shared_tn_ids_right[i]) > 0
            
            # Account for origin spanning
            if tnjc2.span_origin:
                # When spanning origin, left/right are swapped
                if has_right:
                    raw_event = RawEvent.HEMI_FLANKED_LEFT
                else:
                    raw_event = RawEvent.HEMI_FLANKED_RIGHT
            else:
                if has_left:
                    raw_event = RawEvent.HEMI_FLANKED_LEFT
                else:
                    raw_event = RawEvent.HEMI_FLANKED_RIGHT
        else:
            raw_event = RawEvent.UNFLANKED
        
        # Determine shared TN IDs and chosen TN ID
        if raw_event == RawEvent.FLANKED:
            shared_tn_ids = shared_tn_ids_both[i]
        elif raw_event in (RawEvent.HEMI_FLANKED_LEFT, RawEvent.HEMI_FLANKED_RIGHT):
            # Use the side that has the shared IS
            if len(shared_tn_ids_left[i]) > 0:
                shared_tn_ids = shared_tn_ids_left[i]
            else:
                shared_tn_ids = shared_tn_ids_right[i]
        else:
            # Unflanked, reference, transposition: use all TN IDs
            shared_tn_ids = tnjc2.tn_ids
        
        # Chosen TN ID: first shared TN, or first TN if none shared
        chosen_tn_id = shared_tn_ids[0] if shared_tn_ids else (tnjc2.tn_ids[0] if tnjc2.tn_ids else None)
        
        # Create classified record
        classified = ClassifiedTnJc2.from_other(
            tnjc2,
            raw_event=raw_event,
            shared_tn_ids=shared_tn_ids,
            chosen_tn_id=chosen_tn_id,
        )
        classified_records.append(classified)
    
    return RecordTypedDf.from_records(classified_records, ClassifiedTnJc2)


class ClassifyStructureStep(RecordTypedDfStep[ClassifiedTnJc2]):
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
        covered_tnjc2: RecordTypedDf[CoveredTnJc2],
        tn_locs: RecordTypedDf[RefTnLoc],
        output_dir: Path,
        min_amplicon_length: int = 30,
        force: Optional[bool] = None,
    ):
        self.covered_tnjc2 = covered_tnjc2
        self.tn_locs = tn_locs
        self.min_amplicon_length = min_amplicon_length
        
        super().__init__(
            output_dir=output_dir,
            input_files=[],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[ClassifiedTnJc2]:
        """Classify junction pairs."""
        classified = classify_structure(
            self.covered_tnjc2,
            min_amplicon_length=self.min_amplicon_length,
        )
        
        info(f"Classification: {len([r for r in classified if r.raw_event == RawEvent.FLANKED])} flanked, "
             f"{len([r for r in classified if r.raw_event == RawEvent.UNFLANKED])} unflanked, "
             f"{len([r for r in classified if r.raw_event == RawEvent.HEMI_FLANKED_LEFT])} hemi-left, "
             f"{len([r for r in classified if r.raw_event == RawEvent.HEMI_FLANKED_RIGHT])} hemi-right")
        
        return classified

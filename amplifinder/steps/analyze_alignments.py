"""Step 12: Analyze read alignments to synthetic junctions."""

from pathlib import Path
from typing import Optional, List

from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, AnalyzedTnJc2, RawEvent, EventModifier,
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.utils.bam import get_junction_coverage, JunctionReadCounts
from amplifinder.logger import info


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


class AnalyzeAlignmentsStep(RecordTypedDfStep[AnalyzedTnJc2]):
    """Analyze read alignments to classify junction architectures.
    
    Analysis depends on run type:
    - has_ancestor=False: analyze isolate BAM only
    - has_ancestor=True: analyze both isolate and ancestor BAMs, compare patterns
    """

    def __init__(
        self,
        candidates: RecordTypedDf[FilteredTnJc2],
        output_dir: Path,
        anc_output_dir: Optional[Path] = None,
        read_length: int = 150,
        req_overlap: int = 12,
        min_jct_cov: int = 5,
        has_ancestor: bool = False,
        force: Optional[bool] = None,
    ):
        self.candidates = candidates
        self.output_dir = Path(output_dir)
        self.anc_output_dir = Path(anc_output_dir) if anc_output_dir else None
        self.read_length = read_length
        self.req_overlap = req_overlap
        self.min_jct_cov = min_jct_cov
        self.has_ancestor = has_ancestor
        
        # Input files are the BAM files from alignment step
        input_files = []
        for cand in candidates:
            analysis_dir = output_dir / cand.analysis_dir
            input_files.append(analysis_dir / "iso.sorted.bam")
            if has_ancestor:
                # Ancestor BAM is in ancestor folder (as iso.sorted.bam)
                if not self.anc_output_dir:
                    raise ValueError("anc_output_dir must be provided when has_ancestor=True")
                anc_analysis_dir = self.anc_output_dir / cand.analysis_dir
                input_files.append(anc_analysis_dir / "iso.sorted.bam")
        
        super().__init__(
            output_dir=output_dir,
            input_files=input_files,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Analyze alignments for each candidate."""
        junction_length = self.read_length * 2
        
        analyzed_records = []
        for cand in self.candidates:
            analysis_dir = self.output_dir / cand.analysis_dir
            
            # Get isolate junction coverage
            iso_bam = analysis_dir / "iso.sorted.bam"
            if not iso_bam.exists():
                info(f"Skipping {cand.analysis_dir}: no iso.sorted.bam")
                continue
            
            iso_jc_cov = get_junction_coverage(
                iso_bam, junction_length, self.read_length, self.req_overlap
            )
            
            # Classify isolate architecture
            iso_arch = classify_architecture(iso_jc_cov, self.min_jct_cov)
            
            # Get ancestor junction coverage if available
            # Ancestor BAM is stored in ancestor folder as iso.sorted.bam
            anc_jc_cov = None
            anc_arch = None
            if self.has_ancestor:
                if not self.anc_output_dir:
                    raise ValueError("anc_output_dir must be provided when has_ancestor=True")
                anc_analysis_dir = self.anc_output_dir / cand.analysis_dir
                anc_bam = anc_analysis_dir / "iso.sorted.bam"
                
                if anc_bam.exists():
                    anc_jc_cov = get_junction_coverage(
                        anc_bam, junction_length, self.read_length, self.req_overlap
                    )
                    anc_arch = classify_architecture(anc_jc_cov, self.min_jct_cov)
            
            # Classify final event
            event, modifiers = classify_event(iso_arch, anc_arch, self.min_jct_cov)
            
            # Build AnalyzedTnJc2 record from candidate + new fields
            analyzed = AnalyzedTnJc2.from_other(
                cand,
                jc_cov_left=[jc.left for jc in iso_jc_cov],
                jc_cov_right=[jc.right for jc in iso_jc_cov],
                jc_cov_spanning=[jc.spanning for jc in iso_jc_cov],
                anc_jc_cov_left=[jc.left for jc in anc_jc_cov] if anc_jc_cov else None,
                anc_jc_cov_right=[jc.right for jc in anc_jc_cov] if anc_jc_cov else None,
                anc_jc_cov_spanning=[jc.spanning for jc in anc_jc_cov] if anc_jc_cov else None,
                isolate_architecture=iso_arch,
                ancestor_architecture=anc_arch,
                event=event,
                event_modifiers=modifiers,
            )
            analyzed_records.append(analyzed)
        
        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)

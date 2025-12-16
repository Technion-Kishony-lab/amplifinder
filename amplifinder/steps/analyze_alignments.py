"""Step 12: Analyze read alignments to synthetic junctions."""

from pathlib import Path
from typing import Optional, List

from amplifinder.data_types import (
    RecordTypedDF, CandidateTnJc2, AnalyzedTnJc2, RawEvent, EventModifier,
)
from amplifinder.steps.base import Step
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


class AnalyzeAlignmentsStep(Step[RecordTypedDF[AnalyzedTnJc2]]):
    """Analyze read alignments to classify junction architectures.
    
    Analysis depends on run type:
    - has_ancestor=False: analyze isolate BAM only
    - has_ancestor=True: analyze both isolate and ancestor BAMs, compare patterns
    """

    def __init__(
        self,
        candidates: RecordTypedDF[CandidateTnJc2],
        output_dir: Path,
        read_length: int = 150,
        req_overlap: int = 12,
        min_jct_cov: int = 5,
        has_ancestor: bool = False,
        force: Optional[bool] = None,
    ):
        self.candidates = candidates
        self.output_dir = Path(output_dir)
        self.read_length = read_length
        self.req_overlap = req_overlap
        self.min_jct_cov = min_jct_cov
        self.has_ancestor = has_ancestor
        
        self.output_file = output_dir / "tn_jc2_analyzed.csv"
        
        # Input files are the BAM files from alignment step
        input_files = []
        for _, row in candidates.df.iterrows():
            analysis_dir = output_dir / row["analysis_dir"]
            input_files.append(analysis_dir / "iso.sorted.bam")
            if has_ancestor:
                input_files.append(analysis_dir / "anc.sorted.bam")
        
        super().__init__(
            input_files=input_files,
            output_files=[self.output_file],
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDF[AnalyzedTnJc2]:
        """Analyze alignments for each candidate."""
        junction_length = self.read_length * 2
        
        analyzed_records = []
        for _, row in self.candidates.df.iterrows():
            analysis_dir = self.output_dir / row["analysis_dir"]
            
            # Get isolate junction coverage
            iso_bam = analysis_dir / "iso.sorted.bam"
            if not iso_bam.exists():
                info(f"Skipping {row['analysis_dir']}: no iso.sorted.bam")
                continue
            
            iso_jc_cov = get_junction_coverage(
                iso_bam, junction_length, self.read_length, self.req_overlap
            )
            
            # Classify isolate architecture
            iso_arch = classify_architecture(iso_jc_cov, self.min_jct_cov)
            
            # Get ancestor junction coverage if available
            anc_jc_cov = None
            anc_arch = None
            if self.has_ancestor:
                anc_bam = analysis_dir / "anc.sorted.bam"
                if anc_bam.exists():
                    anc_jc_cov = get_junction_coverage(
                        anc_bam, junction_length, self.read_length, self.req_overlap
                    )
                    anc_arch = classify_architecture(anc_jc_cov, self.min_jct_cov)
            
            # Classify final event
            event, modifiers = classify_event(iso_arch, anc_arch, self.min_jct_cov)
            
            # Build AnalyzedTnJc2 record
            analyzed = AnalyzedTnJc2(
                # From CandidateTnJc2
                jc_num_L=row["jc_num_L"],
                jc_num_R=row["jc_num_R"],
                scaf_chr=row["scaf_chr"],
                pos_chr_L=row["pos_chr_L"],
                pos_chr_R=row["pos_chr_R"],
                pos_tn_L=row["pos_tn_L"],
                pos_tn_R=row["pos_tn_R"],
                dir_chr_L=row["dir_chr_L"],
                dir_chr_R=row["dir_chr_R"],
                dir_tn_L=row["dir_tn_L"],
                dir_tn_R=row["dir_tn_R"],
                tn_ids=row["tn_ids"],
                tn_orientations=row["tn_orientations"],
                span_origin=row["span_origin"],
                amplicon_length=row["amplicon_length"],
                complementary_length=row["complementary_length"],
                # From CoveredTnJc2
                ref_name=row["ref_name"],
                iso_name=row["iso_name"],
                anc_name=row.get("anc_name"),
                amplicon_coverage=row["amplicon_coverage"],
                genome_coverage=row["genome_coverage"],
                copy_number=row["copy_number"],
                amplicon_coverage_mode=row["amplicon_coverage_mode"],
                copy_number_ratio=row.get("copy_number_ratio"),
                # From ClassifiedTnJc2
                raw_event=row["raw_event"],
                shared_tn_ids=row["shared_tn_ids"],
                chosen_tn_id=row.get("chosen_tn_id"),
                # From CandidateTnJc2
                analysis_dir=row["analysis_dir"],
                # New fields
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
        
        return RecordTypedDF.from_records(analyzed_records, AnalyzedTnJc2)

    def _save_output(self, output: RecordTypedDF[AnalyzedTnJc2]) -> None:
        """Save analyzed results to CSV."""
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[AnalyzedTnJc2]:
        """Load analyzed results from CSV."""
        return RecordTypedDF.from_csv(self.output_file, AnalyzedTnJc2)


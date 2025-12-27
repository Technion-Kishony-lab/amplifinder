"""Step 12: Analyze read alignments to synthetic junctions."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, AnalyzedTnJc2,
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.steps.analyze_alignments.parse_bam import get_junction_coverage
from amplifinder.steps.analyze_alignments.classify import classify_architecture, classify_event
from amplifinder.logger import info


class AnalyzeTnJc2AlignmentsStep(RecordTypedDfStep[AnalyzedTnJc2]):
    """Analyze read alignments to classify junction architectures.

    Analysis depends on run type:
    - has_ancestor=False: analyze isolate BAM only
    - has_ancestor=True: analyze both isolate and ancestor BAMs, compare patterns
    """

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        output_dir: Path,
        anc_output_dir: Optional[Path] = None,
        read_length: int = 150,
        anc_read_length: Optional[int] = None,
        min_overlap: int = 12,
        min_jct_cov: int = 5,
        has_ancestor: bool = False,
        force: Optional[bool] = None,
    ):
        self.filtered_tnjc2s = filtered_tnjc2s
        self.output_dir = Path(output_dir)
        self.anc_output_dir = Path(anc_output_dir) if anc_output_dir else None
        self.read_length = read_length
        self.anc_read_length = anc_read_length if anc_read_length else read_length
        self.min_overlap = min_overlap
        self.min_jct_cov = min_jct_cov
        self.has_ancestor = has_ancestor

        # Input files are the BAM files from alignment step
        input_files = []
        for filtered_tnjc2 in filtered_tnjc2s:
            analysis_dir = output_dir / filtered_tnjc2.analysis_dir
            input_files.append(analysis_dir / "iso.sorted.bam")
            if has_ancestor:
                # Ancestor BAM is in ancestor folder (as iso.sorted.bam)
                if not self.anc_output_dir:
                    raise ValueError("anc_output_dir must be provided when has_ancestor=True")
                anc_analysis_dir = self.anc_output_dir / filtered_tnjc2.analysis_dir
                input_files.append(anc_analysis_dir / "iso.sorted.bam")

        super().__init__(
            output_dir=output_dir,
            input_files=input_files,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Analyze alignments for each candidate."""
        analyzed_records = []
        for filtered_tnjc2 in self.filtered_tnjc2s:
            analysis_dir = self.output_dir / filtered_tnjc2.analysis_dir

            # Get isolate junction coverage
            iso_bam = analysis_dir / "iso.sorted.bam"
            if not iso_bam.exists():
                info(f"Skipping {filtered_tnjc2.analysis_dir}: no iso.sorted.bam")
                continue

            iso_jc_cov = get_junction_coverage(iso_bam, self.read_length, min_overlap=self.min_overlap)

            # Classify isolate architecture
            iso_arch = classify_architecture(iso_jc_cov, self.min_jct_cov)

            # Get ancestor junction coverage if available
            # Ancestor BAM is stored in ancestor folder as iso.sorted.bam
            anc_jc_cov = None
            anc_arch = None
            if self.has_ancestor:
                if not self.anc_output_dir:
                    raise ValueError("anc_output_dir must be provided when has_ancestor=True")
                anc_analysis_dir = self.anc_output_dir / filtered_tnjc2.analysis_dir
                anc_bam = anc_analysis_dir / "iso.sorted.bam"

                if anc_bam.exists():
                    anc_jc_cov = get_junction_coverage(anc_bam, self.anc_read_length, min_overlap=self.min_overlap)
                    anc_arch = classify_architecture(anc_jc_cov, self.min_jct_cov)

            # Classify final event
            event, modifiers = classify_event(iso_arch, anc_arch, self.min_jct_cov)

            # Build AnalyzedTnJc2 record from candidate + new fields
            analyzed_tnjc2 = AnalyzedTnJc2.from_other(
                filtered_tnjc2,
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
            analyzed_records.append(analyzed_tnjc2)

        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)

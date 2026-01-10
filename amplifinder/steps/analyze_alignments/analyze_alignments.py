"""Step 12: Analyze read alignments to synthetic junctions."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import (
    RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts,
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.steps.analyze_alignments.parse_bam import get_junction_coverage


class AnalyzeTnJc2AlignmentsStep(RecordTypedDfStep[AnalyzedTnJc2]):
    """Analyze read alignments to get junction coverage.

    Analysis depends on run type:
    - has_ancestor=False: analyze isolate BAM only
    - has_ancestor=True: analyze both isolate and ancestor BAMs
    """

    def __init__(
        self,
        synjct_tnjc2s: RecordTypedDf[SynJctsTnJc2],
        output_dir: Path,
        anc_output_dir: Optional[Path] = None,
        iso_read_length: int = 150,
        anc_read_length: Optional[int] = None,
        min_overlap: int = 12,
        force: Optional[bool] = None,
    ):
        self.synjct_tnjc2s = synjct_tnjc2s
        self._iso_output_dir = Path(output_dir)
        self._anc_output_dir = Path(anc_output_dir) if anc_output_dir else None
        self.iso_read_length = iso_read_length
        self.anc_read_length = anc_read_length if anc_read_length else iso_read_length
        self.min_overlap = min_overlap
        self.has_ancestor = self._anc_output_dir is not None

        # Input files are the BAM files from alignment step
        input_files = []
        for synjct_tnjc2 in synjct_tnjc2s:
            iso_bam = synjct_tnjc2.bam_path(self._iso_output_dir)
            input_files.append(iso_bam)
            if self.has_ancestor:
                anc_bam = synjct_tnjc2.bam_path(self._anc_output_dir)
                input_files.append(anc_bam)

        super().__init__(
            output_dir=output_dir,
            input_files=input_files,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Analyze alignments for each candidate."""
        analyzed_records = []
        for synjct_tnjc2 in self.synjct_tnjc2s:
            # Get isolate junction coverage (required)
            jc_cov = self._get_cov(synjct_tnjc2, self._iso_output_dir, self.iso_read_length)

            # Get ancestor junction coverage if available (optional)
            jc_cov_anc = None
            if self.has_ancestor:
                jc_cov_anc = self._get_cov(synjct_tnjc2, self._anc_output_dir, self.anc_read_length)

            analyzed_tnjc2 = AnalyzedTnJc2.from_other(
                synjct_tnjc2,
                jc_cov=jc_cov,
                jc_cov_anc=jc_cov_anc,
            )
            analyzed_records.append(analyzed_tnjc2)

        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)

    def _get_cov(
        self, synjct_tnjc2: SynJctsTnJc2, base_dir: Path, read_length: int
    ) -> dict[JunctionType, JunctionReadCounts]:
        """Load junction coverage for iso/anc BAM and convert to dict."""
        bam = synjct_tnjc2.bam_path(base_dir)
        cov_map = get_junction_coverage(bam, read_length, min_overlap=self.min_overlap)
        return {JunctionType[name]: jc for name, jc in cov_map.items()}

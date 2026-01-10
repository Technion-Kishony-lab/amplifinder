"""Step 12: Analyze read alignments to synthetic junctions."""

import pysam

from pathlib import Path
from typing import Dict, Optional

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts
from amplifinder.steps.base import RecordTypedDfStep


def get_junction_coverage(
    bam_path: Path,
    read_length: int,
    read_length_tolerance: float = 0.1,
    min_overlap: int = 12,
) -> Dict[str, JunctionReadCounts]:
    """Parse BAM and get coverage for all 7 junction types.

    Args:
        bam_path: Path to sorted BAM file
        read_length: Read length for filtering
        read_length_tolerance: Tolerance for read length filtering (default 0.1 = 10%)
        min_overlap: Minimum overlap to count as spanning

    Returns:
        dict mapping reference name -> JunctionReadCounts
    """
    if not bam_path.exists():
        raise FileNotFoundError(f"BAM file not found: {bam_path}")

    read_length_factor = 1 + read_length_tolerance
    min_len = read_length / read_length_factor
    max_len = read_length * read_length_factor

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        ref_lengths = dict(zip(bam.references, bam.lengths))
        junction_points = {name: length // 2 for name, length in ref_lengths.items()}

        # [left, right, spanning] for each junction type
        counts = {name: JunctionReadCounts() for name in ref_lengths.keys()}

        for read in bam.fetch():
            if read.is_unmapped:
                continue

            ref_name = read.reference_name
            assert ref_name in junction_points, f"Reference name {ref_name} not found in junction points"
            length = read.query_alignment_length
            if not (min_len <= length <= max_len):
                continue

            start = read.reference_start
            end = read.reference_end
            junction_point = junction_points[ref_name]

            if end <= junction_point:
                counts[ref_name].left += 1      # left of the junction
            elif start > junction_point:
                counts[ref_name].right += 1     # right of the junction
            elif start <= junction_point - min_overlap and end >= junction_point + min_overlap:
                counts[ref_name].spanning += 1  # spanning the junction
            else:
                counts[ref_name].other += 1     # other (partial overlap)

    return counts



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

"""Step 12: Analyze read alignments to synthetic junctions."""
import pysam

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.visualization.plot_alignments import plot_alignment_coverage


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
            iso_bam = synjct_tnjc2.bam_path(self._iso_output_dir, is_ancestor=False)
            input_files.append(iso_bam)
            if self.has_ancestor:
                anc_bam = synjct_tnjc2.bam_path(self._anc_output_dir, is_ancestor=True)
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
            jc_cov, alignment_data, jct_lengths = self._get_cov(synjct_tnjc2, self._iso_output_dir, self.iso_read_length, is_ancestor=False)
            
            # Get ancestor junction coverage if available (optional)
            jc_cov_anc = None
            alignment_data_anc = None
            if self.has_ancestor:
                jc_cov_anc, alignment_data_anc, jct_lengths_anc = self._get_cov(synjct_tnjc2, self._anc_output_dir, self.anc_read_length, is_ancestor=True)
            
            # Generate combined visualization (isolate above, ancestor below x-axis)
            self._generate_alignment_plot(synjct_tnjc2, alignment_data, jct_lengths, 
                                        alignment_data_anc=alignment_data_anc, is_ancestor=False)

            analyzed_tnjc2 = AnalyzedTnJc2.from_other(
                synjct_tnjc2,
                jc_cov=jc_cov,
                jc_cov_anc=jc_cov_anc,
            )
            analyzed_records.append(analyzed_tnjc2)

        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)

    def _get_cov(
        self,
        synjct_tnjc2: SynJctsTnJc2,
        base_dir: Path,
        avg_read_length: int,
        is_ancestor: bool,
        read_length_tolerance: float = 0.1,
    ) -> tuple[dict[JunctionType, JunctionReadCounts], dict[JunctionType, list[tuple[int, int, str]]], dict[JunctionType, int]]:
        """Parse BAM and get coverage for all 7 junction types.

        Args:
            synjct_tnjc2: Synthetic junction record
            base_dir: Base directory containing BAM file
            avg_read_length: Read length for filtering
            is_ancestor: Whether processing ancestor (True) or isolate (False) data
            read_length_tolerance: Tolerance for read length filtering (default 0.1 = 10%)

        Returns:
            tuple of:
                - dict mapping JunctionType -> JunctionReadCounts
                - dict mapping JunctionType -> list of (start, end, read_type) tuples
                - dict mapping JunctionType -> junction length
        """
        bam_path = synjct_tnjc2.bam_path(base_dir, is_ancestor=is_ancestor)
        if not bam_path.exists():
            raise FileNotFoundError(f"BAM file not found: {bam_path}")

        read_length_factor = 1 + read_length_tolerance
        min_alignment_length = avg_read_length / read_length_factor
        max_alignment_length = avg_read_length * read_length_factor

        counts = {jt: JunctionReadCounts() for jt in JunctionType.sorted()}
        alignment_data = {jt: [] for jt in JunctionType.sorted()}
        jct_lengths = {}

        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            jct_lengths_raw = dict(zip(bam.references, bam.lengths))
            jct_middles = {name: length // 2 for name, length in jct_lengths_raw.items()}
            
            # Store lengths by JunctionType
            for jct_name, length in jct_lengths_raw.items():
                jct_type = JunctionType[jct_name]
                jct_lengths[jct_type] = length

            for read in bam.fetch():
                if read.is_unmapped:
                    continue

                jct_name = read.reference_name
                jct_type = JunctionType[jct_name]
                alignment_length = read.query_alignment_length
                if not (min_alignment_length <= alignment_length <= max_alignment_length):
                    continue

                start = read.reference_start
                end = read.reference_end
                assert start <= end
                junction_point = jct_middles[jct_name]

                # Determine read type and update counts
                if end <= junction_point:
                    read_type = 'left'
                    counts[jct_type].left += 1
                elif start > junction_point:
                    read_type = 'right'
                    counts[jct_type].right += 1
                elif start <= junction_point - self.min_overlap and end >= junction_point + self.min_overlap:
                    read_type = 'spanning'
                    counts[jct_type].spanning += 1
                else:
                    read_type = 'undetermined'
                    counts[jct_type].undetermined += 1

                # Store alignment data for plotting
                alignment_data[jct_type].append((start, end, read_type))

        return counts, alignment_data, jct_lengths

    def _generate_alignment_plot(
        self,
        synjct_tnjc2: SynJctsTnJc2,
        alignment_data: dict[JunctionType, list[tuple[int, int, str]]],
        jct_lengths: dict[JunctionType, int],
        alignment_data_anc: dict[JunctionType, list[tuple[int, int, str]]] | None = None,
        is_ancestor: bool = False,
    ) -> None:
        """Generate PNG plot showing read alignment coverage for all junction types.
        
        Args:
            synjct_tnjc2: Synthetic junction record
            alignment_data: Dict mapping JunctionType to list of (start, end, read_type) tuples for isolate
            jct_lengths: Dict mapping JunctionType to junction length
            alignment_data_anc: Optional dict for ancestor reads (plotted below x-axis)
            is_ancestor: Whether processing ancestor (True) or isolate (False) data
        """
        title = f'Read Alignment Coverage - {synjct_tnjc2.analysis_dir_name(is_ancestor=is_ancestor)}'
        base_dir = self._anc_output_dir if is_ancestor else self._iso_output_dir
        output_path = synjct_tnjc2.alignment_plot_path(base_dir, is_ancestor=is_ancestor)
        
        plot_alignment_coverage(
            alignment_data=alignment_data,
            jct_lengths=jct_lengths,
            title=title,
            output_path=output_path,
            alignment_data_anc=alignment_data_anc,
        )

"""Step 12: Analyze read alignments to synthetic junctions."""
import pysam
import numpy as np

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.visualization.plot_alignments import plot_alignment_coverage


def is_covered(cov: JunctionReadCounts, min_jct_cov: int, 
               jc_len: int, read_len: int, min_overlap_len: int, num_std: int = 3) -> Optional[bool]:
    """Determine if junction is covered based on spanning read statistics.
    
    Args:
        cov: Junction read counts (left, right, spanning)
        min_jct_cov: Minimum expected spanning reads threshold
        jc_len: Junction length
        read_len: Read length
        min_overlap_len: Minimum overlap length for spanning reads
        num_std: Number of standard deviations for coverage threshold (default 3)
    
    Returns:
        True if junction is covered, False if not covered, None if ambiguous
    """
    arm_len = jc_len // 2

    # num options of read-alignments for the left or right side
    num_options_side_aligned_reads = arm_len - read_len  

    # num options of read-alignments spanning the junction
    num_options_spanning_reads = read_len - 2 * min_overlap_len

    ratio_spanning_reads = num_options_spanning_reads / num_options_side_aligned_reads

    # if the junction connects a single-copy region with a multi-copy region, 
    # the junction, if exists, should be covered as expected based on the low-copy region
    min_num_reads_left_right = min(cov.left, cov.right)

    # expected number of spanning reads, and std err
    expected_num_spanning = min_num_reads_left_right * ratio_spanning_reads
    err_expected_num_spanning = np.sqrt(min_num_reads_left_right) * ratio_spanning_reads
    
    num_spanning_reads = cov.spanning
    err_num_spanning_reads = np.sqrt(num_spanning_reads)

    total_err = np.sqrt(err_expected_num_spanning**2 + err_num_spanning_reads**2)

    is_above_minimal_expected = num_spanning_reads >= expected_num_spanning - num_std * total_err
    is_below_min_jct_cov = num_spanning_reads <= min_jct_cov

    if is_above_minimal_expected and not is_below_min_jct_cov:
        return True
    if not is_above_minimal_expected and is_below_min_jct_cov:
        return False
    return None  # ambiguous


class AnalyzeTnJc2AlignmentsStep(RecordTypedDfStep[AnalyzedTnJc2]):
    """Analyze read alignments to get junction coverage."""

    def __init__(
        self,
        synjct_tnjc2s: RecordTypedDf[SynJctsTnJc2],
        output_dir: Path,
        anc_output_dir: Optional[Path] = None,
        iso_read_length: int = 150,
        anc_read_length: Optional[int] = None,
        min_overlap_len: int = 12,
        min_jct_cov: int = 10,
        force: Optional[bool] = None,
    ):
        self.synjct_tnjc2s = synjct_tnjc2s
        self._iso_output_dir = Path(output_dir)
        self._anc_output_dir = Path(anc_output_dir) if anc_output_dir else None
        self.iso_read_length = iso_read_length
        self.anc_read_length = anc_read_length if anc_read_length else iso_read_length
        self.min_overlap_len = min_overlap_len
        self.min_jct_cov = min_jct_cov
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

    def _calculate_jc_calls(
        self,
        jc_cov: dict[JunctionType, JunctionReadCounts],
        jct_lengths: dict[JunctionType, int],
        read_len: int,
    ) -> dict[JunctionType, Optional[bool]]:
        """Calculate junction coverage calls for all junction types.
        
        Args:
            jc_cov: Dict mapping JunctionType to JunctionReadCounts
            jct_lengths: Dict mapping JunctionType to junction length
            read_len: Read length
        
        Returns:
            Dict mapping JunctionType to coverage call (True/False/None)
        """
        return {
            jt: is_covered(
                jc_cov[jt],
                min_jct_cov=self.min_jct_cov,
                jc_len=jct_lengths[jt],
                read_len=read_len,
                min_overlap_len=self.min_overlap_len
            )
            for jt in JunctionType.sorted()
        }

    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Analyze alignments for each candidate.
        
        Returns:
            RecordTypedDf containing AnalyzedTnJc2 records with coverage data and calls
        """
        analyzed_records = []
        for synjct_tnjc2 in self.synjct_tnjc2s:
            # Get isolate junction coverage (required)
            jc_cov, alignment_data, jct_lengths = self._get_cov(synjct_tnjc2, self._iso_output_dir, self.iso_read_length, is_ancestor=False)
            jc_calls = self._calculate_jc_calls(jc_cov, jct_lengths, self.iso_read_length)

            # Get ancestor junction coverage if available (optional)
            jc_cov_anc = None
            alignment_data_anc = None
            jc_calls_anc = None
            if self.has_ancestor:
                jc_cov_anc, alignment_data_anc, jct_lengths_anc = self._get_cov(synjct_tnjc2, self._anc_output_dir, self.anc_read_length, is_ancestor=True)
                jc_calls_anc = self._calculate_jc_calls(jc_cov_anc, jct_lengths_anc, self.anc_read_length)

            # Generate combined visualization (isolate above, ancestor below x-axis)
            self._generate_alignment_plot(synjct_tnjc2, alignment_data, jct_lengths, 
                                          alignment_data_anc=alignment_data_anc, 
                                          jc_calls=jc_calls, jc_calls_anc=jc_calls_anc,
                                          is_ancestor=False)

            analyzed_tnjc2 = AnalyzedTnJc2.from_other(
                synjct_tnjc2,
                jc_cov=jc_cov,
                jc_cov_anc=jc_cov_anc,
                jc_calls=jc_calls,
                jc_calls_anc=jc_calls_anc,
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
            
            # Store lengths by JunctionType
            jct_lengths_raw = dict(zip(bam.references, bam.lengths))
            jct_lengths = {JunctionType[jct_name]: length for jct_name, length in jct_lengths_raw.items()}

            # Count read-alignments of each type at each junction
            for read in bam.fetch():
                if read.is_unmapped:
                    continue

                jct_name = read.reference_name
                jct_type = JunctionType[jct_name]
                alignment_length = read.query_alignment_length
                if not (min_alignment_length <= alignment_length <= max_alignment_length):
                    continue

                start, end = read.reference_start, read.reference_end
                assert start <= end

                read_type = counts[jct_type].add_read(start, end, junction_point=jct_lengths[jct_type] // 2, 
                                                      min_overlap_len=self.min_overlap_len)

                # Store alignment data (for plotting, etc)
                alignment_data[jct_type].append((start, end, read_type))

        return counts, alignment_data, jct_lengths

    def _generate_alignment_plot(
        self,
        synjct_tnjc2: SynJctsTnJc2,
        alignment_data: dict[JunctionType, list[tuple[int, int, str]]],
        jct_lengths: dict[JunctionType, int],
        alignment_data_anc: dict[JunctionType, list[tuple[int, int, str]]] | None = None,
        jc_calls: dict[JunctionType, Optional[bool]] | None = None,
        jc_calls_anc: dict[JunctionType, Optional[bool]] | None = None,
        is_ancestor: bool = False,
    ) -> None:
        """Generate PNG plot showing read alignment coverage for all junction types.
        
        Args:
            synjct_tnjc2: Synthetic junction record
            alignment_data: Dict mapping JunctionType to list of (start, end, read_type) tuples for isolate
            jct_lengths: Dict mapping JunctionType to junction length
            alignment_data_anc: Optional dict for ancestor reads (plotted below x-axis)
            jc_calls: Dict mapping JunctionType to coverage call (True/False/None) for isolate
            jc_calls_anc: Optional dict mapping JunctionType to coverage call for ancestor
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
            jc_calls=jc_calls,
            jc_calls_anc=jc_calls_anc,
        )

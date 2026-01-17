"""Step 12: Analyze read alignments to synthetic junctions."""
import pysam
import numpy as np

from pathlib import Path
from typing import Optional, Dict

from amplifinder.optional_deps import plt

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts, Genome, Side
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.steps.jct_coverage.read_type import get_jct_read_counts
from amplifinder.tools.breseq import load_breseq_coverage
from amplifinder.visualization.plot_coverage import plot_amplicon_coverage
from amplifinder.visualization.plot_alignments import plot_junctions_coverage


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
        genome: Genome,
        iso_breseq_path: Path,
        anc_output_dir: Optional[Path] = None,
        anc_breseq_path: Optional[Path] = None,
        iso_read_length: int = 150,
        anc_read_length: Optional[int] = None,
        min_overlap_len: int = 12,
        min_jct_cov: int = 10,
        create_plots: bool = True,
        force: Optional[bool] = None,
    ):
        self.synjct_tnjc2s = synjct_tnjc2s
        self._iso_output_dir = Path(output_dir)
        self._anc_output_dir = Path(anc_output_dir) if anc_output_dir else None
        self.genome = genome
        self.iso_breseq_path = Path(iso_breseq_path)
        self.anc_breseq_path = Path(anc_breseq_path) if anc_breseq_path else None
        self.iso_read_length = iso_read_length
        self.anc_read_length = anc_read_length if anc_read_length else iso_read_length
        self.min_overlap_len = min_overlap_len
        self.min_jct_cov = min_jct_cov
        self.create_plots = create_plots
        self.has_ancestor = self._anc_output_dir is not None

        # Validate matplotlib availability if plotting requested
        if self.create_plots and plt is None:
            raise ImportError(
                "Plotting is enabled (create_plots=True) but matplotlib is not installed. "
                "Install with: pip install matplotlib, or disable plotting with --no-create-plots"
            )

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

    def _load_coverage_for_plotting(self) -> tuple[Optional[Dict[str, np.ndarray]], Optional[Dict[str, np.ndarray]]]:
        """Load genome coverage data for plotting (if enabled).

        Returns:
            Tuple of (iso_coverage, anc_coverage) dicts, or (None, None) if plotting disabled
        """
        iso_scafs_to_covs = load_breseq_coverage(self.iso_breseq_path)
        anc_scafs_to_covs = None
        if self.has_ancestor and self.anc_breseq_path:
            anc_scafs_to_covs = load_breseq_coverage(self.anc_breseq_path)

        return iso_scafs_to_covs, anc_scafs_to_covs


    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Analyze alignments for each candidate.

        Returns:
            RecordTypedDf containing AnalyzedTnJc2 records with coverage data and calls
        """
        # Load coverage data once for all candidates (if plotting enabled)
        if self.create_plots:
            iso_scafs_to_covs, anc_scafs_to_covs = self._load_coverage_for_plotting()
            print(f'Creating coverage plots (n={len(self.synjct_tnjc2s)}) ', end='', flush=True)

        analyzed_records = []
        for synjct_tnjc2 in self.synjct_tnjc2s:
            # Get isolate junction coverage (required)
            jc_covs, alignment_data, jct_lengths = self._get_cov(
                synjct_tnjc2, self._iso_output_dir, self.iso_read_length, is_ancestor=False
            )
            jc_calls = self._calculate_jc_calls(jc_covs, jct_lengths, self.iso_read_length)

            # Get ancestor junction coverage if available (optional)
            jc_covs_anc = None
            alignment_data_anc = None
            jc_calls_anc = None
            if self.has_ancestor:
                jc_covs_anc, alignment_data_anc, jct_lengths_anc = self._get_cov(
                    synjct_tnjc2, self._anc_output_dir, self.anc_read_length, is_ancestor=True
                )
                jc_calls_anc = self._calculate_jc_calls(jc_covs_anc, jct_lengths_anc, self.anc_read_length)

            analyzed_tnjc2 = AnalyzedTnJc2.from_other(
                synjct_tnjc2,
                jc_covs=jc_covs,
                jc_covs_anc=jc_covs_anc,
                jc_calls=jc_calls,
                jc_calls_anc=jc_calls_anc,
            )
            analyzed_records.append(analyzed_tnjc2)

            if self.create_plots:
                plot_junctions_coverage(
                    jc_lengths=jct_lengths,
                    alignment_data=alignment_data, alignment_data_anc=alignment_data_anc,
                    jc_covs=jc_covs, jc_covs_anc=jc_covs_anc,
                    jc_calls=jc_calls, jc_calls_anc=jc_calls_anc,
                    title=f'Jcts coverage - {synjct_tnjc2.analysis_dir_name(is_ancestor=False)}',
                    output_path=synjct_tnjc2.analysis_dir_path(self._iso_output_dir) / "jct_coverages.png",
                    min_overlap_len=self.min_overlap_len,
                    read_len=self.iso_read_length,
                )
                plot_amplicon_coverage(
                    tnjc2=analyzed_tnjc2,
                    iso_scafs_to_covs=iso_scafs_to_covs,
                    anc_scafs_to_covs=anc_scafs_to_covs,
                    output_path=analyzed_tnjc2.analysis_dir_path(self._iso_output_dir) / "amplicon_coverage.png",
                )
                print('.', end='', flush=True)

        if self.create_plots:
            print('\n', flush=True)

        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)

    def _get_cov(
        self,
        synjct_tnjc2: SynJctsTnJc2,
        base_dir: Path,
        avg_read_length: int,
        is_ancestor: bool,
        read_length_tolerance: float = 0.1,
        max_dist_from_junction: int = 10,
        max_nm_score: int = 3,
        min_as_score: int = -25
    ) -> tuple[
        dict[JunctionType, JunctionReadCounts],
        dict[JunctionType, list[tuple[int, int, str]]],
        dict[JunctionType, int]
    ]:
        """Parse BAM and get coverage for all 7 junction types.

        Args:
            synjct_tnjc2: Synthetic junction record
            base_dir: Base directory containing BAM file
            avg_read_length: Read length for filtering
            is_ancestor: Whether processing ancestor (True) or isolate (False) data
            read_length_tolerance: Tolerance for read length filtering (default 0.1 = 10%)
            max_dist_from_junction: Maximum distance from junction for read classification
            max_nm_score: Maximum NM score threshold (default 3)
            min_as_score: Minimum AS score threshold (default -25)

        Returns:
            tuple of:
                - dict mapping JunctionType -> JunctionReadCounts
                - dict mapping JunctionType -> list of (start, end, read_type) tuples
                - dict mapping JunctionType -> junction length
        """
        bam_path = synjct_tnjc2.bam_path(base_dir, is_ancestor=is_ancestor)
        return get_jct_read_counts(bam_path, avg_read_length, self.min_overlap_len, 
                                   alignment_length_tolerance=read_length_tolerance,
                                   max_dist_from_junction=max_dist_from_junction,
                                   max_nm_score=max_nm_score,
                                   min_as_score=min_as_score)

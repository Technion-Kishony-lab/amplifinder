"""Step 12: Analyze read alignments to synthetic junctions."""
import numpy as np

from pathlib import Path
from typing import Optional, Dict

from amplifinder.config import AlignmentAnalysisParams

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.steps.jct_coverage.read_type import get_jct_read_counts


def get_jct_read_counts_by_tnjc2(
    synjct_tnjc2: SynJctsTnJc2,
    base_dir: Path,
    is_ancestor: bool,
    avg_read_length: int,
    jct_align_params: AlignmentAnalysisParams,
) -> tuple[
    dict[JunctionType, JunctionReadCounts],
    dict[JunctionType, list[tuple[int, int, str]]],
    dict[JunctionType, int]
]:
    """A wrapper around get_jct_read_counts to get the read counts based on a SynJctsTnJc2."""
    bam_path = synjct_tnjc2.bam_path(base_dir, is_ancestor=is_ancestor)
    return get_jct_read_counts(
        bam_path,
        avg_read_length,
        jct_align_params.min_overlap_len,
        alignment_length_tolerance=jct_align_params.read_length_tolerance,
        max_dist_from_junction=jct_align_params.max_dist_from_junction,
        max_nm_score=jct_align_params.max_nm_score,
        min_as_score=jct_align_params.min_as_score,
    )


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
    """Analyze isolate read alignments to get junction coverage."""

    IS_ANCESTOR = False
    COV_FIELD = "jc_covs"
    CALLS_FIELD = "jc_calls"

    def __init__(
        self,
        tnjc2s: RecordTypedDf[SynJctsTnJc2],
        output_dir: Path,
        read_length: int = 150,
        jct_align_params: AlignmentAnalysisParams = None,
        min_jct_cov: int = 10,
        force: Optional[bool] = None,
    ):
        self.tnjc2s = tnjc2s
        self.read_length = read_length
        self.jct_align_params = jct_align_params or AlignmentAnalysisParams()
        self.min_jct_cov = min_jct_cov
        self._output_dir = Path(output_dir)

        # compatibility handle for tests/consumers
        self.synjct_tnjc2s = tnjc2s

        input_files = [tnjc2.bam_path(self._output_dir, is_ancestor=self.IS_ANCESTOR) for tnjc2 in tnjc2s]

        super().__init__(
            output_dir=output_dir,
            input_files=input_files,
            force=force,
        )

    def _calculate_output(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Analyze alignments for each candidate."""
        analyzed_records = []
        for tnjc2 in self.tnjc2s:
            jc_covs, _, jc_lengths = get_jct_read_counts_by_tnjc2(
                synjct_tnjc2=tnjc2,
                base_dir=self._output_dir,
                is_ancestor=self.IS_ANCESTOR,
                avg_read_length=self.read_length,
                jct_align_params=self.jct_align_params,
            )
            jc_calls = {
                jt: is_covered(
                    jc_covs[jt],
                    min_jct_cov=self.min_jct_cov,
                    jc_len=jc_lengths[jt],
                    read_len=self.read_length,
                    min_overlap_len=self.jct_align_params.min_overlap_len
                )
                for jt in JunctionType.sorted()
            }

            analyzed_records.append(AnalyzedTnJc2.from_other(
                tnjc2,
                **{self.COV_FIELD: jc_covs, self.CALLS_FIELD: jc_calls},
            ))

        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)


class AncAnalyzeTnJc2AlignmentsStep(AnalyzeTnJc2AlignmentsStep):
    """Analyze ancestor read alignments."""

    IS_ANCESTOR = True
    COV_FIELD = "jc_covs_anc"
    CALLS_FIELD = "jc_calls_anc"

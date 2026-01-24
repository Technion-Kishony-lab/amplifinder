"""Step 12: Analyze read alignments to synthetic junctions."""
import numpy as np

from pathlib import Path
from typing import Optional

from amplifinder.config import AlignmentClassifyParams, JcCallParams, AlignmentFilterParams
from amplifinder.env import DEBUG
from amplifinder.data_types import (
    RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts, ReadType
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.steps.jct_coverage.alignment_data import AlignmentData
from amplifinder.steps.jct_coverage.classify_alignments import get_jct_read_counts
from amplifinder.steps.jct_coverage.export_bam_indices import write_junction_read_bam_indices


def get_jct_read_counts_by_tnjc2(
    synjct_tnjc2: SynJctsTnJc2,
    base_dir: Path,
    is_ancestor: bool,
    arm_len: int,
    avg_read_length: int,
    alignment_classify_params: AlignmentClassifyParams,
    alignment_filter_params: AlignmentFilterParams,
) -> tuple[
    dict[JunctionType, JunctionReadCounts],
    dict[JunctionType, dict[ReadType, dict[str, AlignmentData]]]
]:
    """A wrapper around get_jct_read_counts to get the read counts based on a SynJctsTnJc2."""
    bam_path = synjct_tnjc2.bam_path(base_dir, is_ancestor=is_ancestor)
    return get_jct_read_counts(bam_path, arm_len, avg_read_length, alignment_classify_params, alignment_filter_params)


def is_covered(cov: JunctionReadCounts, jc_call_params: JcCallParams,
               arm_len: int, read_len: int, min_overlap_len: int) -> Optional[bool]:
    """Determine if junction is covered based on spanning read statistics.

    Args:
        cov: Junction read counts (left, right, spanning)
        jc_call_params: Junction calling parameters. See JcCallParams for details.
        arm_len: Junction arm length
        read_len: Read length
        min_overlap_len: Minimum overlap length for spanning reads

    Returns:
        True if junction is covered, False if not covered, None if ambiguous
    """

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

    is_above_minimal_expected = \
        num_spanning_reads >= expected_num_spanning \
        - jc_call_params.pos_threshold_in_num_std_below_expected * total_err
    is_close_to_zero = num_spanning_reads <= jc_call_params.neg_threshold_abs \
        or num_spanning_reads <= expected_num_spanning * jc_call_params.neg_threshold_rel

    if is_above_minimal_expected and not is_close_to_zero:
        return True
    if not is_above_minimal_expected and is_close_to_zero:
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
        arm_len: int,
        read_length: int = 150,
        alignment_classify_params: AlignmentClassifyParams = None,
        alignment_filter_params: AlignmentFilterParams = None,
        jc_call_params: JcCallParams = None,
        force: Optional[bool] = None,
    ):
        self.tnjc2s = tnjc2s
        self.arm_len = arm_len
        self.read_length = read_length
        self.alignment_classify_params = alignment_classify_params or AlignmentClassifyParams()
        self.alignment_filter_params = alignment_filter_params or AlignmentFilterParams()
        self.jc_call_params = jc_call_params or JcCallParams()
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
            jc_covs, alignment_data = get_jct_read_counts_by_tnjc2(
                synjct_tnjc2=tnjc2,
                base_dir=self._output_dir,
                is_ancestor=self.IS_ANCESTOR,
                arm_len=self.arm_len,
                avg_read_length=self.read_length,
                alignment_classify_params=self.alignment_classify_params,
                alignment_filter_params=self.alignment_filter_params,
            )

            if DEBUG:
                write_junction_read_bam_indices(
                    alignment_data=alignment_data,
                    output_dir=tnjc2.analysis_dir_path(self._output_dir, is_ancestor=self.IS_ANCESTOR),
                )

            jc_calls = {
                jt: is_covered(
                    jc_covs[jt],
                    jc_call_params=self.jc_call_params,
                    arm_len=self.arm_len,
                    read_len=self.read_length,
                    min_overlap_len=self.alignment_classify_params.min_overlap_len
                )
                for jt in JunctionType
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

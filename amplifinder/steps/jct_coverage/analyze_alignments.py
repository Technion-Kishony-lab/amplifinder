"""Step 12: Analyze read alignments to synthetic junctions."""
import numpy as np

from pathlib import Path
from typing import Optional

from amplifinder.config import AlignmentClassifyParams, JcCallParams, AlignmentFilterParams
from amplifinder.env import DEBUG
from amplifinder.data_types import (
    RecordTypedDf, SynJctsTnJc2, AnalyzedTnJc2, JunctionType, JunctionReadCounts
)
from amplifinder.logger import logger
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.steps.jct_coverage.alignment_data import AlignmentData
from amplifinder.steps.jct_coverage.classify_alignments import get_expected_counts, get_jct_read_counts
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
    dict[JunctionType, list[AlignmentData]]
]:
    """A wrapper around get_jct_read_counts to get the read counts based on a SynJctsTnJc2."""
    bam_path = synjct_tnjc2.bam_path(base_dir, is_ancestor=is_ancestor)
    return get_jct_read_counts(bam_path, arm_len, avg_read_length, alignment_classify_params, alignment_filter_params)


def is_covered(cov: JunctionReadCounts, jc_call_params: JcCallParams,
               arm_len: int, read_len: int, alignment_classify_params: AlignmentClassifyParams) -> Optional[bool]:
    """Determine if junction is covered based on spanning read statistics.

    Args:
        cov: Junction read counts (left, right, spanning)
        jc_call_params: Junction calling parameters. See JcCallParams for details.
        arm_len: Junction arm length
        read_len: Read length
        alignment_classify_params: Alignment classification parameters

    Returns:
        True if junction is covered, False if not covered, None if ambiguous
    """

    expected_counts = get_expected_counts(read_len, arm_len, alignment_classify_params)
    assert expected_counts.left == expected_counts.right
    ratio_spanning_reads = expected_counts.spanning / expected_counts.left

    # if the junction connects a single-copy region with a multi-copy region,
    # the junction, if exists, should be covered as expected based on the low-copy region
    min_num_reads_left_right = min(cov.left, cov.right)

    # expected number of spanning reads, and std err
    expected_num_spanning = min_num_reads_left_right * ratio_spanning_reads
    err_expected_num_spanning = np.sqrt(min_num_reads_left_right) * ratio_spanning_reads

    num_spanning_reads = cov.spanning
    err_num_spanning_reads = np.sqrt(num_spanning_reads)

    total_err = np.sqrt(err_expected_num_spanning**2 + err_num_spanning_reads**2)
    num_std = jc_call_params.pos_threshold_in_num_std_below_expected

    is_above_minimal_expected = \
        num_spanning_reads * (1 + jc_call_params.pos_threshold_rel) \
        >= expected_num_spanning - num_std * total_err

    is_close_to_zero = num_spanning_reads <= jc_call_params.neg_threshold_abs \
        or num_spanning_reads <= expected_num_spanning * jc_call_params.neg_threshold_rel

    if is_above_minimal_expected and not is_close_to_zero:
        return True
    if not is_above_minimal_expected and is_close_to_zero:
        return False
    return None  # ambiguous


def print_jc_read_counts_and_calls(
    jc_covs: dict[JunctionType, JunctionReadCounts],
    jc_calls: dict[JunctionType, Optional[bool]]
) -> None:
    """Print junction coverage and calls."""
    # Print header (two rows)
    logger.info(f"{'junction':<23} {'left':>6} {'':>6} {'left':>6} "
          f"{'':>6} {'':>6} {'right':>6} {'':>6} "
          f"{'right':>6} {'':>6}", timestamp=False)
    logger.info(f"{'type':<23} {'far':>6} {'left':>6} {'marg':>6} "
          f"{'span':>6} {'pair':>6} {'marg':>6} {'right':>6} "
          f"{'far':>6} {'call':>6}", timestamp=False)
    logger.info("-" * 86, timestamp=False)

    # Print rows
    for jt in JunctionType:
        cov = jc_covs[jt]
        call_str = str(jc_calls[jt]) if jc_calls[jt] is not None else "None"
        logger.info(f"{jt.name:<23} {cov.left_far:>6} {cov.left:>6} {cov.left_marginal:>6} "
              f"{cov.spanning:>6} {cov.paired:>6} {cov.right_marginal:>6} {cov.right:>6} "
              f"{cov.right_far:>6} {call_str:>6}", timestamp=False)
    logger.info("", timestamp=False)


def warn_if_paired_and_not_spanning(jc_covs: dict[JunctionType, JunctionReadCounts]) -> None:
    """Warn if paired reads are present and no spanning reads are present."""
    for jt in JunctionType:
        if jc_covs[jt].paired > 0 and jc_covs[jt].spanning == 0:
            logger.warning(f"Junction {jt.name} has {jc_covs[jt].paired} paired-end reads "
                    "transversing the junction but no spanning reads.")


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

        # Cache for alignment data (populated during _calculate_output)
        self.alignment_data_cache: dict[str, dict[JunctionType, list[AlignmentData]]] = {}

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
            logger.info(f"\nCASE: '{tnjc2.analysis_dir}'")
            jc_covs, alignment_data = get_jct_read_counts_by_tnjc2(
                synjct_tnjc2=tnjc2,
                base_dir=self._output_dir,
                is_ancestor=self.IS_ANCESTOR,
                arm_len=self.arm_len,
                avg_read_length=self.read_length,
                alignment_classify_params=self.alignment_classify_params,
                alignment_filter_params=self.alignment_filter_params,
            )

            # Cache alignment data for later use (e.g., plotting)
            cache_key = tnjc2.analysis_dir
            self.alignment_data_cache[cache_key] = alignment_data

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
                    alignment_classify_params=self.alignment_classify_params
                )
                for jt in JunctionType
            }

            warn_if_paired_and_not_spanning(jc_covs)

            print_jc_read_counts_and_calls(jc_covs, jc_calls)

            analyzed_records.append(AnalyzedTnJc2.from_other(
                tnjc2,
                **{self.COV_FIELD: jc_covs, self.CALLS_FIELD: jc_calls},
            ))

        return RecordTypedDf.from_records(analyzed_records, AnalyzedTnJc2)

    def run_with_cache(self) -> tuple[RecordTypedDf[AnalyzedTnJc2], dict[str, dict[JunctionType, list[AlignmentData]]]]:
        """Run the step and return both results and alignment data cache."""
        result = self.run()
        return result, self.alignment_data_cache


class AncAnalyzeTnJc2AlignmentsStep(AnalyzeTnJc2AlignmentsStep):
    """Analyze ancestor read alignments."""

    IS_ANCESTOR = True
    COV_FIELD = "jc_covs_anc"
    CALLS_FIELD = "jc_calls_anc"

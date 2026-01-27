"""Step: Plot junction and amplicon coverage after classification."""
from pathlib import Path
from typing import Optional

from amplifinder.config import AlignmentClassifyParams, AlignmentFilterParams
from amplifinder.logger import logger
from amplifinder.optional_deps import plt
from amplifinder.data_types import RecordTypedDf, ClassifiedTnJc2, JunctionType
from amplifinder.steps.base import Step
from amplifinder.steps.read_length import ReadLengths
from amplifinder.steps.jct_coverage.alignment_data import AlignmentData
from amplifinder.tools.breseq import load_breseq_coverage
from amplifinder.visualization.plot_jc_alignments import plot_jc_alignments
from amplifinder.visualization.plot_amp_coverage import plot_amplicon_coverage


class PlotTnJc2CoverageStep(Step):
    """Generate junction and amplicon coverage plots (post-classification)."""

    def __init__(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        output_dir: Path,
        iso_breseq_path: Path,
        read_lengths: ReadLengths,
        anc_output_dir: Optional[Path] = None,
        anc_breseq_path: Optional[Path] = None,
        alignment_classify_params: AlignmentClassifyParams = None,
        alignment_filter_params: AlignmentFilterParams = None,
        iso_alignments: Optional[dict[str, dict[JunctionType, list[AlignmentData]]]] = None,
        anc_alignments: Optional[dict[str, dict[JunctionType, list[AlignmentData]]]] = None,
        force: Optional[bool] = None,
    ):
        self.classified_tnjc2s = classified_tnjc2s
        self._iso_output_dir = Path(output_dir)
        self._anc_output_dir = Path(anc_output_dir) if anc_output_dir else None
        self.iso_breseq_path = Path(iso_breseq_path)
        self.anc_breseq_path = Path(anc_breseq_path) if anc_breseq_path else None
        self.read_lengths = read_lengths
        self.alignment_classify_params = alignment_classify_params or AlignmentClassifyParams()
        self.alignment_filter_params = alignment_filter_params or AlignmentFilterParams()
        self.iso_alignment_data_cache = iso_alignments or {}
        self.anc_alignment_data_cache = anc_alignments or {}
        self.has_ancestor = self._anc_output_dir is not None

        if plt is None:
            raise ImportError(
                "Plotting is enabled (create_plots=True) but matplotlib is not installed. "
                "Install with: pip install matplotlib, or disable plotting with --no-create-plots"
            )

        # input files include:
        # 1. bam files
        input_files = [tnjc2.bam_path(self._iso_output_dir, is_ancestor=False) for tnjc2 in classified_tnjc2s]
        if self.has_ancestor:
            input_files += [tnjc2.bam_path(self._anc_output_dir, is_ancestor=True) for tnjc2 in classified_tnjc2s]

        # 2. breseq coverage dirs
        input_files.append(self.iso_breseq_path)
        if self.has_ancestor and self.anc_breseq_path:
            input_files.append(self.anc_breseq_path)

        # output artifact files (plots to generate)
        artifact_files = []
        for tnjc2 in classified_tnjc2s:
            artifact_files.append(tnjc2.analysis_dir_path(self._iso_output_dir) / "jct_coverages.png")
            artifact_files.append(tnjc2.analysis_dir_path(self._iso_output_dir) / "amplicon_coverage.png")

        super().__init__(input_files=input_files, artifact_files=artifact_files, force=force)

    def _load_coverage_for_plotting(self):
        iso_scafs_to_covs = load_breseq_coverage(self.iso_breseq_path)
        anc_scafs_to_covs = None
        if self.has_ancestor and self.anc_breseq_path:
            anc_scafs_to_covs = load_breseq_coverage(self.anc_breseq_path)
        return iso_scafs_to_covs, anc_scafs_to_covs

    def _generate_artifacts(self) -> None:
        iso_scafs_to_covs, anc_scafs_to_covs = self._load_coverage_for_plotting()
        logger.print_progress(f'Creating coverage plots (n={len(self.classified_tnjc2s)}) ', end='')

        for tnjc2 in self.classified_tnjc2s:
            jct_cov_path = tnjc2.analysis_dir_path(self._iso_output_dir) / "jct_coverages.png"
            amp_cov_path = tnjc2.analysis_dir_path(self._iso_output_dir) / "amplicon_coverage.png"

            # Skip if both plots already exist
            if jct_cov_path.exists() and amp_cov_path.exists():
                logger.print_progress('-', end='')
                continue

            # Get cached alignment data or re-read BAM if not cached
            cache_key = tnjc2.analysis_dir
            jc_to_alignments = self.iso_alignment_data_cache[cache_key]

            jc_covs = tnjc2.jc_covs
            jc_calls = tnjc2.jc_calls

            jc_covs_anc = None
            jc_to_alignments_anc = None
            jc_calls_anc = None
            if self.has_ancestor:
                jc_to_alignments_anc = self.anc_alignment_data_cache[cache_key]
                jc_covs_anc = tnjc2.jc_covs_anc
                jc_calls_anc = tnjc2.jc_calls_anc

            if not jct_cov_path.exists():
                plot_jc_alignments(
                    jc_to_alignments=jc_to_alignments,
                    jc_to_alignments_anc=jc_to_alignments_anc,
                    jc_covs=jc_covs,
                    jc_covs_anc=jc_covs_anc,
                    jc_calls=jc_calls,
                    jc_calls_anc=jc_calls_anc,
                    jc_arm_len_iso=self.read_lengths.jc_arm_len_iso,
                    jc_arm_len_anc=self.read_lengths.jc_arm_len_anc,
                    read_len_iso=self.read_lengths.read_len_iso,
                    read_len_anc=self.read_lengths.read_len_anc,
                    title=f'Jcts coverage - {tnjc2.analysis_dir_name(is_ancestor=False)}',
                    output_path=jct_cov_path,
                    alignment_classify_params=self.alignment_classify_params,
                )

            if not amp_cov_path.exists():
                plot_amplicon_coverage(
                    tnjc2=tnjc2,
                    iso_scafs_to_covs=iso_scafs_to_covs,
                    anc_scafs_to_covs=anc_scafs_to_covs,
                    output_path=amp_cov_path,
                )

            logger.print_progress('.', end='')

        logger.print_progress('')  # Newline

    def _artifact_labels(self) -> list[str]:
        """Summarize outputs as count."""
        n = len(self.classified_tnjc2s)
        return [f"{n} amplicon_coverage.png and {n} jct_coverages.png"]

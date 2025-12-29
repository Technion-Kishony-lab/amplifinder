"""Step 7: Calculate amplicon coverage."""

from pathlib import Path
from typing import Optional

import numpy as np

from amplifinder.data_types import RecordTypedDf, RawTnJc2, CoveredTnJc2, Genome, AverageMethod
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.tools.breseq import load_breseq_coverage
from amplifinder.utils import start_timer, end_timer

from .coverage import get_scaffold_coverage, calc_coverage_stats, calc_scaffold_coverages_and_stats, mean_positive
from .statistics import calc_distribution_mode


class CalcTnJc2AmpliconCoverageStep(RecordTypedDfStep[CoveredTnJc2]):
    """Calculate amplicon coverage for TnJc2 candidates.

    Coverage calculation depends on run type:
    - anc_breseq_path=None: raw coverage only
    - anc_breseq_path=set: normalized coverage (iso/anc ratio)
    """

    def __init__(
        self,
        raw_tnjc2s: RecordTypedDf[RawTnJc2],
        genome: Genome,
        output_dir: Path,
        ref_name: str,
        iso_breseq_path: Path,
        iso_name: str,
        anc_breseq_path: Optional[Path] = None,
        anc_name: Optional[str] = None,
        ncp_min: float = 0.1,
        ncp_max: float = 1000.0,
        ncp_n: int = 150,
        average_method: AverageMethod = AverageMethod.MEDIAN,
        force: Optional[bool] = None,
    ):
        self.raw_tnjc2s = raw_tnjc2s
        self.genome = genome
        self.ref_name = ref_name
        self.iso_breseq_path = Path(iso_breseq_path)
        self.iso_name = iso_name
        self.anc_breseq_path = Path(anc_breseq_path) if anc_breseq_path else None
        self.anc_name = anc_name
        self.ncp_min = ncp_min
        self.ncp_max = ncp_max
        self.ncp_n = ncp_n
        self.average_method = average_method

        input_files = [iso_breseq_path]
        if anc_breseq_path:
            input_files.append(anc_breseq_path)

        super().__init__(
            output_dir=output_dir,
            input_files=input_files,
            force=force,
        )

    @property
    def has_ancestor(self) -> bool:
        """True if ancestor comparison (iso/anc ratio) should be performed."""
        return self.anc_breseq_path is not None

    def _get_profiler_functions(self) -> list:
        """Functions to profile when profiling is enabled."""
        return [
            get_scaffold_coverage,
            calc_coverage_stats,
            calc_scaffold_coverages_and_stats,
            calc_distribution_mode,
            self._calc_candidate_coverage,
        ]

    def _calculate_output(self) -> RecordTypedDf[CoveredTnJc2]:
        """Calculate coverage for each TnJc2 candidate."""

        # Pre-calculate scaffold coverage arrays and medians for all unique scaffolds
        unique_scaffolds = list(set(raw_tnjc2.scaf for raw_tnjc2 in self.raw_tnjc2s))

        # Load isolate coverage
        with self.print_timer(f"loading iso coverage from {self.iso_breseq_path} ..."):
            iso_cov = load_breseq_coverage(self.iso_breseq_path, self.genome.name)
        with self.print_timer(f"calculating iso scaffold stats ({len(unique_scaffolds)} scaffolds) ..."):
            iso_scaf_covs, iso_scaf_stats = calc_scaffold_coverages_and_stats(iso_cov, unique_scaffolds, self.genome, self.average_method)

        # Load ancestor coverage if provided
        if self.has_ancestor:
            with self.print_timer(f"loading anc coverage from {self.anc_breseq_path} ..."):
                anc_cov = load_breseq_coverage(self.anc_breseq_path, self.genome.name)
            with self.print_timer(f"calculating anc scaffold stats ({len(unique_scaffolds)} scaffolds) ..."):
                anc_scaf_covs, anc_scaf_stats = calc_scaffold_coverages_and_stats(anc_cov, unique_scaffolds, self.genome, self.average_method)
        else:
            anc_cov = None
            anc_scaf_covs = {}
            anc_scaf_stats = {}

        # Process each candidate
        with self.print_timer("Calculating coverage for each candidate..."):
            covered_records = []
            for raw_tnjc2 in self.raw_tnjc2s:
                covered = self._calc_candidate_coverage(
                    raw_tnjc2,
                    iso_scaf_covs[raw_tnjc2.scaf],
                    iso_scaf_stats[raw_tnjc2.scaf],
                    anc_scaf_covs[raw_tnjc2.scaf] if self.has_ancestor else None,
                    anc_scaf_stats[raw_tnjc2.scaf] if self.has_ancestor else None,
                )
                covered_records.append(covered)
                self.print(".", end="")

        df = RecordTypedDf.from_records(covered_records, CoveredTnJc2)
        return df

    def _calc_candidate_coverage(
        self,
        raw_tnjc2: RawTnJc2,
        iso_scaf_cov: np.ndarray,
        iso_scaf_stat: float,
        anc_scaf_cov: Optional[np.ndarray],
        anc_scaf_stat: Optional[float],
    ) -> CoveredTnJc2:
        """Calculate coverage for a single candidate."""
        # Get coverage in amplicon region (using scaffold-relative positions)
        # pos_scaf_L and pos_scaf_R are 1-based inclusive (from BLAST/junction positions)
        start = raw_tnjc2.pos_scaf_L  # 1-based inclusive
        end = raw_tnjc2.pos_scaf_R  # 1-based inclusive

        iso_region_cov = self.genome.slice_in_range(raw_tnjc2.scaf, iso_scaf_cov, start, end)

        # Calculate raw copy number using scaffold-specific statistic
        iso_mean_cov = mean_positive(iso_region_cov)
        copy_number = iso_mean_cov / iso_scaf_stat if iso_scaf_stat > 0 else 0.0

        amplicon_coverage = None
        copy_number_ratio = None
        amplicon_coverage_mode = None
        anc_region_cov = None

        if self.has_ancestor and anc_scaf_cov is not None and anc_scaf_stat is not None:
            # Get ancestor coverage in same region (using scaffold-relative positions)
            anc_region_cov = self.genome.slice_in_range(raw_tnjc2.scaf, anc_scaf_cov, start, end)

            anc_mean_cov = mean_positive(anc_region_cov)
            anc_copy_number = anc_mean_cov / anc_scaf_stat if anc_scaf_stat > 0 else 0.0

            # Normalized coverage = iso/anc
            amplicon_coverage = copy_number / anc_copy_number if anc_copy_number > 0 else None
            copy_number_ratio = amplicon_coverage

            # Calculate normalized copy number using scaffold-specific statistics
            cp = iso_region_cov.astype(float) / iso_scaf_stat if iso_scaf_stat > 0 else iso_region_cov.astype(float)
            anc_cp = anc_region_cov.astype(float) / anc_scaf_stat if anc_scaf_stat > 0 else anc_region_cov.astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                ncp = cp / anc_cp
            ncp[~np.isfinite(ncp)] = np.nan

            # Calculate distribution mode
            amplicon_coverage_mode = calc_distribution_mode(ncp, x_min=self.ncp_min, x_max=self.ncp_max, n_bins=self.ncp_n, is_log=True)
        else:
            # Raw coverage
            amplicon_coverage = copy_number

            # Calculate raw copy number using scaffold-specific statistic
            ncp = iso_region_cov.astype(float) / iso_scaf_stat if iso_scaf_stat > 0 else iso_region_cov.astype(float)

            # Calculate distribution mode (raw)
            amplicon_coverage_mode = calc_distribution_mode(ncp, x_min=self.ncp_min, x_max=self.ncp_max, n_bins=self.ncp_n, is_log=True)

        # Calculate coverage statistics
        iso_amplicon_cov_stats = calc_coverage_stats(iso_region_cov, average_method=self.average_method)
        scaf_coverage = iso_scaf_stat
        anc_amplicon_cov_stats = (calc_coverage_stats(anc_region_cov, average_method=self.average_method) if anc_region_cov is not None else None)
        anc_scaf_cov_stats = anc_scaf_stat if self.has_ancestor else None

        # Build CoveredTnJc2 from RawTnJc2 + new coverage fields
        return CoveredTnJc2.from_other(
            raw_tnjc2,
            iso_amplicon_coverage=iso_amplicon_cov_stats,
            iso_scaf_coverage=scaf_coverage,
            anc_amplicon_coverage=anc_amplicon_cov_stats,
            anc_scaf_coverage=anc_scaf_cov_stats,
            copy_number=copy_number,
            copy_number_vs_anc=copy_number_ratio,
            amplicon_coverage=amplicon_coverage,
            amplicon_coverage_mode=amplicon_coverage_mode,
        )

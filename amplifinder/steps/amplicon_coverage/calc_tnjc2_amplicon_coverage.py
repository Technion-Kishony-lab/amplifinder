"""Step 7: Calculate amplicon coverage."""

from pathlib import Path
from typing import Optional

import numpy as np

from amplifinder.data_types import RecordTypedDf, RawTnJc2, CoveredTnJc2, Genome, AverageMethod

from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.tools.breseq import load_breseq_coverage

from .coverage import get_scaffold_coverage, calc_average, calc_scaffold_coverages_and_averages
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
            calc_average,
            calc_scaffold_coverages_and_averages,
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
            iso_scaf_covs, iso_scaf_avgs = calc_scaffold_coverages_and_averages(
                iso_cov, unique_scaffolds, self.genome, self.average_method)

        # Load ancestor coverage if provided
        if self.has_ancestor:
            with self.print_timer(f"loading anc coverage from {self.anc_breseq_path} ..."):
                anc_cov = load_breseq_coverage(self.anc_breseq_path, self.genome.name)
            with self.print_timer(f"calculating anc scaffold stats ({len(unique_scaffolds)} scaffolds) ..."):
                anc_scaf_covs, anc_scaf_avgs = calc_scaffold_coverages_and_averages(
                    anc_cov, unique_scaffolds, self.genome, self.average_method)
        else:
            anc_cov = None
            anc_scaf_covs = {}
            anc_scaf_avgs = {}

        # Process each raw_tnjc2s
        with self.print_timer("Calculating coverage for each raw_tnjc2 ..."):
            covered_records = []
            for raw_tnjc2 in self.raw_tnjc2s:
                covered = self._calc_candidate_coverage(
                    raw_tnjc2,
                    iso_scaf_covs[raw_tnjc2.scaf],
                    iso_scaf_avgs[raw_tnjc2.scaf],
                    anc_scaf_covs[raw_tnjc2.scaf] if self.has_ancestor else None,
                    anc_scaf_avgs[raw_tnjc2.scaf] if self.has_ancestor else None,
                )
                covered_records.append(covered)
                self.print(".", end="")

        df = RecordTypedDf.from_records(covered_records, CoveredTnJc2)
        return df

    def calc_average(self, coverage: np.ndarray) -> float:
        """Calculate average coverage."""
        return calc_average(coverage, average_method=self.average_method)

    def _calc_candidate_coverage(
        self,
        raw_tnjc2: RawTnJc2,
        iso_scaf_cov: np.ndarray,
        iso_scaf_avg: float,
        anc_scaf_cov: Optional[np.ndarray],
        anc_scaf_avg: Optional[float],
    ) -> CoveredTnJc2:
        """Calculate coverage for a single candidate."""
        # Get segment scaffold for this amplicon
        seg_scaf = raw_tnjc2.get_segment_scaffold(self.genome)

        iso_amplicon_cov = seg_scaf.slice(seq=iso_scaf_cov)

        # Remove zero coverage regions
        if self.has_ancestor:
            # we mask by the ancestor coverage.
            # if the iso_cov=0 and anc_cov>0 it is meaningful (a deletion)
            anc_amplicon_cov = seg_scaf.slice(seq=anc_scaf_cov)
            mask = anc_amplicon_cov > 0
        else:
            mask = iso_amplicon_cov > 0

        iso_amplicon_cov = iso_amplicon_cov[mask]
        iso_amplicon_avg = self.calc_average(iso_amplicon_cov)

        # Calculate isolate amplicon coverage and copy number
        if self.has_ancestor:
            anc_amplicon_cov = anc_amplicon_cov[mask]
            anc_amplicon_avg = self.calc_average(anc_amplicon_cov)
            norm_cov = (iso_amplicon_cov / anc_amplicon_cov) / (iso_scaf_avg / anc_scaf_avg)
            avg_norm_cov = self.calc_average(norm_cov)
        else:
            avg_norm_cov = iso_amplicon_avg / iso_scaf_avg
            anc_amplicon_avg = None

        # Build CoveredTnJc2 from RawTnJc2 + new coverage fields
        return CoveredTnJc2.from_other(
            raw_tnjc2,
            iso_scaf_avg=iso_scaf_avg,
            iso_amplicon_avg=iso_amplicon_avg,
            anc_scaf_avg=anc_scaf_avg,
            anc_amplicon_avg=anc_amplicon_avg,
            avg_norm_cov=avg_norm_cov,
        )

"""Step 7: Calculate amplicon coverage."""

from pathlib import Path
from typing import Optional

import numpy as np

from amplifinder.data_types import (
    RecordTypedDf, RawTnJc2, CoveredTnJc2, Genome,
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.tools.breseq import load_breseq_coverage
from amplifinder.utils.coverage import (
    get_coverage_in_range,
    calc_genome_coverage,
    calc_distribution_mode,
    get_scaffold_coverage,
    calc_coverage_stats,
)


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
        iso_breseq_path: Path,
        output_dir: Path,
        ref_name: str,
        iso_name: str,
        anc_breseq_path: Optional[Path] = None,
        anc_name: Optional[str] = None,
        ncp_limit1: float = -1,
        ncp_limit2: float = 3,
        ncp_n: int = 150,
        force: Optional[bool] = None,
    ):
        self.raw_tnjc2s = raw_tnjc2s
        self.genome = genome
        self.iso_breseq_path = Path(iso_breseq_path)
        self.ref_name = ref_name
        self.iso_name = iso_name
        self.anc_breseq_path = Path(anc_breseq_path) if anc_breseq_path else None
        self.anc_name = anc_name
        self.ncp_limit1 = ncp_limit1
        self.ncp_limit2 = ncp_limit2
        self.ncp_n = ncp_n
        
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
        """True if ancestor comparison should be performed."""
        return self.anc_breseq_path is not None

    def _calculate_output(self) -> RecordTypedDf[CoveredTnJc2]:
        """Calculate coverage for each TnJc2 candidate."""
        # Load isolate coverage
        iso_cov = load_breseq_coverage(self.iso_breseq_path, self.genome.name)
        iso_genome_median = calc_genome_coverage(iso_cov)
        
        # Load ancestor coverage if provided
        anc_cov = None
        anc_genome_median = None
        if self.has_ancestor:
            anc_cov = load_breseq_coverage(self.anc_breseq_path, self.genome.name)
            anc_genome_median = calc_genome_coverage(anc_cov)
        
        # Process each candidate
        covered_records = []
        for raw_tnjc2 in self.raw_tnjc2s:
            covered = self._calc_candidate_coverage(
                raw_tnjc2, iso_cov, iso_genome_median, anc_cov, anc_genome_median
            )
            covered_records.append(covered)
        
        return RecordTypedDf.from_records(covered_records, CoveredTnJc2)

    def _calc_candidate_coverage(
        self,
        raw_tnjc2: RawTnJc2,
        iso_cov: np.ndarray,
        iso_genome_median: float,
        anc_cov: Optional[np.ndarray],
        anc_genome_median: Optional[float],
    ) -> CoveredTnJc2:
        """Calculate coverage for a single candidate."""
        # Get coverage in amplicon region
        start = raw_tnjc2.pos_scaf_L
        end = raw_tnjc2.pos_scaf_R
        span_origin = raw_tnjc2.span_origin
        
        iso_region_cov = get_coverage_in_range(
            iso_cov, start, end, span_origin, self.genome.length
        )
        
        # Calculate raw copy number
        iso_mean_cov = float(np.mean(iso_region_cov[iso_region_cov > 0])) if np.any(iso_region_cov > 0) else 0.0
        copy_number = iso_mean_cov / iso_genome_median if iso_genome_median > 0 else 0.0
        
        amplicon_coverage = None
        copy_number_ratio = None
        amplicon_coverage_mode = None
        
        if self.has_ancestor:
            # Get ancestor coverage in same region
            anc_region_cov = get_coverage_in_range(
                anc_cov, start, end, span_origin, self.genome.length
            )
            anc_mean_cov = float(np.mean(anc_region_cov[anc_region_cov > 0])) if np.any(anc_region_cov > 0) else 0.0
            anc_copy_number = anc_mean_cov / anc_genome_median if anc_genome_median > 0 else 0.0
            
            # Normalized coverage = iso/anc
            amplicon_coverage = copy_number / anc_copy_number if anc_copy_number > 0 else None
            copy_number_ratio = amplicon_coverage
            
            # Calculate normalized copy number
            cp = iso_region_cov.astype(float) / iso_genome_median if iso_genome_median > 0 else iso_region_cov.astype(float)
            anc_cp = anc_region_cov.astype(float) / anc_genome_median if anc_genome_median > 0 else anc_region_cov.astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                ncp = cp / anc_cp
            ncp[~np.isfinite(ncp)] = np.nan
            
            # Calculate distribution mode
            amplicon_coverage_mode = calc_distribution_mode(
                ncp,
                x_min=10**self.ncp_limit1,
                x_max=10**self.ncp_limit2,
                n_bins=self.ncp_n,
                is_log=True,
            )
        else:
            # Raw coverage
            amplicon_coverage = copy_number
            
            # Calculate raw copy number
            ncp = iso_region_cov.astype(float) / iso_genome_median if iso_genome_median > 0 else iso_region_cov.astype(float)
            
            # Calculate distribution mode (raw)
            amplicon_coverage_mode = calc_distribution_mode(
                ncp,
                x_min=10**self.ncp_limit1,
                x_max=10**self.ncp_limit2,
                n_bins=self.ncp_n,
                is_log=True,
            )
        
        # Calculate scaffold-specific coverage statistics
        scaf_cov = get_scaffold_coverage(iso_cov, raw_tnjc2.scaf, self.genome)
        scaf_coverage = calc_coverage_stats(scaf_cov)
        
        # Build CoveredTnJc2 from RawTnJc2 + new coverage fields
        return CoveredTnJc2.from_other(
            raw_tnjc2,
            amplicon_coverage=amplicon_coverage,
            scaf_coverage=scaf_coverage,
            copy_number=copy_number,
            amplicon_coverage_mode=amplicon_coverage_mode,
            copy_number_ratio=copy_number_ratio,
        )

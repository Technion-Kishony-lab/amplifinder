"""Step 7: Calculate amplicon coverage."""

from pathlib import Path
from typing import Optional

import numpy as np

from amplifinder.data_types import (
    RecordTypedDf, TnJc2, CoveredTnJc2, Genome,
)
from amplifinder.steps.base import RecordTypedDfStep
from amplifinder.tools.breseq import load_breseq_coverage
from amplifinder.utils.coverage import (
    get_coverage_in_range,
    calc_genome_coverage,
    calc_copy_number_distribution,
)


class CalcAmpliconCoverageStep(RecordTypedDfStep[CoveredTnJc2]):
    """Calculate amplicon coverage for TnJc2 candidates.
    
    Coverage calculation depends on run type:
    - anc_breseq_path=None: raw coverage only
    - anc_breseq_path=set: normalized coverage (iso/anc ratio)
    """

    def __init__(
        self,
        tnjc2: RecordTypedDf[TnJc2],
        genome: Genome,
        iso_breseq_path: Path,
        output_dir: Path,
        ref_name: str,
        iso_name: str,
        anc_breseq_path: Optional[Path] = None,
        anc_name: Optional[str] = None,
        min_amplicon_length: int = 30,
        max_amplicon_length: int = 1_000_000,
        ncp_limit1: float = -1,
        ncp_limit2: float = 3,
        ncp_n: int = 150,
        force: Optional[bool] = None,
    ):
        self.tnjc2 = tnjc2
        self.genome = genome
        self.iso_breseq_path = Path(iso_breseq_path)
        self.ref_name = ref_name
        self.iso_name = iso_name
        self.anc_breseq_path = Path(anc_breseq_path) if anc_breseq_path else None
        self.anc_name = anc_name
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
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
        for tnjc in self.tnjc2:
            covered = self._calc_candidate_coverage(
                tnjc, iso_cov, iso_genome_median, anc_cov, anc_genome_median
            )
            covered_records.append(covered)
        
        return RecordTypedDf.from_records(covered_records, CoveredTnJc2)

    def _calc_candidate_coverage(
        self,
        tnjc: TnJc2,
        iso_cov: np.ndarray,
        iso_genome_median: float,
        anc_cov: Optional[np.ndarray],
        anc_genome_median: Optional[float],
    ) -> CoveredTnJc2:
        """Calculate coverage for a single candidate."""
        amplicon_length = tnjc.amplicon_length
        
        # Default values for candidates outside length range
        amplicon_coverage = None
        copy_number = None
        copy_number_ratio = None
        amplicon_coverage_mode = None
        
        # Only calculate for valid amplicon lengths
        if self.min_amplicon_length < amplicon_length < self.max_amplicon_length:
            # Get coverage in amplicon region
            start = tnjc.pos_chr_L
            end = tnjc.pos_chr_R
            span_origin = tnjc.span_origin
            
            iso_region_cov = get_coverage_in_range(
                iso_cov, start, end, span_origin, self.genome.length
            )
            
            # Calculate raw copy number
            iso_mean_cov = float(np.mean(iso_region_cov[iso_region_cov > 0])) if np.any(iso_region_cov > 0) else 0.0
            copy_number = iso_mean_cov / iso_genome_median if iso_genome_median > 0 else 0.0
            
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
                
                # Calculate distribution mode
                _, amplicon_coverage_mode = calc_copy_number_distribution(
                    iso_region_cov, iso_genome_median,
                    anc_region_cov, anc_genome_median,
                    self.ncp_limit1, self.ncp_limit2, self.ncp_n,
                )
            else:
                # Raw coverage
                amplicon_coverage = copy_number
                
                # Calculate distribution mode (raw)
                _, amplicon_coverage_mode = calc_copy_number_distribution(
                    iso_region_cov, iso_genome_median,
                    ncp_limit1=self.ncp_limit1, ncp_limit2=self.ncp_limit2, ncp_n=self.ncp_n,
                )
        
        # Build CoveredTnJc2 from TnJc2 + new coverage fields
        return CoveredTnJc2.from_other(
            tnjc,
            ref_name=self.ref_name,
            iso_name=self.iso_name,
            anc_name=self.anc_name,
            amplicon_coverage=amplicon_coverage,
            genome_coverage=iso_genome_median,
            copy_number=copy_number,
            amplicon_coverage_mode=amplicon_coverage_mode,
            copy_number_ratio=copy_number_ratio,
        )

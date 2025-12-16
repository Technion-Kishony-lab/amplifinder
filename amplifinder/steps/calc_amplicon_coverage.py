"""Step 7: Calculate amplicon coverage."""

from pathlib import Path
from typing import Optional

import numpy as np

from amplifinder.data_types import (
    RecordTypedDF, TnJc2, CoveredTnJc2, Genome,
)
from amplifinder.steps.base import Step
from amplifinder.utils.coverage import (
    load_breseq_coverage,
    get_coverage_in_range,
    calc_genome_coverage,
    calc_copy_number_distribution,
)


class CalcAmpliconCoverageStep(Step[RecordTypedDF[CoveredTnJc2]]):
    """Calculate amplicon coverage for TnJc2 candidates.
    
    Coverage calculation depends on run type:
    - anc_breseq_path=None: raw coverage only
    - anc_breseq_path=set: normalized coverage (iso/anc ratio)
    """

    def __init__(
        self,
        tnjc2: RecordTypedDF[TnJc2],
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
        self.output_dir = Path(output_dir)
        self.ref_name = ref_name
        self.iso_name = iso_name
        self.anc_breseq_path = Path(anc_breseq_path) if anc_breseq_path else None
        self.anc_name = anc_name
        self.min_amplicon_length = min_amplicon_length
        self.max_amplicon_length = max_amplicon_length
        self.ncp_limit1 = ncp_limit1
        self.ncp_limit2 = ncp_limit2
        self.ncp_n = ncp_n
        
        self.output_file = output_dir / "tn_jc2_covered.csv"
        
        input_files = [iso_breseq_path]
        if anc_breseq_path:
            input_files.append(anc_breseq_path)
        
        super().__init__(
            input_files=input_files,
            output_files=[self.output_file],
            force=force,
        )

    @property
    def has_ancestor(self) -> bool:
        """True if ancestor comparison should be performed."""
        return self.anc_breseq_path is not None

    def _calculate_output(self) -> RecordTypedDF[CoveredTnJc2]:
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
        for _, row in self.tnjc2.df.iterrows():
            covered = self._calc_candidate_coverage(
                row, iso_cov, iso_genome_median, anc_cov, anc_genome_median
            )
            covered_records.append(covered)
        
        return RecordTypedDF.from_records(covered_records, CoveredTnJc2)

    def _calc_candidate_coverage(
        self,
        row,
        iso_cov: np.ndarray,
        iso_genome_median: float,
        anc_cov: Optional[np.ndarray],
        anc_genome_median: Optional[float],
    ) -> CoveredTnJc2:
        """Calculate coverage for a single candidate."""
        amplicon_length = row["amplicon_length"]
        
        # Default values for candidates outside length range
        amplicon_coverage = float("nan")
        copy_number = float("nan")
        copy_number_ratio = None
        amplicon_coverage_mode = float("nan")
        
        # Only calculate for valid amplicon lengths
        if self.min_amplicon_length < amplicon_length < self.max_amplicon_length:
            # Get coverage in amplicon region
            start = row["pos_chr_L"]
            end = row["pos_chr_R"]
            span_origin = row["span_origin"]
            
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
                amplicon_coverage = copy_number / anc_copy_number if anc_copy_number > 0 else float("nan")
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
        
        # Build CoveredTnJc2 from existing TnJc2 fields + new coverage fields
        return CoveredTnJc2(
            # From TnJc2
            jc_num_L=row["jc_num_L"],
            jc_num_R=row["jc_num_R"],
            scaf_chr=row["scaf_chr"],
            pos_chr_L=row["pos_chr_L"],
            pos_chr_R=row["pos_chr_R"],
            pos_tn_L=row["pos_tn_L"],
            pos_tn_R=row["pos_tn_R"],
            dir_chr_L=row["dir_chr_L"],
            dir_chr_R=row["dir_chr_R"],
            dir_tn_L=row["dir_tn_L"],
            dir_tn_R=row["dir_tn_R"],
            tn_ids=row["tn_ids"],
            tn_orientations=row["tn_orientations"],
            span_origin=row["span_origin"],
            amplicon_length=row["amplicon_length"],
            complementary_length=row["complementary_length"],
            # New coverage fields
            ref_name=self.ref_name,
            iso_name=self.iso_name,
            anc_name=self.anc_name,
            amplicon_coverage=amplicon_coverage,
            genome_coverage=iso_genome_median,
            copy_number=copy_number,
            amplicon_coverage_mode=amplicon_coverage_mode,
            copy_number_ratio=copy_number_ratio,
        )

    def _save_output(self, output: RecordTypedDF[CoveredTnJc2]) -> None:
        """Save covered TnJc2 to CSV."""
        output.to_csv(self.output_file)

    def load_outputs(self) -> RecordTypedDF[CoveredTnJc2]:
        """Load covered TnJc2 from CSV."""
        return RecordTypedDF.from_csv(self.output_file, CoveredTnJc2)


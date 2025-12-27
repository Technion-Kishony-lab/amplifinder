"""Coverage analysis utilities."""
from __future__ import annotations

from typing import Optional, List, Dict, TYPE_CHECKING

import numpy as np

from amplifinder.data_types import Coverage
from amplifinder.steps.amplicon_coverage.statistics import calc_distribution_mode

if TYPE_CHECKING:
    from amplifinder.data_types.genome import Genome


def get_scaffold_coverage(
    cov: np.ndarray,
    scaffold_name: str,
    genome: Genome,
    ) -> np.ndarray:
    """Get coverage array for a specific scaffold.
    
    Args:
        cov: Full genome coverage array (concatenated by scaffold order, 0-based indexing)
        scaffold_name: Name of scaffold to extract
        genome: Genome object with scaffold information
    
    Returns:
        Coverage array for the specified scaffold (0-based indexing)
    
    Note:
        Uses 0-based positions from scaffold_ranges property for array slicing.
    """
    start, end = genome.scaffold_ranges[scaffold_name]
    return cov[start:end]


def get_coverage_in_range(
    cov: np.ndarray, 
    start: int, 
    end: int, 
    span_origin: bool = False,
    genome_length: Optional[int] = None,
) -> np.ndarray:
    """Get coverage values in range [start, end] (1-based, inclusive).
    
    Args:
        cov: Full coverage array (0-based indexing)
        start: Start position (1-based, inclusive) - converted to 0-based internally
        end: End position (1-based, inclusive) - converted to 0-based exclusive internally
        span_origin: If True, region spans circular genome origin
        genome_length: Required if span_origin=True
    
    Returns:
        Coverage values in the specified range
    
    Note:
        This function accepts 1-based genomic coordinates (matching BLAST/GenBank)
        and converts them to 0-based array indices internally.
    """
    if not span_origin:
        return cov[start - 1:end]
    if genome_length is None:
        genome_length = len(cov)
    # Region spans origin: concatenate end-to-origin and origin-to-start
    return np.concatenate([cov[start - 1:], cov[:end]])


def calc_coverage_stats(cov: np.ndarray, include_zeros: bool = False) -> Coverage:
    """Calculate mean, median and mode of coverage array.
    
    Args:
        cov: Coverage values
        include_zeros: If False (default), exclude zero values from statistics
    
    Returns:
        Coverage namedtuple with mean, median, mode
    """
    if not include_zeros:
        cov = cov[cov > 0]
    
    if len(cov) == 0:
        return Coverage(mean=0.0, median=0.0, mode=0.0)
    
    return Coverage(
        mean=float(np.mean(cov)), 
        median=float(np.median(cov)), 
        mode=calc_distribution_mode(cov, is_log=False, skip_first_bin=False)
    )


def calc_median_coverage(cov: np.ndarray) -> float:
    """Calculate median coverage for normalization (excluding zeros).
    
    Args:
        cov: Coverage array (can be genome-wide, scaffold-specific, or any region)
    
    Returns:
        Median coverage (excluding zeros)
    """
    return calc_coverage_stats(cov, include_zeros=False).median


def calc_scaffold_coverage_median(
    cov: np.ndarray,
    scaffold_name: str,
    genome: Genome,
) -> float:
    """Calculate median coverage for a specific scaffold.
    
    Args:
        cov: Full genome coverage array (concatenated by scaffold order)
        scaffold_name: Name of scaffold to extract
        genome: Genome object with scaffold information
    
    Returns:
        Median coverage for the scaffold (excluding zeros)
    """
    scaf_cov = get_scaffold_coverage(cov, scaffold_name, genome)
    return calc_median_coverage(scaf_cov)


def calc_scaffold_coverages_and_medians(
    cov: np.ndarray,
    scaffold_names: List[str],
    genome: Genome,
) -> tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """Calculate coverage arrays and medians for multiple scaffolds.
    
    Args:
        cov: Full genome coverage array (concatenated by scaffold order)
        scaffold_names: List of scaffold names to process
        genome: Genome object with scaffold information
    
    Returns:
        Tuple of (scaffold_coverage_arrays, scaffold_medians) dictionaries
    """
    scaf_covs = {}
    scaf_medians = {}
    for scaf in scaffold_names:
        scaf_cov = get_scaffold_coverage(cov, scaf, genome)
        scaf_covs[scaf] = scaf_cov
        scaf_medians[scaf] = calc_median_coverage(scaf_cov)
    return scaf_covs, scaf_medians

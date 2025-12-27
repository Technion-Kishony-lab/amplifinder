"""Coverage analysis utilities."""
from __future__ import annotations

from typing import Optional, List, Dict, TYPE_CHECKING

import numpy as np

from amplifinder.data_types import Coverage

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


def calc_distribution_mode(
    x: np.ndarray,
    x_min: Optional[float] = None,
    x_max: Optional[float] = None,
    n_bins: Optional[int] = None,
    is_log: bool = True,
    skip_first_bin: bool = True,
) -> float:
    """Calculate mode of distribution.
    
    Args:
        x: Data array
        x_min: Minimum threshold value. If None, uses min of data
        x_max: Maximum threshold value. If None, uses max of data
        n_bins: Number of histogram bins. If None, uses Freedman-Diaconis rule
        is_log: If True, work in log space (take log of x)
        skip_first_bin: If True, skip first bin when finding mode (default: True)
    
    Returns:
        Mode value
    """
    # Remove NaN and infinite values
    x_clean = x[np.isfinite(x)]
    
    if len(x_clean) == 0:
        return np.nan
    
    # Set defaults from data if not provided
    if x_min is None:
        x_min = np.min(x_clean)
    if x_max is None:
        x_max = np.max(x_clean)
    
    # Transform to log space if needed
    if is_log:
        x_clean = np.log10(x_clean)
        x_min = np.log10(x_min)
        x_max = np.log10(x_max)
    
    # Clamp to limits
    x_clamped = np.clip(x_clean, x_min, x_max)
    
    # Auto-determine bins using Freedman-Diaconis rule if not specified
    if n_bins is None:
        n = len(x_clamped)
        q75, q25 = np.percentile(x_clamped, [75, 25])
        iqr = q75 - q25
        if iqr > 0:
            bin_width = 2 * iqr / (n ** (1/3))
            n_bins = max(10, int((x_max - x_min) / bin_width))
        else:
            # Fallback for constant data: use square root rule
            n_bins = max(10, min(50, int(np.sqrt(n))))
    
    # Build histogram (always linear bins after transformation)
    edges = np.linspace(x_min, x_max, n_bins)
    hist, _ = np.histogram(x_clamped, bins=edges)
    
    # Find mode
    if skip_first_bin and len(hist) > 1:
        max_idx = np.argmax(hist[1:]) + 1  # offset by 1 since we skipped first
    else:
        max_idx = np.argmax(hist)
    mode = (edges[max_idx] + edges[max_idx + 1]) / 2  # arithmetic mean
    
    # Transform back from log space if needed
    if is_log:
        mode = 10**mode
    
    return float(mode)


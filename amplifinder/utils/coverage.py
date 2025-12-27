"""Coverage analysis utilities."""
from __future__ import annotations

from typing import Optional, List, TYPE_CHECKING

import numpy as np

from amplifinder.data_types import Coverage

if TYPE_CHECKING:
    from amplifinder.data_types.genome import Genome


def get_coverage_in_range(
    cov: np.ndarray, 
    start: int, 
    end: int, 
    span_origin: bool = False,
    genome_length: Optional[int] = None,
) -> np.ndarray:
    """Get coverage values in range [start, end] (1-based, inclusive).
    
    Args:
        cov: Full coverage array
        start: Start position (1-based)
        end: End position (1-based)
        span_origin: If True, region spans circular genome origin
        genome_length: Required if span_origin=True
    
    Returns:
        Coverage values in the specified range
    """
    if span_origin:
        if genome_length is None:
            genome_length = len(cov)
        # Region spans origin: concatenate end-to-origin and origin-to-start
        return np.concatenate([cov[start - 1:], cov[:end]])
    else:
        return cov[start - 1:end]


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


def calc_genome_coverage(cov: np.ndarray) -> float:
    """Calculate genome-wide median coverage for normalization.
    
    Args:
        cov: Full genome coverage array
    
    Returns:
        Median coverage (excluding zeros)
    """
    return calc_coverage_stats(cov, include_zeros=False).median


def get_scaffold_coverage(
    cov: np.ndarray,
    scaffold_name: str,
    genome: Genome,
) -> np.ndarray:
    """Get coverage array for a specific scaffold.
    
    Args:
        cov: Full genome coverage array (concatenated by scaffold order)
        scaffold_name: Name of scaffold to extract
        genome: Genome object with scaffold information
    
    Returns:
        Coverage array for the specified scaffold
    """
    start, end = genome.get_scaffold_coverage_range(scaffold_name)
    return cov[start:end]


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


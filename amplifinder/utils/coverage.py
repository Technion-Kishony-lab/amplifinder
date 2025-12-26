"""Coverage analysis utilities."""

from typing import Optional

import numpy as np

from amplifinder.data_types import Coverage


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
    
    mean = float(np.mean(cov))
    median = float(np.median(cov))
    
    # Mode calculation using histogram (numpy-only approach)
    # Works for len >= 1: histogram creates appropriate bins and mode is bin midpoint
    hist, bin_edges = np.histogram(cov, bins='auto')
    max_bin = np.argmax(hist)
    mode = float((bin_edges[max_bin] + bin_edges[max_bin + 1]) / 2)
    
    return Coverage(mean=mean, median=median, mode=mode)


def calc_genome_coverage(cov: np.ndarray) -> float:
    """Calculate genome-wide median coverage for normalization.
    
    Args:
        cov: Full genome coverage array
    
    Returns:
        Median coverage (excluding zeros)
    """
    cov_nonzero = cov[cov > 0]
    return float(np.median(cov_nonzero)) if len(cov_nonzero) > 0 else 0.0


def calc_copy_number_distribution(
    cov: np.ndarray,
    genome_median: float,
    anc_cov: Optional[np.ndarray] = None,
    anc_genome_median: Optional[float] = None,
    ncp_limit1: float = -1,
    ncp_limit2: float = 3,
    ncp_n: int = 150,
) -> tuple[np.ndarray, float]:
    """Calculate copy number distribution and mode.
    
    Based on MATLAB calc_coverage_ISJC2.m
    
    Args:
        cov: Isolate coverage in region
        genome_median: Isolate genome median coverage
        anc_cov: Ancestor coverage in region (optional)
        anc_genome_median: Ancestor genome median coverage (optional)
        ncp_limit1: Log10 lower limit for distribution
        ncp_limit2: Log10 upper limit for distribution
        ncp_n: Number of histogram bins
    
    Returns:
        (histogram_counts, mode_value)
    """
    # Calculate copy number
    cp = cov.astype(float) / genome_median
    
    if anc_cov is not None and anc_genome_median is not None:
        # Normalized copy number
        anc_cp = anc_cov.astype(float) / anc_genome_median
        with np.errstate(divide='ignore', invalid='ignore'):
            ncp = cp / anc_cp
        ncp[~np.isfinite(ncp)] = np.nan
    else:
        # Raw copy number
        ncp = cp
    
    # Clamp to limits
    ncp = np.clip(ncp, 10**ncp_limit1, 10**ncp_limit2)
    
    # Build histogram
    edges = np.logspace(ncp_limit1, ncp_limit2, ncp_n)
    hist, _ = np.histogram(ncp[~np.isnan(ncp)], bins=edges)
    
    # Find mode (skip first bin to avoid deletions)
    if len(hist) > 1:
        max_idx = np.argmax(hist[1:]) + 1  # offset by 1 since we skipped first
        mode = np.sqrt(edges[max_idx] * edges[max_idx + 1])  # geometric mean of bin edges
    else:
        mode = 1.0
    
    return hist, float(mode)


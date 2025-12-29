"""Coverage analysis utilities."""
from __future__ import annotations

from typing import Optional, List, Dict, TYPE_CHECKING, Union

import numpy as np

from amplifinder.data_types import AverageMethod
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


def calc_coverage_stats(cov: np.ndarray, average_method: Union[AverageMethod, str] = AverageMethod.MEDIAN, include_zeros: bool = False) -> float:
    """Calculate coverage statistic based on average_method.

    Args:
        cov: Coverage values
        average_method: Method to use (AverageMethod enum or string)
        include_zeros: If False (default), exclude zero values from statistics

    Returns:
        Coverage statistic value (float)
    """
    if not include_zeros:
        cov = cov[cov > 0]

    if len(cov) == 0:
        return 0.0

    # Convert string to enum if needed
    if isinstance(average_method, str):
        average_method = AverageMethod(average_method)

    if average_method == AverageMethod.MEAN:
        return float(np.mean(cov))
    elif average_method == AverageMethod.MEDIAN:
        return float(np.median(cov))
    elif average_method == AverageMethod.MODE:
        return calc_distribution_mode(cov, is_log=False, skip_first_bin=False)
    else:
        raise ValueError(f"Invalid average_method: {average_method}")


def mean_positive(cov: np.ndarray) -> float:
    """Calculate mean of positive values in coverage array.

    Efficiently computes mean of values > 0.

    Args:
        cov: Coverage values

    Returns:
        Mean of positive values, or 0.0 if no positive values exist
    """
    mask = cov > 0
    n_positive = np.sum(mask)
    if n_positive == 0:
        return 0.0
    return float(cov.sum() / n_positive)


def calc_scaffold_coverages_and_stats(
    cov: np.ndarray,
    scaffold_names: List[str],
    genome: Genome,
    average_method: Union[AverageMethod, str] = AverageMethod.MEDIAN,
) -> tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """Calculate coverage arrays and statistics for multiple scaffolds.

    Args:
        cov: Full genome coverage array (concatenated by scaffold order)
        scaffold_names: List of scaffold names to process
        genome: Genome object with scaffold information
        average_method: Method to use ('median', 'mode', or 'mean')

    Returns:
        Tuple of (scaffold_coverage_arrays, scaffold_stats) dictionaries
    """
    scaf_covs = {}
    scaf_stats = {}
    for scaf in scaffold_names:
        scaf_cov = get_scaffold_coverage(cov, scaf, genome)
        scaf_covs[scaf] = scaf_cov
        scaf_stats[scaf] = calc_coverage_stats(scaf_cov, average_method=average_method)
    return scaf_covs, scaf_stats

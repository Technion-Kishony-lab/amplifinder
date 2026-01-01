"""Coverage analysis utilities."""
from __future__ import annotations

from typing import List, Dict, TYPE_CHECKING, Union

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


def calc_average(cov: np.ndarray, average_method: AverageMethod = AverageMethod.MEDIAN) -> float:
    """Calculate coverage statistic based on average_method.

    Args:
        cov: Coverage values
        average_method: Method to use (AverageMethod enum or string)
        exclude_zeros: If False (default), include zero values in statistics

    Returns:
        Coverage statistic value (float)
    """
    if len(cov) == 0:
        return np.nan

    if average_method == AverageMethod.MEAN:
        return float(np.mean(cov))
    elif average_method == AverageMethod.MEDIAN:
        return float(np.median(cov))
    elif average_method == AverageMethod.MODE:
        return calc_distribution_mode(cov, is_log=True)
    else:
        raise ValueError(f"Invalid average_method: {average_method}")


def calc_scaffold_coverages_and_averages(
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
        scaf_stats[scaf] = calc_average(scaf_cov, average_method=average_method)
    return scaf_covs, scaf_stats

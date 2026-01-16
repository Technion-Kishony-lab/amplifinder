"""Coverage analysis utilities."""
from __future__ import annotations

from typing import List, Dict, Union

import numpy as np

from amplifinder.data_types import AverageMethod
from amplifinder.steps.amplicon_coverage.statistics import calc_distribution_mode


def calc_average(cov: np.ndarray, average_method: AverageMethod = AverageMethod.MEDIAN) -> float:
    """Calculate coverage statistic based on average_method.

    Args:
        cov: Coverage values array
        average_method: Method to use (AverageMethod enum or string)

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
    scafs_to_covs: Dict[str, np.ndarray],
    scaffold_names: List[str],
    average_method: Union[AverageMethod, str] = AverageMethod.MEDIAN,
) -> tuple[Dict[str, np.ndarray], Dict[str, float]]:
    """Calculate coverage arrays and statistics for multiple scaffolds.

    Args:
        scafs_to_covs: Dictionary mapping scaffold names to coverage arrays
        scaffold_names: List of scaffold names to process
        average_method: Method to use ('median', 'mode', or 'mean')

    Returns:
        Tuple of (scaffold_coverage_arrays, scaffold_stats) dictionaries
    """
    scaf_covs = {}
    scaf_stats = {}
    for scaf in scaffold_names:
        scaf_cov = scafs_to_covs[scaf]
        scaf_covs[scaf] = scaf_cov
        scaf_stats[scaf] = calc_average(scaf_cov, average_method=average_method)
    return scaf_covs, scaf_stats

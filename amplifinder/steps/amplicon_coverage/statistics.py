"""Statistical utilities for coverage analysis."""

from typing import Optional

import numpy as np


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
            bin_width = 2 * iqr / (n ** (1 / 3))
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

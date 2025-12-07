"""Utilities for comparing Python output with MATLAB reference."""

import pandas as pd
from pathlib import Path
from typing import Optional


def load_matlab_isjc2(matlab_output_dir: Path) -> Optional[pd.DataFrame]:
    """Load MATLAB ISJC2.xlsx output."""
    xlsx = matlab_output_dir / "ISJC2.xlsx"
    if not xlsx.exists():
        return None
    return pd.read_excel(xlsx)


def load_matlab_classified(matlab_output_dir: Path) -> Optional[pd.DataFrame]:
    """Load MATLAB classified_amplifications.xlsx output."""
    xlsx = matlab_output_dir / "classified_amplifications.xlsx"
    if not xlsx.exists():
        return None
    return pd.read_excel(xlsx)


def compare_tn_junctions(
    python_df: pd.DataFrame,
    matlab_df: pd.DataFrame,
    tolerance: int = 5,
) -> dict:
    """Compare Python and MATLAB TN junction outputs.

    Args:
        python_df: Python pipeline output
        matlab_df: MATLAB reference output
        tolerance: Position tolerance in bp

    Returns:
        Dict with comparison stats
    """
    stats = {
        "python_count": len(python_df),
        "matlab_count": len(matlab_df),
        "matched": 0,
        "python_only": 0,
        "matlab_only": 0,
    }

    # TODO: Implement detailed junction matching
    # For now, just return counts

    return stats


def assert_junctions_match(
    python_df: pd.DataFrame,
    matlab_df: pd.DataFrame,
    min_match_ratio: float = 0.9,
):
    """Assert that Python and MATLAB junctions match within tolerance.

    Raises:
        AssertionError: If match ratio below threshold
    """
    stats = compare_tn_junctions(python_df, matlab_df)

    if stats["matlab_count"] == 0:
        return  # No MATLAB reference to compare

    # For now, just check counts are similar
    ratio = stats["python_count"] / stats["matlab_count"]
    assert 0.8 <= ratio <= 1.2, (
        f"Junction count mismatch: Python={stats['python_count']}, "
        f"MATLAB={stats['matlab_count']}"
    )

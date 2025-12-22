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


def match_junctions(
    python_df: pd.DataFrame,
    matlab_df: pd.DataFrame,
    pos_tolerance: int = 5,
):
    """Match Python and MATLAB junctions by position (1-to-1 matching).
    
    Args:
        python_df: Python pipeline output
        matlab_df: MATLAB reference output
        pos_tolerance: Position tolerance in bp
    
    Returns:
        matches: List of (python_idx, matlab_idx) tuples
        python_matched: Set of matched Python indices
        matlab_matched: Set of matched MATLAB indices
    """
    matches = []
    python_matched = set()
    matlab_matched = set()
    
    # Build candidate matches (may have multiple candidates per junction)
    candidates = []
    for i, matlab_row in matlab_df.iterrows():
        matlab_pos1 = matlab_row['Positions_in_chromosome_1']
        matlab_pos2 = matlab_row['Positions_in_chromosome_2']
        
        for j, python_row in python_df.iterrows():
            # Parse Python positions (format: "15386-16732")
            py_pos_str = python_row['Positions_in_chromosome']
            if isinstance(py_pos_str, str) and '-' in py_pos_str:
                py_positions = py_pos_str.split('-')
                py_pos1, py_pos2 = int(py_positions[0]), int(py_positions[1])
            else:
                continue  # Skip if format is wrong
            
            # Check if positions match within tolerance
            pos1_diff = abs(py_pos1 - matlab_pos1)
            pos2_diff = abs(py_pos2 - matlab_pos2)
            if pos1_diff <= pos_tolerance and pos2_diff <= pos_tolerance:
                total_diff = pos1_diff + pos2_diff
                candidates.append((j, i, total_diff))
    
    # Sort by distance (best matches first)
    candidates.sort(key=lambda x: x[2])
    
    # Greedy matching: assign best matches first, ensuring 1-to-1
    for py_idx, matlab_idx, distance in candidates:
        if py_idx not in python_matched and matlab_idx not in matlab_matched:
            matches.append((py_idx, matlab_idx))
            python_matched.add(py_idx)
            matlab_matched.add(matlab_idx)
    
    return matches, python_matched, matlab_matched


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
    matches, py_matched, matlab_matched = match_junctions(
        python_df, matlab_df, tolerance
    )
    
    stats = {
        "python_count": len(python_df),
        "matlab_count": len(matlab_df),
        "matched": len(matches),
        "python_only": len(python_df) - len(py_matched),
        "matlab_only": len(matlab_df) - len(matlab_matched),
    }
    
    return stats


def compare_isjc2_outputs(
    python_df: pd.DataFrame,
    matlab_df: pd.DataFrame,
    pos_tolerance: int = 5,
    length_tolerance: int = 10,
):
    """Compare ISJC2 outputs with 1-to-1 matching requirement.
    
    Args:
        python_df: Python ISJC2 output
        matlab_df: MATLAB ISJC2 output
        pos_tolerance: Position tolerance in bp
        length_tolerance: Amplicon length tolerance in bp
    
    Raises:
        AssertionError: If 1-to-1 matching fails or properties don't match
    """
    matches, py_matched, matlab_matched = match_junctions(
        python_df, matlab_df, pos_tolerance
    )
    
    # Require 1-to-1 match: every junction must have exactly one match
    assert len(matches) == len(matlab_df), (
        f"Not all MATLAB junctions matched: {len(matches)}/{len(matlab_df)} matched. "
        f"Missing MATLAB junctions: {len(matlab_df) - len(matlab_matched)}"
    )
    
    assert len(matches) == len(python_df), (
        f"Not all Python junctions matched: {len(matches)}/{len(python_df)} matched. "
        f"Extra Python junctions: {len(python_df) - len(py_matched)}"
    )
    
    # Verify no duplicates
    assert len(py_matched) == len(matches), "Duplicate Python junction matches"
    assert len(matlab_matched) == len(matches), "Duplicate MATLAB junction matches"
    
    # Compare matched junctions
    for py_idx, matlab_idx in matches:
        py_row = python_df.iloc[py_idx]
        matlab_row = matlab_df.iloc[matlab_idx]
        
        # Compare amplicon length
        py_len = py_row['amplicon_length']
        matlab_len = matlab_row['amplicon_length']
        assert abs(py_len - matlab_len) <= length_tolerance, (
            f"Length mismatch at positions {py_row.get('Positions_in_chromosome', 'unknown')}: "
            f"Python={py_len}, MATLAB={matlab_len}"
        )
        
        # Compare IS elements (order-independent)
        py_is_str = str(py_row.get('IS_element', ''))
        matlab_is_str = str(matlab_row.get('IS_element', ''))
        py_is = set(x.strip() for x in py_is_str.split(',') if x.strip())
        matlab_is = set(x.strip() for x in matlab_is_str.split(',') if x.strip())
        assert py_is == matlab_is, (
            f"IS elements mismatch at positions {py_row.get('Positions_in_chromosome', 'unknown')}: "
            f"Python={py_is}, MATLAB={matlab_is}"
        )

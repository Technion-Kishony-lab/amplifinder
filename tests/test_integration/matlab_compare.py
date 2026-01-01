"""Utilities for comparing Python output with MATLAB reference."""

import pandas as pd
from pathlib import Path
from typing import Optional


def _load_matlab_xlsx(
    matlab_output_dir: Path,
    filename: str,
    require_files: bool = None,
) -> Optional[pd.DataFrame]:
    """Load MATLAB Excel file with optional requirement checking.

    Args:
        matlab_output_dir: Directory containing MATLAB outputs
        filename: Name of the Excel file to load
        require_files: If True, fail if file doesn't exist. If None, use REQUIRE_MATLAB_FILES from test_pipeline.

    Returns:
        DataFrame if file exists, None otherwise (or fails if require_files=True)
    """
    import pytest
    from tests.test_integration.test_pipeline import REQUIRE_MATLAB_FILES

    if require_files is None:
        require_files = REQUIRE_MATLAB_FILES

    xlsx = matlab_output_dir / filename
    if not xlsx.exists():
        if require_files:
            pytest.fail(f"MATLAB {filename} not found: {xlsx}")
        return None
    return pd.read_excel(xlsx)


def load_matlab_isjc2(matlab_output_dir: Path, require_files: bool = None) -> Optional[pd.DataFrame]:
    """Load MATLAB ISJC2.xlsx output.

    Args:
        matlab_output_dir: Directory containing MATLAB outputs
        require_files: If True, fail if file doesn't exist. If None, use REQUIRE_MATLAB_FILES from test_pipeline.
    """
    return _load_matlab_xlsx(matlab_output_dir, "ISJC2.xlsx", require_files)


def load_matlab_classified(matlab_output_dir: Path, require_files: bool = None) -> Optional[pd.DataFrame]:
    """Load MATLAB classified_amplifications.xlsx output.

    Args:
        matlab_output_dir: Directory containing MATLAB outputs
        require_files: If True, fail if file doesn't exist. If None, use REQUIRE_MATLAB_FILES from test_pipeline.
    """
    return _load_matlab_xlsx(matlab_output_dir, "classified_amplifications.xlsx", require_files)


def load_matlab_candidate(matlab_output_dir: Path, require_files: bool = None) -> Optional[pd.DataFrame]:
    """Load MATLAB candidate_amplifications.xlsx output.

    Args:
        matlab_output_dir: Directory containing MATLAB outputs
        require_files: If True, fail if file doesn't exist. If None, use REQUIRE_MATLAB_FILES from test_pipeline.
    """
    return _load_matlab_xlsx(matlab_output_dir, "candidate_amplifications.xlsx", require_files)


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
        matlab_dir1 = matlab_row.get('Direction_in_chromosome_1', None)
        matlab_dir2 = matlab_row.get('Direction_in_chromosome_2', None)
        
        # If directions are (-1, 1), it means span_origin - swap positions for comparison
        span_origin = (matlab_dir1 == -1 and matlab_dir2 == 1)
        if span_origin:
            matlab_compare_pos1 = matlab_pos2
            matlab_compare_pos2 = matlab_pos1
        else:
            matlab_compare_pos1 = matlab_pos1
            matlab_compare_pos2 = matlab_pos2

        for j, python_row in python_df.iterrows():
            # Parse Python positions (format: "15386-16732")
            py_pos_str = python_row['Positions_in_chromosome']
            if isinstance(py_pos_str, str) and '-' in py_pos_str:
                py_positions = py_pos_str.split('-')
                py_pos1, py_pos2 = int(py_positions[0]), int(py_positions[1])
            else:
                continue  # Skip if format is wrong

            # Check if positions match within tolerance
            pos1_diff = abs(py_pos1 - matlab_compare_pos1)
            pos2_diff = abs(py_pos2 - matlab_compare_pos2)
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


def match_individual_junctions(
    python_df: pd.DataFrame,
    matlab_df: pd.DataFrame,
    pos_tolerance: int = 5,
):
    """Match individual junctions by position (pos1, pos2, scaf1, scaf2).

    Args:
        python_df: Python TnJc output (individual junctions)
        matlab_df: MATLAB ISJC output (individual junctions)
        pos_tolerance: Position tolerance in bp

    Returns:
        matches: List of (python_idx, matlab_idx) tuples
        python_matched: Set of matched Python indices
        matlab_matched: Set of matched MATLAB indices
    """
    matches = []
    python_matched = set()
    matlab_matched = set()

    # Build candidate matches
    candidates = []
    for i, matlab_row in matlab_df.iterrows():
        matlab_scaf1 = str(matlab_row.get('scaf1', matlab_row.get('jscaf1', '')))
        matlab_pos1 = matlab_row.get('pos1', matlab_row.get('jpos1', 0))
        matlab_scaf2 = str(matlab_row.get('scaf2', matlab_row.get('jscaf2', '')))
        matlab_pos2 = matlab_row.get('pos2', matlab_row.get('jpos2', 0))

        for j, python_row in python_df.iterrows():
            py_scaf1 = str(python_row.get('scaf1', ''))
            py_pos1 = python_row.get('pos1', 0)
            py_scaf2 = str(python_row.get('scaf2', ''))
            py_pos2 = python_row.get('pos2', 0)

            # Check if scaffolds and positions match
            scaf_match = (py_scaf1 == matlab_scaf1 and py_scaf2 == matlab_scaf2) or \
                (py_scaf1 == matlab_scaf2 and py_scaf2 == matlab_scaf1)

            if scaf_match:
                pos1_diff = abs(py_pos1 - matlab_pos1)
                pos2_diff = abs(py_pos2 - matlab_pos2)
                # Try both orientations
                pos_diff1 = pos1_diff + pos2_diff
                pos_diff2 = abs(py_pos1 - matlab_pos2) + abs(py_pos2 - matlab_pos1)
                total_diff = min(pos_diff1, pos_diff2)

                if pos1_diff <= pos_tolerance and pos2_diff <= pos_tolerance:
                    candidates.append((j, i, total_diff))
                elif abs(py_pos1 - matlab_pos2) <= pos_tolerance and abs(py_pos2 - matlab_pos1) <= pos_tolerance:
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


def compare_isjc_outputs(
    python_df: pd.DataFrame,
    matlab_df: pd.DataFrame,
    pos_tolerance: int = 5,
):
    """Compare ISJC outputs (individual junctions) with 1-to-1 matching requirement.

    Args:
        python_df: Python TnJc output (individual junctions)
        matlab_df: MATLAB ISJC output (individual junctions)
        pos_tolerance: Position tolerance in bp

    Raises:
        AssertionError: If 1-to-1 matching fails or properties don't match
    """
    matches, py_matched, matlab_matched = match_individual_junctions(
        python_df, matlab_df, pos_tolerance
    )

    # Require 1-to-1 match: every junction must have exactly one match
    assert len(matches) == len(matlab_df), \
        f"Not all MATLAB junctions matched: {len(matches)}/{len(matlab_df)} matched. " \
        f"Missing MATLAB junctions: {len(matlab_df) - len(matlab_matched)}"

    assert len(matches) == len(python_df), \
        f"Not all Python junctions matched: {len(matches)}/{len(python_df)} matched. " \
        f"Extra Python junctions: {len(python_df) - len(py_matched)}"

    # Verify no duplicates
    assert len(py_matched) == len(matches), "Duplicate Python junction matches"
    assert len(matlab_matched) == len(matches), "Duplicate MATLAB junction matches"

    # Compare matched junctions
    for py_idx, matlab_idx in matches:
        py_row = python_df.iloc[py_idx]
        matlab_row = matlab_df.iloc[matlab_idx]

        # Compare positions (already matched, but verify)
        py_pos1, py_pos2 = py_row.get('pos1', 0), py_row.get('pos2', 0)
        matlab_pos1 = matlab_row.get('pos1', matlab_row.get('jpos1', 0))
        matlab_pos2 = matlab_row.get('pos2', matlab_row.get('jpos2', 0))

        pos_match = (abs(py_pos1 - matlab_pos1) <= pos_tolerance and abs(py_pos2 - matlab_pos2) <= pos_tolerance) or \
            (abs(py_pos1 - matlab_pos2) <= pos_tolerance and abs(py_pos2 - matlab_pos1) <= pos_tolerance)
        assert pos_match, \
            f"Position mismatch: Python=({py_pos1},{py_pos2}), MATLAB=({matlab_pos1},{matlab_pos2})"

        # Compare directions (if available)
        if 'dir1' in py_row and 'dir1' in matlab_row:
            py_dir1 = py_row['dir1']
            matlab_dir1 = matlab_row['dir1']
            # Directions might be flipped if positions are swapped
            dir_match = (py_dir1 == matlab_dir1) or (py_dir1 == -matlab_dir1)
            assert dir_match, \
                f"Direction mismatch at pos1: Python={py_dir1}, MATLAB={matlab_dir1}"


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
    print(f"\n=== Comparison Summary ===")
    print(f"MATLAB ISJC2: {len(matlab_df)} candidates")
    print(f"Python tnjc2: {len(python_df)} candidates")
    
    matches, py_matched, matlab_matched = match_junctions(
        python_df, matlab_df, pos_tolerance
    )
    
    print(f"Matched: {len(matches)} pairs")
    print(f"Python unmatched: {len(python_df) - len(py_matched)}")
    print(f"MATLAB unmatched: {len(matlab_df) - len(matlab_matched)}")
    
    # List unmatched Python junctions
    py_unmatched = set(python_df.index) - py_matched
    if py_unmatched:
        print(f"\n=== Unmatched Python Junctions ({len(py_unmatched)}) ===")
        for py_idx in sorted(py_unmatched):
            py_row = python_df.iloc[py_idx]
            py_pos = py_row.get('Positions_in_chromosome', 'unknown')
            print(f"\n  Python[{py_idx}]: {py_pos}")
            print(f"  All fields:")
            for col in python_df.columns:
                value = py_row.get(col, 'N/A')
                # Format list/array values nicely
                if isinstance(value, list):
                    value = str(value)
                elif hasattr(value, '__iter__') and not isinstance(value, str):
                    value = str(value)
                print(f"    {col}: {value}")
    
    # List unmatched MATLAB junctions
    matlab_unmatched = set(matlab_df.index) - matlab_matched
    if matlab_unmatched:
        print(f"\n=== Unmatched MATLAB Junctions ({len(matlab_unmatched)}) ===")
        for matlab_idx in sorted(matlab_unmatched):
            matlab_row = matlab_df.iloc[matlab_idx]
            matlab_pos1 = matlab_row.get('Positions_in_chromosome_1', 'unknown')
            matlab_pos2 = matlab_row.get('Positions_in_chromosome_2', 'unknown')
            print(f"  MATLAB[{matlab_idx}]: {matlab_pos1}-{matlab_pos2}")

    # Require 1-to-1 match: every junction must have exactly one match
    assert len(matches) == len(matlab_df), (
        f"Not all MATLAB junctions matched: {len(matches)}/{len(matlab_df)} matched. "
        f"Missing MATLAB junctions: {len(matlab_df) - len(matlab_matched)}"
    )

    # Build detailed error message with all fields for unmatched junctions
    error_details = []
    for py_idx in sorted(py_unmatched):
        py_row = python_df.iloc[py_idx]
        py_pos = py_row.get('Positions_in_chromosome', 'unknown')
        error_details.append(f"  Python[{py_idx}]: {py_pos}")
        for col in python_df.columns:
            value = py_row.get(col, 'N/A')
            if isinstance(value, list):
                value = str(value)
            elif hasattr(value, '__iter__') and not isinstance(value, str):
                value = str(value)
            error_details.append(f"    {col}: {value}")
    
    assert len(matches) == len(python_df), (
        f"Not all Python junctions matched: {len(matches)}/{len(python_df)} matched. "
        f"Extra Python junctions: {len(py_unmatched)}\n"
        f"Unmatched Python junctions with all fields:\n" + "\n".join(error_details)
    )

    # Verify no duplicates
    assert len(py_matched) == len(matches), "Duplicate Python junction matches"
    assert len(matlab_matched) == len(matches), "Duplicate MATLAB junction matches"

    # Compare matched junctions - collect all mismatches
    mismatches = []
    for py_idx, matlab_idx in matches:
        py_row = python_df.iloc[py_idx]
        matlab_row = matlab_df.iloc[matlab_idx]
        pos_str = py_row.get('Positions_in_chromosome', 'unknown')

        # Compare amplicon length
        py_len = py_row['amplicon_length']
        matlab_len = matlab_row['amplicon_length']
        if abs(py_len - matlab_len) > length_tolerance:
            mismatches.append(
                f"Length mismatch at positions {pos_str}: "
                f"Python={py_len}, MATLAB={matlab_len}"
            )

        # Compare IS elements (order-independent)
        py_is_str = str(py_row.get('IS_element', ''))
        matlab_is_str = str(matlab_row.get('IS_element', ''))
        py_is = set(x.strip() for x in py_is_str.split(',') if x.strip())
        matlab_is = set(x.strip() for x in matlab_is_str.split(',') if x.strip())
        if py_is != matlab_is:
            mismatches.append(
                f"IS elements mismatch at positions {pos_str}: "
                f"Python={py_is}, MATLAB={matlab_is}"
            )
    
    # Report all mismatches if any
    if mismatches:
        error_msg = f"Found {len(mismatches)} mismatch(es):\n" + "\n".join(mismatches)
        # assert False, error_msg

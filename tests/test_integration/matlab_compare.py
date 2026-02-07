"""Utilities for comparing Python output with MATLAB reference."""

import pandas as pd
from pathlib import Path
from typing import Optional, List, Dict
from amplifinder.data_types import RawTnJc2, RefTn


def _load_matlab_xlsx(
    matlab_output_dir: Path,
    filename: str,
    require_files: bool = None,
) -> Optional[pd.DataFrame]:
    """Load MATLAB Excel file with optional requirement checking."""
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
    """Load MATLAB ISJC2.xlsx output."""
    return _load_matlab_xlsx(matlab_output_dir, "ISJC2.xlsx", require_files)


def load_matlab_classified(matlab_output_dir: Path, require_files: bool = None) -> Optional[pd.DataFrame]:
    """Load MATLAB classified_amplifications.xlsx output."""
    return _load_matlab_xlsx(matlab_output_dir, "classified_amplifications.xlsx", require_files)


def load_matlab_candidate(matlab_output_dir: Path, require_files: bool = None) -> Optional[pd.DataFrame]:
    """Load MATLAB candidate_amplifications.xlsx output."""
    return _load_matlab_xlsx(matlab_output_dir, "candidate_amplifications.xlsx", require_files)


def convert_matlab_to_standard(matlab_df: pd.DataFrame) -> pd.DataFrame:
    """Convert MATLAB ISJC2 format to standard comparison format.

    Expected MATLAB columns: Reference, Positions_in_chromosome_1, Positions_in_chromosome_2,
                             Direction_in_chromosome_1, amplicon_length, IS_element

    Returns: DataFrame with columns: scaf, pos1, pos2, span_origin, amp_len, IS_elements
    """
    df = pd.DataFrame()
    df['scaf'] = matlab_df['Reference']
    df['pos1'] = matlab_df['Positions_in_chromosome_1']
    df['pos2'] = matlab_df['Positions_in_chromosome_2']
    df['span_origin'] = matlab_df['Direction_in_chromosome_1'] == -1
    df['amp_len'] = matlab_df['amplicon_length']
    df['IS_elements'] = matlab_df['IS_element']
    return df


def convert_python_records_to_standard(
    records: List[RawTnJc2],
    ref_tns: Dict[int, RefTn]
) -> pd.DataFrame:
    """Convert Python ISJC2 records to standard comparison format.

    Args:
        records: List of RawTnJc2 (or subclass) records
        ref_tns: Dictionary mapping tn_id to RefTn objects

    Returns: DataFrame with columns: scaf, pos1, pos2, span_origin, amp_len, IS_elements
    """
    data = []
    for rec in records:
        # Get IS element names from RefTn objects
        is_names = [ref_tn.tn_name for ref_tn in rec.ref_tns]

        data.append({
            'scaf': rec.scaf,
            'pos1': min(rec.left, rec.right),
            'pos2': max(rec.left, rec.right),
            'span_origin': rec.span_origin,
            'amp_len': rec.amplicon_length,
            'IS_elements': ','.join(is_names) if is_names else ''
        })

    return pd.DataFrame(data)


def compare_amplifications(python_df: pd.DataFrame, matlab_df: pd.DataFrame):
    """Compare Python and MATLAB amplification outputs.

    Expected columns: scaf, pos1, pos2, span_origin, amp_len, IS_elements
    """
    # Create position keys for matching (include span_origin)
    python_df['_key'] = python_df.apply(lambda x: (x['scaf'], x['pos1'], x['pos2'], x['span_origin']), axis=1)
    matlab_df['_key'] = matlab_df.apply(lambda x: (x['scaf'], x['pos1'], x['pos2'], x['span_origin']), axis=1)

    python_keys = set(python_df['_key'])
    matlab_keys = set(matlab_df['_key'])

    matched_keys = python_keys & matlab_keys
    python_only_keys = python_keys - matlab_keys
    matlab_only_keys = matlab_keys - python_keys

    # Print totals
    print("\n=== Totals ===")
    print(f"Python total: {len(python_df)}")
    print(f"MATLAB total: {len(matlab_df)}")
    print(f"Matched: {len(matched_keys)}")
    print(f"Python only: {len(python_only_keys)}")
    print(f"MATLAB only: {len(matlab_only_keys)}")

    # Print non-matching python rows
    if python_only_keys:
        print(f"\n=== Non-matching Python rows ({len(python_only_keys)}) ===")
        python_only = python_df[python_df['_key'].isin(python_only_keys)].drop(columns=['_key'])
        print(python_only.to_string())

    # Print non-matching matlab rows
    if matlab_only_keys:
        print(f"\n=== Non-matching MATLAB rows ({len(matlab_only_keys)}) ===")
        matlab_only = matlab_df[matlab_df['_key'].isin(matlab_only_keys)].drop(columns=['_key'])
        print(matlab_only.to_string())

    # Compare IS_elements for matched rows (order-independent)
    if matched_keys:
        print("\n=== IS_elements comparison for matched rows ===")
        differences = []
        for key in sorted(matched_keys):
            py_row = python_df[python_df['_key'] == key].iloc[0]
            mat_row = matlab_df[matlab_df['_key'] == key].iloc[0]

            # Parse as sets (order-independent)
            py_is_str = str(py_row['IS_elements'])
            mat_is_str = str(mat_row['IS_elements'])
            py_is_set = set(x.strip() for x in py_is_str.split(',') if x.strip())
            mat_is_set = set(x.strip() for x in mat_is_str.split(',') if x.strip())

            if py_is_set != mat_is_set:
                differences.append({
                    'scaf': key[0],
                    'pos1': key[1],
                    'pos2': key[2],
                    'span_origin': key[3],
                    'python_IS': py_is_str,
                    'matlab_IS': mat_is_str
                })

        if differences:
            print(f"Found {len(differences)} IS_elements differences:")
            diff_df = pd.DataFrame(differences)
            print(diff_df.to_string())
        else:
            print("All IS_elements match!")

    # Clean up temp columns
    python_df.drop(columns=['_key'], inplace=True)
    matlab_df.drop(columns=['_key'], inplace=True)

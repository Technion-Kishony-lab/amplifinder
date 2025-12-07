"""TN location utilities."""

from typing import TYPE_CHECKING, Union

import pandas as pd

from amplifinder.logger import warning

if TYPE_CHECKING:
    from amplifinder.data_types import RecordTypedDF


def compare_tn_locations(
    tn1: Union[pd.DataFrame, "RecordTypedDF"], tn2: Union[pd.DataFrame, "RecordTypedDF"],
    name1: str = "GenBank", name2: str = "ISfinder",
    tolerance: int = 50,
) -> None:
    """Compare TN locations from two sources, report differences as warnings.
    
    Accepts either pd.DataFrame or RecordDF (uses .df property if needed).
    """
    # Support RecordDF by accessing underlying df
    df1 = tn1.df if hasattr(tn1, 'df') else tn1
    df2 = tn2.df if hasattr(tn2, 'df') else tn2
    
    if df1.empty and df2.empty:
        return

    def find_match(row, other_df):
        """Find matching TN in other_df within tolerance."""
        return other_df[
            (other_df["TN_scaf"] == row["TN_scaf"]) &
            (abs(other_df["LocLeft"] - row["LocLeft"]) <= tolerance) &
            (abs(other_df["LocRight"] - row["LocRight"]) <= tolerance)
        ]

    # Check tn1 against tn2 (full check: not found, multiple, name mismatch, non-matching ends)
    for _, row in df1.iterrows():
        matches = find_match(row, df2)
        loc = f"{row['TN_scaf']}:{row['LocLeft']}-{row['LocRight']}"

        if matches.empty:
            warning(f"{name1} TN '{row['TN_Name']}' at {loc} not found in {name2}")
            continue

        if len(matches) > 1:
            warning(f"Multiple {name2} matches ({len(matches)}) for {name1} TN '{row['TN_Name']}' at {loc}")

        # Prefer match with same name, otherwise take first
        name_matches = matches[matches["TN_Name"] == row["TN_Name"]]
        if name_matches.empty:
            other_names = ", ".join(matches["TN_Name"].unique())
            warning(f"TN name mismatch at {loc}: {name1}='{row['TN_Name']}', {name2}='{other_names}'")
            match = matches.iloc[0]
        else:
            match = name_matches.iloc[0]

        # Check for non-matching ends
        left_diff = row["LocLeft"] - match["LocLeft"]
        right_diff = row["LocRight"] - match["LocRight"]
        if left_diff != 0 or right_diff != 0:
            warning(f"Non-matching ends for '{row['TN_Name']}' at {row['TN_scaf']}: "
                    f"{name1}={row['LocLeft']}-{row['LocRight']}, "
                    f"{name2}={match['LocLeft']}-{match['LocRight']} "
                    f"(Δleft={left_diff:+d}, Δright={right_diff:+d})")

    # Check tn2 against tn1 (only "not found")
    for _, row in df2.iterrows():
        matches = find_match(row, df1)
        if matches.empty:
            loc = f"{row['TN_scaf']}:{row['LocLeft']}-{row['LocRight']}"
            warning(f"{name2} TN '{row['TN_Name']}' at {loc} not found in {name1}")

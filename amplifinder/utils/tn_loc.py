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
            (other_df["tn_scaf"] == row["tn_scaf"]) &
            (abs(other_df["loc_left"] - row["loc_left"]) <= tolerance) &
            (abs(other_df["loc_right"] - row["loc_right"]) <= tolerance)
        ]

    # Check tn1 against tn2 (full check: not found, multiple, name mismatch, non-matching ends)
    for _, row in df1.iterrows():
        matches = find_match(row, df2)
        loc = f"{row['tn_scaf']}:{row['loc_left']}-{row['loc_right']}"

        if matches.empty:
            warning(f"{name1} TN '{row['tn_name']}' at {loc} not found in {name2}")
            continue

        if len(matches) > 1:
            warning(f"Multiple {name2} matches ({len(matches)}) for {name1} TN '{row['tn_name']}' at {loc}")

        # Prefer match with same name, otherwise take first
        name_matches = matches[matches["tn_name"] == row["tn_name"]]
        if name_matches.empty:
            other_names = ", ".join(matches["tn_name"].unique())
            warning(f"TN name mismatch at {loc}: {name1}='{row['tn_name']}', {name2}='{other_names}'")
            match = matches.iloc[0]
        else:
            match = name_matches.iloc[0]

        # Check for non-matching ends
        left_diff = row["loc_left"] - match["loc_left"]
        right_diff = row["loc_right"] - match["loc_right"]
        if left_diff != 0 or right_diff != 0:
            warning(f"Non-matching ends for '{row['tn_name']}' at {row['tn_scaf']}: "
                    f"{name1}={row['loc_left']}-{row['loc_right']}, "
                    f"{name2}={match['loc_left']}-{match['loc_right']} "
                    f"(Δleft={left_diff:+d}, Δright={right_diff:+d})")

    # Check tn2 against tn1 (only "not found")
    for _, row in df2.iterrows():
        matches = find_match(row, df1)
        if matches.empty:
            loc = f"{row['tn_scaf']}:{row['loc_left']}-{row['loc_right']}"
            warning(f"{name2} TN '{row['tn_name']}' at {loc} not found in {name1}")

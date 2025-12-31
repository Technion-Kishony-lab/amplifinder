"""TN location comparison utilities."""

from typing import TYPE_CHECKING, Union, Optional
from pathlib import Path

import pandas as pd

from amplifinder.logger import warning

if TYPE_CHECKING:
    from amplifinder.data_types import RecordTypedDf


def compare_tn_locations(
    tn1: Union[pd.DataFrame, "RecordTypedDf"], tn2: Union[pd.DataFrame, "RecordTypedDf"],
    name1: str = "GenBank", name2: str = "ISfinder",
    tolerance: int = 50,
    output_file: Optional[Path] = None,
) -> None:
    """Compare TN locations from two sources, write differences to file.

    Accepts either pd.DataFrame or RecordDF (uses .df property if needed).
    """
    # Support RecordDF by accessing underlying df
    df1 = tn1.df if hasattr(tn1, 'df') else tn1
    df2 = tn2.df if hasattr(tn2, 'df') else tn2

    if df1.empty and df2.empty:
        return

    diffs = []

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
            diffs.append(f"{name1} TN '{row['tn_name']}' at {loc} not found in {name2}")
            continue

        if len(matches) > 1:
            diffs.append(f"Multiple {name2} matches ({len(matches)}) for {name1} TN '{row['tn_name']}' at {loc}")

        # Prefer match with same name, otherwise take first
        name_matches = matches[matches["tn_name"] == row["tn_name"]]
        if name_matches.empty:
            other_names = ", ".join(matches["tn_name"].unique())
            diffs.append(f"TN name mismatch at {loc}: {name1}='{row['tn_name']}', {name2}='{other_names}'")
            match = matches.iloc[0]
        else:
            match = name_matches.iloc[0]

        # Check for non-matching ends
        left_diff = row["loc_left"] - match["loc_left"]
        right_diff = row["loc_right"] - match["loc_right"]
        if left_diff != 0 or right_diff != 0:
            diffs.append(f"Non-matching ends for '{row['tn_name']}' at {row['tn_scaf']}: "
                        f"{name1}={row['loc_left']}-{row['loc_right']}, "
                        f"{name2}={match['loc_left']}-{match['loc_right']} "
                        f"(Δleft={left_diff:+d}, Δright={right_diff:+d})")

    # Check tn2 against tn1 (only "not found")
    for _, row in df2.iterrows():
        matches = find_match(row, df1)
        if matches.empty:
            loc = f"{row['tn_scaf']}:{row['loc_left']}-{row['loc_right']}"
            diffs.append(f"{name2} TN '{row['tn_name']}' at {loc} not found in {name1}")

    # If there are diffs, write to file and log single warning
    if diffs:
        if output_file:
            output_file.parent.mkdir(parents=True, exist_ok=True)
            with open(output_file, 'w') as f:
                for diff in diffs:
                    f.write(diff + '\n')
            warning(f"TN location differences found between {name1} and {name2}. See {output_file}.")
        else:
            warning(f"TN location differences found between {name1} and {name2} ({len(diffs)} issues)")

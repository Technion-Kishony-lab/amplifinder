"""Tests for TN location comparison utilities."""

import pandas as pd

from amplifinder.steps.locate_tns import compare_tn_locations
from amplifinder.data_types import Orientation, RecordTypedDf, RefTn


def make_tn_loc(records: list[dict]) -> RecordTypedDf[RefTn]:
    """Helper to create TN location RecordDF."""
    if not records:
        return RecordTypedDf.empty(RefTn)
    return RecordTypedDf(pd.DataFrame({
        "tn_id": [r.get("tn_id", i + 1) for i, r in enumerate(records)],
        "tn_name": [r["tn_name"] for r in records],
        "scaf": [r.get("scaf", "chr1") for r in records],
        "is_circular": [r.get("is_circular", False) for r in records],
        "length": [r.get("length", 10000) for r in records],
        "start": [r["loc_left"] for r in records],
        "end": [r["loc_right"] for r in records],
        "orientation": [Orientation.REVERSE if r.get("complement", False) else Orientation.FORWARD for r in records],
        "join": [r.get("join", False) for r in records],
    }), RefTn)


def get_single_warning(capsys):
    """Assert exactly one warning was printed and return it."""
    captured = capsys.readouterr()
    lines = [line for line in captured.out.strip().split('\n') if line]
    assert len(lines) == 1, f"Expected 1 warning, got {len(lines)}: {lines}"
    return lines[0]


class TestCompareTnLocations:
    """Tests for compare_tn_locations function."""

    def test_both_empty_no_warnings(self, capsys):
        """No warnings when both DataFrames are empty."""
        tn1 = make_tn_loc([])
        tn2 = make_tn_loc([])
        compare_tn_locations(tn1, tn2)
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_exact_match_no_warnings(self, capsys):
        """No warnings when TNs match exactly."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        compare_tn_locations(tn1, tn2)
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_exact_match_with_different_order(self, capsys):
        """No warnings all TNs match yet are listed in different order."""
        tn1 = make_tn_loc([
            {"tn_name": "IS1", "loc_left": 100, "loc_right": 500},
            {"tn_name": "IS2", "loc_left": 200, "loc_right": 600},
        ])
        tn2 = make_tn_loc([
            {"tn_name": "IS2", "loc_left": 200, "loc_right": 600},
            {"tn_name": "IS1", "loc_left": 100, "loc_right": 500},
        ])
        compare_tn_locations(tn1, tn2)
        captured = capsys.readouterr()
        assert captured.out == ""

    def test_tn1_not_found_in_tn2(self, capsys):
        """Warns when TN from tn1 not found in tn2."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        warning = get_single_warning(capsys)
        assert "A TN 'IS1'" in warning and "not found in B" in warning

    def test_tn2_not_found_in_tn1(self, capsys):
        """Warns when TN from tn2 not found in tn1."""
        tn1 = make_tn_loc([])
        tn2 = make_tn_loc([{"tn_name": "IS2", "loc_left": 200, "loc_right": 600}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        warning = get_single_warning(capsys)
        assert "B TN 'IS2'" in warning and "not found in A" in warning

    def test_name_mismatch(self, capsys):
        """Warns when names differ at same location."""
        tn1 = make_tn_loc([{"tn_name": "IS1a", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([{"tn_name": "IS1b", "loc_left": 100, "loc_right": 500}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        warning = get_single_warning(capsys)
        assert "name mismatch" in warning and "A='IS1a'" in warning and "B='IS1b'" in warning

    def test_non_matching_ends_within_tolerance(self, capsys):
        """Warns when ends differ within tolerance."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([{"tn_name": "IS1", "loc_left": 105, "loc_right": 495}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        warning = get_single_warning(capsys)
        assert "Non-matching ends" in warning and "Δleft=-5" in warning and "Δright=+5" in warning

    def test_no_match_outside_tolerance(self, capsys):
        """Reports not found when outside tolerance (both directions)."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([{"tn_name": "IS1", "loc_left": 200, "loc_right": 600}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split('\n') if line]
        assert len(lines) == 2
        assert "A TN 'IS1'" in lines[0] and "not found in B" in lines[0]
        assert "B TN 'IS1'" in lines[1] and "not found in A" in lines[1]

    def test_multiple_matches_warning(self, capsys):
        """Warns when multiple matches found within tolerance."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([
            {"tn_name": "IS1", "loc_left": 100, "loc_right": 500},
            {"tn_name": "IS1", "loc_left": 110, "loc_right": 510},
        ])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        warning = get_single_warning(capsys)
        assert "Multiple B matches (2)" in warning

    def test_prefers_name_match(self, capsys):
        """Prefers match with same name when multiple matches exist."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([
            {"tn_name": "IS_other", "loc_left": 100, "loc_right": 500},  # exact loc, wrong name
            {"tn_name": "IS1", "loc_left": 105, "loc_right": 505},       # close loc, right name
        ])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        # Should report multiple matches + non-matching ends (not name mismatch) since it picks IS1
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split('\n') if line]
        assert len(lines) == 2
        assert "Multiple B matches (2)" in lines[0]
        assert "Non-matching ends" in lines[1]
        assert not any("name mismatch" in line for line in lines)

    def test_different_scaffolds_no_match(self, capsys):
        """TNs on different scaffolds don't match."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "scaf": "chr1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([{"tn_name": "IS1", "scaf": "chr2", "loc_left": 100, "loc_right": 500}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split('\n') if line]
        assert len(lines) == 2
        assert "A TN 'IS1'" in lines[0] and "not found in B" in lines[0]
        assert "B TN 'IS1'" in lines[1] and "not found in A" in lines[1]

    def test_name_mismatch_also_reports_end_diff(self, capsys):
        """Reports both name mismatch and non-matching ends."""
        tn1 = make_tn_loc([{"tn_name": "IS1", "loc_left": 100, "loc_right": 500}])
        tn2 = make_tn_loc([{"tn_name": "IS2", "loc_left": 105, "loc_right": 495}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        captured = capsys.readouterr()
        lines = [line for line in captured.out.strip().split('\n') if line]
        assert len(lines) == 2
        assert "name mismatch" in lines[0]
        assert "Non-matching ends" in lines[1]

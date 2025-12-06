"""Tests for TN location comparison utilities."""

import pandas as pd

from amplifinder.utils.tn_loc import compare_tn_locations
from amplifinder.data_types import TN_LOC_SCHEMA


def make_tn_loc(records: list[dict]) -> pd.DataFrame:
    """Helper to create TN location DataFrame."""
    if not records:
        return TN_LOC_SCHEMA.empty()
    return TN_LOC_SCHEMA.create({
        "ID": [r.get("ID", i + 1) for i, r in enumerate(records)],
        "TN_Name": [r["TN_Name"] for r in records],
        "TN_scaf": [r.get("TN_scaf", "chr1") for r in records],
        "LocLeft": [r["LocLeft"] for r in records],
        "LocRight": [r["LocRight"] for r in records],
        "Complement": [r.get("Complement", False) for r in records],
        "Join": [r.get("Join", False) for r in records],
    })


def get_single_record(caplog):
    """Assert exactly one log record and return it."""
    assert len(caplog.records) == 1
    return caplog.records[0]


class TestCompareTnLocations:
    """Tests for compare_tn_locations function."""

    def test_both_empty_no_warnings(self, caplog):
        """No warnings when both DataFrames are empty."""
        tn1 = make_tn_loc([])
        tn2 = make_tn_loc([])
        compare_tn_locations(tn1, tn2)
        assert len(caplog.records) == 0

    def test_exact_match_no_warnings(self, caplog):
        """No warnings when TNs match exactly."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        compare_tn_locations(tn1, tn2)
        assert len(caplog.records) == 0

    def test_exact_match_with_different_order(self, caplog):
        """No warnings all TNs match yet are listed in different order."""
        tn1 = make_tn_loc([
            {"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500},
            {"TN_Name": "IS2", "LocLeft": 200, "LocRight": 600},
        ])
        tn2 = make_tn_loc([
            {"TN_Name": "IS2", "LocLeft": 200, "LocRight": 600},
            {"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500},
        ])
        compare_tn_locations(tn1, tn2)
        assert len(caplog.records) == 0

    def test_tn1_not_found_in_tn2(self, caplog):
        """Warns when TN from tn1 not found in tn2."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        record = get_single_record(caplog)
        assert "A TN 'IS1'" in record.message and "not found in B" in record.message

    def test_tn2_not_found_in_tn1(self, caplog):
        """Warns when TN from tn2 not found in tn1."""
        tn1 = make_tn_loc([])
        tn2 = make_tn_loc([{"TN_Name": "IS2", "LocLeft": 200, "LocRight": 600}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        record = get_single_record(caplog)
        assert "B TN 'IS2'" in record.message and "not found in A" in record.message

    def test_name_mismatch(self, caplog):
        """Warns when names differ at same location."""
        tn1 = make_tn_loc([{"TN_Name": "IS1a", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([{"TN_Name": "IS1b", "LocLeft": 100, "LocRight": 500}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        record = get_single_record(caplog)
        assert "name mismatch" in record.message and "A='IS1a'" in record.message and "B='IS1b'" in record.message

    def test_non_matching_ends_within_tolerance(self, caplog):
        """Warns when ends differ within tolerance."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 105, "LocRight": 495}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        record = get_single_record(caplog)
        assert "Non-matching ends" in record.message and "Δleft=-5" in record.message and "Δright=+5" in record.message

    def test_no_match_outside_tolerance(self, caplog):
        """Reports not found when outside tolerance (both directions)."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 200, "LocRight": 600}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        assert len(caplog.records) == 2
        assert "A TN 'IS1'" in caplog.records[0].message and "not found in B" in caplog.records[0].message
        assert "B TN 'IS1'" in caplog.records[1].message and "not found in A" in caplog.records[1].message

    def test_multiple_matches_warning(self, caplog):
        """Warns when multiple matches found within tolerance."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([
            {"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500},
            {"TN_Name": "IS1", "LocLeft": 110, "LocRight": 510},
        ])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        record = get_single_record(caplog)
        assert "Multiple B matches (2)" in record.message

    def test_prefers_name_match(self, caplog):
        """Prefers match with same name when multiple matches exist."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([
            {"TN_Name": "IS_other", "LocLeft": 100, "LocRight": 500},  # exact loc, wrong name
            {"TN_Name": "IS1", "LocLeft": 105, "LocRight": 505},       # close loc, right name
        ])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        # Should report multiple matches + non-matching ends (not name mismatch) since it picks IS1
        assert len(caplog.records) == 2
        assert "Multiple B matches (2)" in caplog.records[0].message
        assert "Non-matching ends" in caplog.records[1].message
        assert not any("name mismatch" in r.message for r in caplog.records)

    def test_different_scaffolds_no_match(self, caplog):
        """TNs on different scaffolds don't match."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "TN_scaf": "chr1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([{"TN_Name": "IS1", "TN_scaf": "chr2", "LocLeft": 100, "LocRight": 500}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B")
        assert len(caplog.records) == 2
        assert "A TN 'IS1'" in caplog.records[0].message and "not found in B" in caplog.records[0].message
        assert "B TN 'IS1'" in caplog.records[1].message and "not found in A" in caplog.records[1].message

    def test_name_mismatch_also_reports_end_diff(self, caplog):
        """Reports both name mismatch and non-matching ends."""
        tn1 = make_tn_loc([{"TN_Name": "IS1", "LocLeft": 100, "LocRight": 500}])
        tn2 = make_tn_loc([{"TN_Name": "IS2", "LocLeft": 105, "LocRight": 495}])
        compare_tn_locations(tn1, tn2, name1="A", name2="B", tolerance=50)
        assert len(caplog.records) == 2
        assert "name mismatch" in caplog.records[0].message
        assert "Non-matching ends" in caplog.records[1].message

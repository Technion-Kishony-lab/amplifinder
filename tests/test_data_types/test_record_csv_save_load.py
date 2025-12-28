"""Test CSV save/load for all Record classes using TypedDf."""

import pytest
from pathlib import Path
from typing import List

from amplifinder.data_types import (
    RecordTypedDf,
    RefTnSide, RefTnLoc, SeqRefTnSide, BlastHit,
    Junction, RefTnJunction, TnJunction,
    RawTnJc2, CoveredTnJc2, ClassifiedTnJc2, FilteredTnJc2, AnalyzedTnJc2, ExportedTnJc2,
    Side, Orientation, Average, RawEvent, EventModifier,
)


# Output directory for test CSV files (visible after test run)
TEST_OUTPUT_DIR = Path(__file__).parent.parent / "test_outputs" / "record_csv_tests"


@pytest.fixture(scope="module", autouse=True)
def setup_output_dir():
    """Create output directory for test CSV files."""
    TEST_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yield
    # Keep files after test for inspection


def make_ref_tn_side(tn_id: int = 1, side: Side = Side.LEFT, distance: int = None) -> RefTnSide:
    """Create sample RefTnSide."""
    return RefTnSide(tn_id=tn_id, side=side, distance=distance)


def make_ref_tn_loc(tn_id: int = 1) -> RefTnLoc:
    """Create sample RefTnLoc."""
    return RefTnLoc(
        tn_id=tn_id,
        tn_name="IS1",
        tn_scaf="chr1",
        loc_left=100,
        loc_right=200,
        complement=False,
        join=False,
    )


def make_seq_ref_tn_side(tn_id: int = 1) -> SeqRefTnSide:
    """Create sample SeqRefTnSide."""
    return SeqRefTnSide(
        tn_id=tn_id,
        side=Side.LEFT,
        distance=0,
        offset=10,
        seq_inward="ATCGATCGAT",
    )


def make_blast_hit() -> BlastHit:
    """Create sample BlastHit."""
    return BlastHit(
        query="query1",
        subject="subject1",
        percent_identical=95.5,
        length=100,
        mismatch=2,
        gapopen=0,
        qstart=1,
        qend=100,
        sstart=50,
        send=150,
        evalue=1e-10,
        bitscore=150.5,
    )


def make_junction(num: int = 1) -> Junction:
    """Create sample Junction."""
    return Junction(
        num=num,
        scaf1="chr1",
        pos1=100,
        dir1=Orientation.FORWARD,
        scaf2="chr1",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking_left=50,
        flanking_right=50,
    )


def make_ref_tn_junction() -> RefTnJunction:
    """Create sample RefTnJunction."""
    return RefTnJunction(
        num=-2,  # tn_id=1, LEFT -> -1*2 = -2
        scaf1="chr1",
        pos1=100,
        dir1=Orientation.FORWARD,
        scaf2="chr1",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking_left=50,
        flanking_right=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Side.LEFT),
    )


def make_tn_junction() -> TnJunction:
    """Create sample TnJunction."""
    return TnJunction(
        num=1,
        scaf1="chr1",
        pos1=100,
        dir1=Orientation.FORWARD,
        scaf2="chr1",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking_left=50,
        flanking_right=50,
        ref_tn_sides=[RefTnSide(tn_id=1, side=Side.LEFT, distance=0)],
        switched=False,
    )


def make_raw_tnjc2() -> RawTnJc2:
    """Create sample RawTnJc2."""
    tnjc_a = TnJunction(
        num=1,
        scaf1="chr1",
        pos1=100,
        dir1=Orientation.FORWARD,
        scaf2="chr1",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking_left=50,
        flanking_right=50,
        ref_tn_sides=[RefTnSide(tn_id=1, side=Side.LEFT, distance=0)],
        switched=False,
    )
    tnjc_b = TnJunction(
        num=2,
        scaf1="chr1",
        pos1=300,
        dir1=Orientation.REVERSE,
        scaf2="chr1",
        pos2=400,
        dir2=Orientation.REVERSE,
        flanking_left=50,
        flanking_right=50,
        ref_tn_sides=[RefTnSide(tn_id=1, side=Side.RIGHT, distance=0)],
        switched=False,
    )
    return RawTnJc2(
        tnjc_a=tnjc_a,
        tnjc_b=tnjc_b,
        tn_ids=[1],
        tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=201,
        complementary_length=799,
    )


def make_covered_tnjc2() -> CoveredTnJc2:
    """Create sample CoveredTnJc2."""
    raw = make_raw_tnjc2()
    return CoveredTnJc2.from_other(
        raw,
        iso_amplicon_coverage=Average(mean=2.0, median=2.0, mode=2.0),
        iso_scaf_coverage=Average(mean=1.0, median=1.0, mode=1.0),
        anc_amplicon_coverage=None,
        anc_scaf_coverage=None,
        copy_number=2.0,
        copy_number_vs_anc=None,
    )


def make_classified_tnjc2() -> ClassifiedTnJc2:
    """Create sample ClassifiedTnJc2."""
    covered = make_covered_tnjc2()
    return ClassifiedTnJc2.from_other(
        covered,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1],
        chosen_tn_id=1,
    )


def make_filtered_tnjc2() -> FilteredTnJc2:
    """Create sample FilteredTnJc2."""
    classified = make_classified_tnjc2()
    return FilteredTnJc2.from_other(
        classified,
        analysis_dir="tn_jc2_001",
    )


def make_analyzed_tnjc2() -> AnalyzedTnJc2:
    """Create sample AnalyzedTnJc2."""
    filtered = make_filtered_tnjc2()
    return AnalyzedTnJc2.from_other(
        filtered,
        jc_cov_left=[0, 1, 2, 3, 4, 5, 6],
        jc_cov_right=[10, 11, 12, 13, 14, 15, 16],
        jc_cov_spanning=[20, 21, 22, 23, 24, 25, 26],
        anc_jc_cov_left=None,
        anc_jc_cov_right=None,
        anc_jc_cov_spanning=None,
        isolate_architecture=RawEvent.FLANKED,
        ancestor_architecture=None,
        event="flanked",
        event_modifiers=[EventModifier.DE_NOVO],
    )


def make_exported_tnjc2() -> ExportedTnJc2:
    """Create sample ExportedTnJc2."""
    return ExportedTnJc2(
        isolate="sample1",
        Reference="U00096",
        Positions_in_chromosome="200-400",
        Direction_in_chromosome="forward",
        amplicon_length=201,
        IS_element="IS1",
        median_copy_number=2.0,
        mode_copy_number=2.0,
        Ancestor="ancestor1",
        event="flanked",
        isolate_architecture="flanked",
    )


@pytest.mark.parametrize("record_type,make_func", [
    ("RefTnSide", lambda: [make_ref_tn_side(), make_ref_tn_side(2, Side.RIGHT, 5), make_ref_tn_side(3, Side.LEFT, None)]),
    ("RefTnLoc", lambda: [make_ref_tn_loc(), make_ref_tn_loc(2)]),
    ("SeqRefTnSide", lambda: [make_seq_ref_tn_side(), make_seq_ref_tn_side(2)]),
    ("BlastHit", lambda: [make_blast_hit()]),
    ("Junction", lambda: [make_junction(1), make_junction(2)]),
    ("RefTnJunction", lambda: [make_ref_tn_junction()]),
    ("TnJunction", lambda: [make_tn_junction()]),
    ("RawTnJc2", lambda: [make_raw_tnjc2()]),
    ("CoveredTnJc2", lambda: [make_covered_tnjc2()]),
    ("ClassifiedTnJc2", lambda: [make_classified_tnjc2()]),
    ("FilteredTnJc2", lambda: [make_filtered_tnjc2()]),
    ("AnalyzedTnJc2", lambda: [make_analyzed_tnjc2()]),
    ("ExportedTnJc2", lambda: [make_exported_tnjc2()]),
])
def test_record_csv_save_load(record_type, make_func):
    """Test CSV save/load round-trip for each Record class."""
    # Get record class and create sample records
    record_classes = {
        "RefTnSide": RefTnSide,
        "RefTnLoc": RefTnLoc,
        "SeqRefTnSide": SeqRefTnSide,
        "BlastHit": BlastHit,
        "Junction": Junction,
        "RefTnJunction": RefTnJunction,
        "TnJunction": TnJunction,
        "RawTnJc2": RawTnJc2,
        "CoveredTnJc2": CoveredTnJc2,
        "ClassifiedTnJc2": ClassifiedTnJc2,
        "FilteredTnJc2": FilteredTnJc2,
        "AnalyzedTnJc2": AnalyzedTnJc2,
        "ExportedTnJc2": ExportedTnJc2,
    }
    
    record_cls = record_classes[record_type]
    sample_records = make_func()
    
    # Create CSV path
    csv_path = TEST_OUTPUT_DIR / f"{record_type.lower()}_test.csv"
    
    # Save to CSV
    df = RecordTypedDf.from_records(sample_records, record_cls)
    df.to_csv(csv_path)
    
    # Verify file exists
    assert csv_path.exists(), f"CSV file not created: {csv_path}"
    
    # Load from CSV
    df2 = RecordTypedDf.from_csv(csv_path, record_cls)
    
    # Verify record count
    assert len(df2) == len(sample_records), f"Record count mismatch for {record_type}"
    
    # Verify each record matches
    loaded_records = list(df2)
    for i, (original, loaded) in enumerate(zip(sample_records, loaded_records)):
        # Compare using model_dump for comprehensive comparison
        orig_dict = original.model_dump()
        loaded_dict = loaded.model_dump()
        
        # Check all fields match
        for key in orig_dict.keys():
            assert key in loaded_dict, f"Missing field {key} in loaded record {i} for {record_type}"
            assert orig_dict[key] == loaded_dict[key], (
                f"Field {key} mismatch in record {i} for {record_type}: "
                f"expected {orig_dict[key]}, got {loaded_dict[key]}"
            )


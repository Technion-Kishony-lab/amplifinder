"""Test CSV save/load for all Record classes using TypedDf."""

import pytest
from pathlib import Path

from amplifinder.data_types import (
    RecordTypedDf,
    RefTnSide, OffsetRefTnSide, RefTn, BlastHit,
    Junction, RefTnJunction, TnJunction,
    RawTnJc2, CoveredTnJc2, SingleLocusLinkedTnJc2, SynJctsTnJc2, AnalyzedTnJc2, ExportedTnJc2,
    Side, Orientation, RawEvent, EventModifier, SeqScaffold,
)


# Output directory for test CSV files (visible after test run)
TEST_OUTPUT_DIR = Path(__file__).parent.parent / "test_outputs" / "record_csv_tests"


@pytest.fixture(scope="module", autouse=True)
def setup_output_dir():
    """Create output directory for test CSV files."""
    TEST_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    yield
    # Keep files after test for inspection


def make_ref_tn_side(tn_id: int = 1, side: Side = Side.START) -> RefTnSide:
    """Create sample RefTnSide."""
    return RefTnSide(tn_id=tn_id, side=side)


def make_offset_ref_tn_side(tn_id: int = 1, side: Side = Side.START, offset: int = 0) -> OffsetRefTnSide:
    """Create sample OffsetRefTnSide."""
    return OffsetRefTnSide(tn_id=tn_id, side=side, offset=offset)


def make_ref_tn(tn_id: int = 1) -> RefTn:
    """Create sample RefTn."""
    return RefTn(
        tn_id=tn_id,
        tn_name="IS1",
        scaf="chr1",
        is_circular=False,
        length=1000,
        start=100,
        end=200,
        orientation=Orientation.FORWARD,
        join=False,
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
        scaf1="chr1",
        pos1=100,
        dir1=Orientation.FORWARD,
        scaf2="chr1",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking1=50,
        flanking2=50,
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
        flanking1=50,
        flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Side.START),
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
        flanking1=50,
        flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.START, offset=0)],
        swapped=False,
    )


def make_raw_tnjc2() -> RawTnJc2:
    """Create sample RawTnJc2."""
    from amplifinder.data_types import OffsetRefTnSide, Side, Scaffold
    
    tn_jc_S = TnJunction(
        num=1,
        scaf1="chr1",
        pos1=100,
        dir1=Orientation.FORWARD,
        scaf2="chr1",
        pos2=200,
        dir2=Orientation.FORWARD,
        flanking1=50,
        flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.START, offset=0)],
        swapped=False,
    )
    tn_jc_E = TnJunction(
        num=2,
        scaf1="chr1",
        pos1=300,
        dir1=Orientation.REVERSE,
        scaf2="chr1",
        pos2=400,
        dir2=Orientation.REVERSE,
        flanking1=50,
        flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.END, offset=0)],
        swapped=False,
    )
    scaffold = SeqScaffold(scaf="chr1", seq="A"*1000, is_circular=False)
    return RawTnJc2(
        tnjc_left=tn_jc_S,
        tnjc_right=tn_jc_E,
        scaffold=scaffold,
    )


def make_covered_tnjc2() -> CoveredTnJc2:
    """Create sample CoveredTnJc2."""
    raw = make_raw_tnjc2()
    return CoveredTnJc2.from_other(
        raw,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )


def make_classified_tnjc2() -> SingleLocusLinkedTnJc2:
    """Create sample SingleLocusLinkedTnJc2."""
    covered = make_covered_tnjc2()
    return SingleLocusLinkedTnJc2.from_other(
        covered,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1],
        chosen_tn_id=1,
    )


def make_syn_jcts_tnjc2() -> SynJctsTnJc2:
    """Create sample SynJctsTnJc2."""
    classified = make_classified_tnjc2()
    return SynJctsTnJc2.from_other(
        classified,
        analysis_dir="tn_jc2_001",
        analysis_dir_anc="tn_jc2_001_anc",
    )


def make_analyzed_tnjc2() -> AnalyzedTnJc2:
    """Create sample AnalyzedTnJc2."""
    syn_jcts = make_syn_jcts_tnjc2()
    return AnalyzedTnJc2.from_other(
        syn_jcts,
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
    ("RefTnSide", lambda: [make_ref_tn_side(), make_ref_tn_side(2, Side.END)]),
    ("OffsetRefTnSide", lambda: [make_offset_ref_tn_side(), make_offset_ref_tn_side(2, Side.END, 5),
                                 make_offset_ref_tn_side(3, Side.START, 10)]),
    ("RefTn", lambda: [make_ref_tn(), make_ref_tn(2)]),
    ("BlastHit", lambda: [make_blast_hit()]),
    ("Junction", lambda: [make_junction(1), make_junction(2)]),
    ("RefTnJunction", lambda: [make_ref_tn_junction()]),
    ("TnJunction", lambda: [make_tn_junction()]),
    # Skip RawTnJc2 and derivatives - they have CSV_EXPORT_FIELDS that only exports derived properties
    # and can't be round-tripped through CSV (complex nested objects not exported)
    # ("RawTnJc2", lambda: [make_raw_tnjc2()]),
    # ("CoveredTnJc2", lambda: [make_covered_tnjc2()]),
    # ("SingleLocusLinkedTnJc2", lambda: [make_classified_tnjc2()]),
    # ("SynJctsTnJc2", lambda: [make_syn_jcts_tnjc2()]),
    # ("AnalyzedTnJc2", lambda: [make_analyzed_tnjc2()]),
    ("ExportedTnJc2", lambda: [make_exported_tnjc2()]),
])
def test_record_csv_save_load(record_type, make_func):
    """Test CSV save/load round-trip for each Record class."""
    # Get record class and create sample records
    record_classes = {
        "RefTnSide": RefTnSide,
        "OffsetRefTnSide": OffsetRefTnSide,
        "RefTn": RefTn,
        "BlastHit": BlastHit,
        "Junction": Junction,
        "RefTnJunction": RefTnJunction,
        "TnJunction": TnJunction,
        "RawTnJc2": RawTnJc2,
        "CoveredTnJc2": CoveredTnJc2,
        "SingleLocusLinkedTnJc2": SingleLocusLinkedTnJc2,
        "SynJctsTnJc2": SynJctsTnJc2,
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

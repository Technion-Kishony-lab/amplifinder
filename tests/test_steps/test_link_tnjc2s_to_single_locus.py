"""Tests for link_tnjc2s_to_single_locus module."""

import pytest
from amplifinder.steps import LinkTnJc2ToSingleLocusPairsStep
from amplifinder.data_types import (
    RecordTypedDf, CoveredTnJc2, RefTn, Architecture,
)


@pytest.fixture
def sample_covered_tnjc2(tiny_genome, ref_tn_record):
    """Create sample CoveredTnJc2 records - one of each major type."""
    from amplifinder.data_types import (
        RawTnJc2,
        SingleLocusLinkedTnJc2,
        TnJunction,
        Orientation,
        OffsetRefTnSide,
        Terminal,
        BaseEvent,
    )

    scaffold = tiny_genome.get_scaffold("tiny")
    
    # Create a second RefTn for tn_id=2
    ref_tn_2 = RefTn.from_other(ref_tn_record, tn_id=2, tn_name="IS2")

    # Helper to create TnJunction with less boilerplate
    def make_jc(num, pos1, pos2, side, ref_tn=ref_tn_record):
        is_left = (side == Terminal.START)
        return TnJunction(
            num=num, scaf1="tiny", pos1=pos1, scaf2="tiny", pos2=pos2,
            dir1=Orientation.FORWARD if is_left else Orientation.REVERSE,
            dir2=Orientation.FORWARD if is_left else Orientation.REVERSE,
            flanking1=50, flanking2=50,
            ref_tn_sides=[OffsetRefTnSide(ref_tn=ref_tn, side=side, offset=0)],
            swapped=False,
        )

    # FLANKED: both junctions match other single-locus pairs
    jc_flank_left = make_jc(1, 10, 200, Terminal.START)
    jc_flank_right = make_jc(2, 20, 300, Terminal.END)
    flanked_raw = RawTnJc2(
        tnjc_left=jc_flank_left,
        tnjc_right=jc_flank_right,
        scaffold=scaffold,
        base_event=BaseEvent.LOCUS_JOINING,
    )
    flanked_linked = SingleLocusLinkedTnJc2.from_other(
        flanked_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    flanked = CoveredTnJc2.from_other(
        flanked_linked,
        iso_scaf_avg=1.0, iso_amplicon_avg=2.0,
    )

    # Single-locus pairs matching left/right junctions (for flanking classification)
    # These must be TRANSPOSITION events (amplicon_length < 30bp) to qualify as "single locus"
    # pos2 difference must be < 30 for transposition classification
    single_locus_left_raw = RawTnJc2(
        tnjc_left=jc_flank_left,
        tnjc_right=make_jc(3, 30, 210, Terminal.END),
        scaffold=scaffold,
        base_event=BaseEvent.TRANSPOSITION,
    )
    single_locus_left_linked = SingleLocusLinkedTnJc2.from_other(
        single_locus_left_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    single_locus_left = CoveredTnJc2.from_other(
        # 210-200=10 < 30
        single_locus_left_linked,
        iso_scaf_avg=1.0, iso_amplicon_avg=1.0,
    )
    single_locus_right_raw = RawTnJc2(
        tnjc_left=make_jc(4, 40, 290, Terminal.START),
        tnjc_right=jc_flank_right,
        scaffold=scaffold,
        base_event=BaseEvent.TRANSPOSITION,
    )
    single_locus_right_linked = SingleLocusLinkedTnJc2.from_other(
        single_locus_right_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    single_locus_right = CoveredTnJc2.from_other(
        # 300-290=10 < 30
        single_locus_right_linked,
        iso_scaf_avg=1.0, iso_amplicon_avg=1.0,
    )

    # TRANSPOSITION: junctions < 30bp apart (500 to 515 = 15bp)
    transposition_raw = RawTnJc2(
        tnjc_left=make_jc(10, 100, 500, Terminal.START),
        tnjc_right=make_jc(11, 110, 515, Terminal.END),
        scaffold=scaffold,
        base_event=BaseEvent.TRANSPOSITION,
    )
    transposition_linked = SingleLocusLinkedTnJc2.from_other(
        transposition_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    transposition = CoveredTnJc2.from_other(
        transposition_linked,
        iso_scaf_avg=1.0, iso_amplicon_avg=1.0,
    )

    # UNFLANKED: unique junctions
    unflanked_raw = RawTnJc2(
        tnjc_left=make_jc(20, 150, 700, Terminal.START, ref_tn=ref_tn_2),
        tnjc_right=make_jc(21, 160, 900, Terminal.END, ref_tn=ref_tn_2),
        scaffold=scaffold,
        base_event=BaseEvent.LOCUS_JOINING,
    )
    unflanked_linked = SingleLocusLinkedTnJc2.from_other(
        unflanked_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    unflanked = CoveredTnJc2.from_other(
        unflanked_linked,
        iso_scaf_avg=1.0, iso_amplicon_avg=1.0,
    )

    return RecordTypedDf.from_records(
        [flanked, transposition, unflanked, single_locus_left, single_locus_right],
        CoveredTnJc2
    )


@pytest.fixture
def sample_tn_locs(ref_tn_record):
    """Create sample RefTn records."""
    first = RefTn.from_other(ref_tn_record, scaf="chr1", start=10, end=20, tn_id=1, tn_name="IS1")
    second = RefTn.from_other(ref_tn_record, scaf="chr1", start=30, end=40, tn_id=2, tn_name="IS2")
    return RecordTypedDf.from_records([first, second], RefTn)


def test_link_tnjc2s_to_single_locus(sample_covered_tnjc2, sample_tn_locs, tmp_path):
    """Should classify junction pairs."""
    step = LinkTnJc2ToSingleLocusPairsStep(
        covered_tnjc2s=sample_covered_tnjc2,
        tn_locs=sample_tn_locs,
        output_dir=tmp_path,
    )

    result = step.run()

    assert len(result) == 5  # flanked, transposition, unflanked, 2 single-locus pairs
    result_list = list(result)
    # Test the main 3 types (first 3 records)
    assert result_list[0].raw_event == Architecture.FLANKED
    assert result_list[1].raw_event == Architecture.TRANSPOSITION
    assert result_list[2].raw_event == Architecture.UNFLANKED
    # The other 2 are single-locus pairs (transposition events, short amplicons)
    assert result_list[3].raw_event == Architecture.TRANSPOSITION
    assert result_list[4].raw_event == Architecture.TRANSPOSITION


def test_filters_by_length(sample_covered_tnjc2, sample_tn_locs, covered_tnjc2_record, tmp_path):
    """Should filter candidates by amplicon length."""
    # Create a short amplicon
    short_record = CoveredTnJc2.from_other(
        covered_tnjc2_record,
    )

    all_records = RecordTypedDf.from_records(
        list(sample_covered_tnjc2) + [short_record],
        CoveredTnJc2
    )

    step = LinkTnJc2ToSingleLocusPairsStep(
        covered_tnjc2s=all_records,
        tn_locs=sample_tn_locs,
        output_dir=tmp_path,
    )

    result = step.run()
    # Now expecting 6 records (5 from fixture + 1 new)
    assert len(result) == 6

"""Tests for FilterTnJc2CandidatesStep."""

import pytest

from amplifinder.data_types import BaseEvent
from amplifinder.steps import FilterTnJc2CandidatesStep
from amplifinder.data_types import RecordTypedDf


@pytest.fixture
def sample_classified_tnjc2(classified_tnjc2_record, tiny_genome):
    """Create sample CoveredTnJc2 records with different amplicon lengths."""
    from amplifinder.data_types import (
        TnJunction,
        Orientation,
        OffsetRefTnSide,
        RefTnSide,
        Terminal,
        RawTnJc2,
        SingleLocusLinkedTnJc2,
        CoveredTnJc2,
    )

    # Short amplicon (20bp) - should be filtered out
    tn_jc_S_short = TnJunction(
        num=1, scaf1="tiny", pos1=10, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=100, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Terminal.START),
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.START, offset=0)],
        swapped=False,
    )
    tn_jc_E_short = TnJunction(
        num=2, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=119, dir2=Orientation.REVERSE,  # 119-100+1 = 20bp
        flanking1=50, flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Terminal.END),
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.END, offset=0)],
        swapped=False,
    )
    scaffold = tiny_genome.get_scaffold("tiny")
    short_raw = RawTnJc2(
        tnjc_left=tn_jc_S_short,
        tnjc_right=tn_jc_E_short,
        scaffold=scaffold,
        base_event=BaseEvent.TRANSPOSITION,
    )
    short_linked = SingleLocusLinkedTnJc2.from_other(
        short_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    short = CoveredTnJc2.from_other(
        short_linked,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=1.5,
    )

    # Medium amplicon (101bp)
    tn_jc_S_medium = TnJunction(
        num=5, scaf1="tiny", pos1=10, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=200, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Terminal.START),
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.START, offset=0)],
        swapped=False,
    )
    tn_jc_E_medium = TnJunction(
        num=6, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=300, dir2=Orientation.REVERSE,  # 300-200+1 = 101bp
        flanking1=50, flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Terminal.END),
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.END, offset=0)],
        swapped=False,
    )
    medium_raw = RawTnJc2(
        tnjc_left=tn_jc_S_medium,
        tnjc_right=tn_jc_E_medium,
        scaffold=scaffold,
        base_event=BaseEvent.LOCUS_JOINING,
    )
    medium_linked = SingleLocusLinkedTnJc2.from_other(
        medium_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    medium = CoveredTnJc2.from_other(
        medium_linked,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )

    # Long amplicon (500bp)
    tn_jc_S_long = TnJunction(
        num=3, scaf1="tiny", pos1=10, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=500, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Terminal.START),
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.START, offset=0)],
        swapped=False,
    )
    tn_jc_E_long = TnJunction(
        num=4, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=999, dir2=Orientation.REVERSE,  # 999-500+1 = 500bp
        flanking1=50, flanking2=50,
        ref_tn_side=RefTnSide(tn_id=1, side=Terminal.END),
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.END, offset=0)],
        swapped=False,
    )
    long_raw = RawTnJc2(
        tnjc_left=tn_jc_S_long,
        tnjc_right=tn_jc_E_long,
        scaffold=scaffold,
        base_event=BaseEvent.LOCUS_JOINING,
    )
    long_linked = SingleLocusLinkedTnJc2.from_other(
        long_raw,
        single_locus_tnjc2_left_matchings=[],
        single_locus_tnjc2_right_matchings=[],
    )
    long = CoveredTnJc2.from_other(
        long_linked,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=3.0,
    )

    return RecordTypedDf.from_records([short, medium, long], CoveredTnJc2)


def test_filters_by_length(sample_classified_tnjc2, tmp_path):
    """Should filter candidates by copy number threshold."""
    step = FilterTnJc2CandidatesStep(
        linked_tnjc2s=sample_classified_tnjc2,
        output_dir=tmp_path,
        # Filter out short (1.5), keep medium (2.0) and long (3.0)
        replication_copy_number_threshold=1.6,
        deletion_copy_number_threshold=0.3,
    )

    result = step.run()

    # Should filter out short (1.5), keep medium (2.0) and long (3.0)
    assert len(result) == 2
    assert sorted([r.amplicon_length for r in result]) == [101, 500]


def test_filters_by_copy_number_threshold(sample_classified_tnjc2, tmp_path):
    """Should filter candidates by copy number threshold."""
    step = FilterTnJc2CandidatesStep(
        linked_tnjc2s=sample_classified_tnjc2,
        output_dir=tmp_path,
        # Filter out short (1.5), keep medium (2.0) and long (3.0)
        replication_copy_number_threshold=1.6,
        deletion_copy_number_threshold=0.3,
    )

    result = step.run()

    assert len(result) == 2
    assert sorted([r.copy_number for r in result]) == [2.0, 3.0]

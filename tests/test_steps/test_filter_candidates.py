"""Tests for FilterTnJc2CandidatesStep."""

from amplifinder.data_types.enums import BaseRawEvent
import pytest
from amplifinder.steps import FilterTnJc2CandidatesStep
from amplifinder.data_types import (
    RecordTypedDf, ClassifiedTnJc2, RawEvent,
)


@pytest.fixture
def sample_classified_tnjc2(classified_tnjc2_record, tiny_genome):
    """Create sample ClassifiedTnJc2 records with different amplicon lengths."""
    from amplifinder.data_types import TnJunction, Orientation, OffsetRefTnSide, Side, RawTnJc2
    
    # Short amplicon (20bp) - should be filtered out
    tn_jc_S_short = TnJunction(
        num=1, scaf1="tiny", pos1=10, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=100, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.START, offset=0)],
        swapped=False,
    )
    tn_jc_E_short = TnJunction(
        num=2, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=119, dir2=Orientation.REVERSE,  # 119-100+1 = 20bp
        flanking1=50, flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.END, offset=0)],
        swapped=False,
    )
    scaffold = tiny_genome.get_scaffold("tiny")
    short_raw = RawTnJc2(tnjc_left=tn_jc_S_short, tnjc_right=tn_jc_E_short, scaffold=scaffold)
    short = ClassifiedTnJc2.from_other(
        short_raw,
        base_raw_event=BaseRawEvent.TRANSPOSITION,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=1.5,
    )

    # Medium amplicon (101bp) - from base fixture
    medium = ClassifiedTnJc2.from_other(
        classified_tnjc2_record,
        base_raw_event=BaseRawEvent.LOCUS_JOINING,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=2.0,
    )

    # Long amplicon (500bp)
    tn_jc_S_long = TnJunction(
        num=3, scaf1="tiny", pos1=10, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=500, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.START, offset=0)],
        swapped=False,
    )
    tn_jc_E_long = TnJunction(
        num=4, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=999, dir2=Orientation.REVERSE,  # 999-500+1 = 500bp
        flanking1=50, flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Side.END, offset=0)],
        swapped=False,
    )
    long_raw = RawTnJc2(tnjc_left=tn_jc_S_long, tnjc_right=tn_jc_E_long, scaffold=scaffold)
    long = ClassifiedTnJc2.from_other(
        long_raw,
        base_raw_event=BaseRawEvent.LOCUS_JOINING,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=3.0,
    )

    return RecordTypedDf.from_records([short, medium, long], ClassifiedTnJc2)


def test_filters_by_length(sample_classified_tnjc2, tmp_path):
    """Should filter candidates by amplicon length."""
    step = FilterTnJc2CandidatesStep(
        classified_tnjc2s=sample_classified_tnjc2,
        output_dir=tmp_path,
        min_amplicon_length=50,  # Filter out the 20bp amplicon
        max_amplicon_length=1_000_000,
    )

    result = step.run()

    # Should filter out short (20bp), keep medium (101bp) and long (500bp)
    assert len(result) == 2
    result_list = list(result)
    assert result_list[0].amplicon_length == 101
    assert result_list[1].amplicon_length == 500


def test_keeps_all_when_no_length_filter(sample_classified_tnjc2, tmp_path):
    """Should keep all candidates when length filter is permissive."""
    step = FilterTnJc2CandidatesStep(
        classified_tnjc2s=sample_classified_tnjc2,
        output_dir=tmp_path,
        min_amplicon_length=10,  # Keep all
        max_amplicon_length=1_000_000,
    )

    result = step.run()

    assert len(result) == 3  # All pass length filter
    result_list = list(result)
    
    # Check all have expected amplicon lengths
    lengths = sorted([r.amplicon_length for r in result_list])
    assert lengths == [20, 101, 500]

"""Tests for CreateSyntheticJunctionsStep."""

from amplifinder.steps import CreateSyntheticJunctionsStep
from amplifinder.data_types import (
    RecordTypedDf, FilteredTnJc2, Orientation,
)
from amplifinder.utils.file_utils import ensure_parent_dir


def test_creates_junctions_fasta(tiny_genome, filtered_tnjc2_record, ref_tns_indexed, tmp_path):
    """Should create junctions.fasta file."""
    filtered_tnjc2s = RecordTypedDf.from_records([filtered_tnjc2_record], FilteredTnJc2)

    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        ref_tns=ref_tns_indexed,
        output_dir=tmp_path,
        read_length=150,
    )

    step.run()

    # Check that junctions.fasta was created
    junctions_file = tmp_path / "junctions" / filtered_tnjc2_record.analysis_dir / "junctions.fasta"
    assert junctions_file.exists()

    # Check that it contains 7 junction sequences
    with open(junctions_file) as f:
        content = f.read()
        # Should have 7 headers (one per junction type)
        assert content.count(">") == 7


def test_handles_missing_tn(tiny_genome, ref_tns_indexed, tmp_path):
    """Should skip candidates with missing TN (chosen_tn_id=None)."""
    from amplifinder.data_types import FilteredTnJc2, BaseRawEvent, TnJunction, RawTnJc2, Orientation, Side
    
    # Create TnJunctions with no ref_tn_sides (empty list) so tn_ids will be empty
    tn_jc_left = TnJunction(
        num=100, scaf1="tiny", pos1=200, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=300, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_sides=[],  # No TN sides
        swapped=False,
    )
    tn_jc_right = TnJunction(
        num=101, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=500, dir2=Orientation.REVERSE,
        flanking1=50, flanking2=50,
        ref_tn_sides=[],  # No TN sides
        swapped=False,
    )
    scaffold = tiny_genome.get_scaffold("tiny")
    raw_tnjc2 = RawTnJc2(tnjc_left=tn_jc_left, tnjc_right=tn_jc_right, scaffold=scaffold)
    
    candidate_no_tn = FilteredTnJc2.from_other(
        raw_tnjc2,
        iso_scaf_avg=1.0,
        iso_amplicon_avg=1.0,
        base_raw_event=BaseRawEvent.LOCUS_JOINING,
        analysis_dir="jc_200_300_no_tn",
    )

    filtered_tnjc2s = RecordTypedDf.from_records([candidate_no_tn], FilteredTnJc2)

    # Create empty junctions.fasta to satisfy step output requirements
    junctions_file = tmp_path / "junctions" / candidate_no_tn.analysis_dir / "junctions.fasta"
    ensure_parent_dir(junctions_file)
    junctions_file.write_text("")

    step = CreateSyntheticJunctionsStep(
        filtered_tnjc2s=filtered_tnjc2s,
        genome=tiny_genome,
        ref_tns=ref_tns_indexed,
        output_dir=tmp_path,
    )

    # Should not raise error, just skip creating content
    result = step.run()
    # Step returns None (it creates files but doesn't return data)

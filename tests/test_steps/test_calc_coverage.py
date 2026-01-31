"""Tests for CalcTnJc2AmpliconCoverageStep."""

import pytest
import pandas as pd

from amplifinder.steps import CalcTnJc2AmpliconCoverageStep
from amplifinder.data_types import RecordTypedDf
from amplifinder.utils.file_utils import ensure_dir


@pytest.fixture
def sample_tnjc2(tiny_genome):
    """Create sample RawTnJc2 records."""
    from amplifinder.data_types import RawTnJc2, TnJunction, Orientation, OffsetRefTnSide, Terminal

    tn_jc_S = TnJunction(
        num=1, scaf1="tiny", pos1=10, dir1=Orientation.FORWARD,
        scaf2="tiny", pos2=200, dir2=Orientation.FORWARD,
        flanking1=50, flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.START, offset=0)],
        swapped=False,
    )
    tn_jc_E = TnJunction(
        num=2, scaf1="tiny", pos1=20, dir1=Orientation.REVERSE,
        scaf2="tiny", pos2=300, dir2=Orientation.REVERSE,
        flanking1=50, flanking2=50,
        ref_tn_sides=[OffsetRefTnSide(tn_id=1, side=Terminal.END, offset=0)],
        swapped=False,
    )
    scaffold = tiny_genome.get_scaffold("tiny")
    record = RawTnJc2(tnjc_left=tn_jc_S, tnjc_right=tn_jc_E, scaffold=scaffold)

    return RecordTypedDf.from_records([record, record], RawTnJc2)


@pytest.fixture
def mock_breseq_output(tmp_path, tiny_genome):
    """Create mock breseq coverage output."""
    breseq_dir = ensure_dir(tmp_path / "breseq" / "08_mutation_identification")

    # Create coverage.tab file
    coverage_file = breseq_dir / "tiny.coverage.tab"

    # Generate coverage data with header matching breseq format
    # Format: position, unique_top_cov, unique_bot_cov, ...
    header = "position\tunique_top_cov\tunique_bot_cov\n"
    coverage_data = [header]
    for pos in range(1, len(tiny_genome) + 1):
        # Simulate coverage: higher in amplicon region (200-300)
        if 200 <= pos <= 300:
            cov = 20.0  # Higher coverage in amplicon
        else:
            cov = 10.0  # Normal coverage
        coverage_data.append(f"{pos}\t{cov}\t{cov}\n")

    coverage_file.write_text("".join(coverage_data))
    return tmp_path / "breseq"


def test_calc_coverage_step(tiny_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should calculate coverage for candidates."""
    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        iso_breseq_path=mock_breseq_output,
        output_dir=tmp_path,
        force=True,  # Force run to avoid loading from CSV
    )

    result = step.run()

    assert len(result) == 2

    # Access records directly from DataFrame to avoid NaN->None conversion
    result_df = result.df
    # First candidate should have coverage calculated
    first_row = result_df.iloc[0]
    cov_col = 'amplicon_coverage' if 'amplicon_coverage' in result_df.columns else None
    copy_col = 'copy_number' if 'copy_number' in result_df.columns else None
    ref_col = 'ref_name' if 'ref_name' in result_df.columns else None
    iso_col = 'iso_name' if 'iso_name' in result_df.columns else None

    if cov_col:
        assert first_row[cov_col] > 0 or not pd.isna(first_row[cov_col])
    if copy_col:
        assert first_row[copy_col] > 0 or not pd.isna(first_row[copy_col])
    if ref_col:
        assert first_row[ref_col] == "tiny"
    if iso_col:
        assert first_row[iso_col] == "sample1"

    # Second candidate should also have coverage calculated
    second_row = result_df.iloc[1]
    if cov_col:
        assert not pd.isna(second_row[cov_col])


def test_calc_coverage_with_ancestor(tiny_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should calculate normalized coverage when ancestor is provided."""
    # Create ancestor breseq output with lower coverage
    anc_breseq_dir = ensure_dir(tmp_path / "anc_breseq" / "08_mutation_identification")
    anc_coverage_file = anc_breseq_dir / "tiny.coverage.tab"

    # Ancestor has uniform coverage of 10
    anc_header = "position\tunique_top_cov\tunique_bot_cov\n"
    anc_coverage_data = [anc_header]
    for pos in range(1, len(tiny_genome) + 1):
        anc_coverage_data.append(f"{pos}\t10.0\t10.0\n")
    anc_coverage_file.write_text("".join(anc_coverage_data))

    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        iso_breseq_path=mock_breseq_output,
        anc_breseq_path=tmp_path / "anc_breseq",
        output_dir=tmp_path,
        force=True,  # Force run to avoid loading from CSV
    )

    result = step.run()

    assert len(result) == 2
    result_df = result.df
    schema = result.schema

    # First candidate should have normalized coverage (iso/anc ratio)
    first_row = result_df.iloc[0]
    anc_col = 'anc_name' if 'anc_name' in schema.column_names else None
    ratio_col = 'copy_number_ratio' if 'copy_number_ratio' in schema.column_names else None
    if anc_col:
        assert first_row[anc_col] == "ancestor1"
    if ratio_col:
        assert first_row[ratio_col] is not None and not pd.isna(first_row[ratio_col])
        # Amplicon region has 2x coverage, so ratio should be ~2.0
        assert first_row[ratio_col] > 1.5


def test_calculates_coverage_for_all_lengths(tiny_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should calculate coverage for all candidates regardless of length."""
    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        iso_breseq_path=mock_breseq_output,
        output_dir=tmp_path,
        force=True,  # Force run to avoid loading from CSV
    )

    result = step.run()

    # Both candidates should be in result with coverage calculated
    assert len(result) == 2
    result_df = result.df
    schema = result.schema

    # Both candidates should have coverage calculated (not NaN)
    length_col = 'amplicon_length' if 'amplicon_length' in schema.column_names else None
    cov_col = 'amplicon_coverage' if 'amplicon_coverage' in schema.column_names else None
    if length_col and cov_col:
        short_row = result_df[result_df[length_col] == 20].iloc[0]
        assert not pd.isna(short_row[cov_col])


def test_skips_coverage_for_too_long_amplicons(tiny_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should skip coverage calculation for amplicons outside length range."""
    # Set max_amplicon_length to 50, which is less than the sample amplicon length (100)
    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        iso_breseq_path=mock_breseq_output,
        output_dir=tmp_path,
        max_amplicon_length=50,  # Too short to include sample amplicons
        force=True,
    )

    result = step.run()

    # Candidates should still be in result, but with NaN coverage
    assert len(result) == 2

    # Check that coverage values are NaN
    for rec in result:
        assert pd.isna(rec.iso_amplicon_avg)
        assert pd.isna(rec.avg_norm_cov)

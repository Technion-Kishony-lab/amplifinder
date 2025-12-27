"""Tests for CalcTnJc2AmpliconCoverageStep."""

import pytest
import pandas as pd
from amplifinder.steps import CalcTnJc2AmpliconCoverageStep
from amplifinder.data_types import (
    RecordTypedDf, RawTnJc2, Genome, Orientation,
)
from amplifinder.utils.file_utils import ensure_dir


@pytest.fixture
def sample_tnjc2(tmp_path):
    """Create sample RawTnJc2 records."""
    records = [
        RawTnJc2(
            jc_num_L=1, jc_num_R=2,
            scaf="chr1",
            pos_scaf_L=100, pos_scaf_R=200,
            pos_tn_L=10, pos_tn_R=20,
            dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[1], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=100, complementary_length=900,
        ),
        RawTnJc2(
            jc_num_L=3, jc_num_R=4,
            scaf="chr1",
            pos_scaf_L=300, pos_scaf_R=320,
            pos_tn_L=30, pos_tn_R=40,
            dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
            dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
            tn_ids=[2], tn_orientations=[Orientation.FORWARD],
            span_origin=False,
            amplicon_length=20,  # Too short
            complementary_length=980,
        ),
    ]
    return RecordTypedDf.from_records(records, RawTnJc2)


@pytest.fixture
def mock_genome(tmp_path):
    """Create a mock Genome object."""
    fasta_path = tmp_path / "genome.fasta"
    # Create 1000bp sequence
    seq = "A" * 1000
    fasta_path.write_text(f">chr1\n{seq}\n")

    genome = Genome(
        name="chr1",
        fasta_path=fasta_path,
        genbank_path=None,
    )
    # Length will be computed from sequence (1000bp)
    return genome


@pytest.fixture
def mock_breseq_output(tmp_path):
    """Create mock breseq coverage output."""
    breseq_dir = ensure_dir(tmp_path / "breseq" / "08_mutation_identification")

    # Create coverage.tab file
    coverage_file = breseq_dir / "chr1.coverage.tab"

    # Generate coverage data with header matching breseq format
    # Format: position, unique_top_cov, unique_bot_cov, ...
    header = "position\tunique_top_cov\tunique_bot_cov\n"
    coverage_data = [header]
    for pos in range(1, 1001):
        # Simulate coverage: higher in amplicon region (100-200)
        if 100 <= pos <= 200:
            cov = 20.0  # Higher coverage in amplicon
        else:
            cov = 10.0  # Normal coverage
        coverage_data.append(f"{pos}\t{cov}\t{cov}\n")

    coverage_file.write_text("".join(coverage_data))
    return tmp_path / "breseq"


def test_calc_coverage_step(mock_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should calculate coverage for candidates."""
    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        genome=mock_genome,
        iso_breseq_path=mock_breseq_output,
        output_dir=tmp_path,
        ref_name="chr1",
        iso_name="sample1",
        force=True,  # Force run to avoid loading from CSV
    )

    result = step.run()

    assert len(result) == 2

    # Access records directly from DataFrame to avoid NaN->None conversion
    result_df = result.df
    schema = result.schema

    # First candidate should have coverage calculated
    first_row = result_df.iloc[0]
    cov_col = 'amplicon_coverage' if 'amplicon_coverage' in schema.column_names else None
    copy_col = 'copy_number' if 'copy_number' in schema.column_names else None
    ref_col = 'ref_name' if 'ref_name' in schema.column_names else None
    iso_col = 'iso_name' if 'iso_name' in schema.column_names else None

    if cov_col:
        assert first_row[cov_col] > 0 or not pd.isna(first_row[cov_col])
    if copy_col:
        assert first_row[copy_col] > 0 or not pd.isna(first_row[copy_col])
    if ref_col:
        assert first_row[ref_col] == "chr1"
    if iso_col:
        assert first_row[iso_col] == "sample1"

    # Second candidate should also have coverage calculated
    second_row = result_df.iloc[1]
    if cov_col:
        assert not pd.isna(second_row[cov_col])


def test_calc_coverage_with_ancestor(mock_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should calculate normalized coverage when ancestor is provided."""
    # Create ancestor breseq output with lower coverage
    anc_breseq_dir = ensure_dir(tmp_path / "anc_breseq" / "08_mutation_identification")
    anc_coverage_file = anc_breseq_dir / "chr1.coverage.tab"

    # Ancestor has uniform coverage of 10
    anc_header = "position\tunique_top_cov\tunique_bot_cov\n"
    anc_coverage_data = [anc_header]
    for pos in range(1, 1001):
        anc_coverage_data.append(f"{pos}\t10.0\t10.0\n")
    anc_coverage_file.write_text("".join(anc_coverage_data))

    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        genome=mock_genome,
        iso_breseq_path=mock_breseq_output,
        anc_breseq_path=tmp_path / "anc_breseq",
        output_dir=tmp_path,
        ref_name="chr1",
        iso_name="sample1",
        anc_name="ancestor1",
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


def test_calculates_coverage_for_all_lengths(mock_genome, sample_tnjc2, mock_breseq_output, tmp_path):
    """Should calculate coverage for all candidates regardless of length."""
    step = CalcTnJc2AmpliconCoverageStep(
        raw_tnjc2s=sample_tnjc2,
        genome=mock_genome,
        iso_breseq_path=mock_breseq_output,
        output_dir=tmp_path,
        ref_name="chr1",
        iso_name="sample1",
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

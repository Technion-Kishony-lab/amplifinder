"""Tests for ISEScan tools."""

import os
import pytest

# Add conda to PATH before imports
conda_path = "/opt/anaconda3/bin"
if conda_path not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"{os.environ.get('PATH', '')}:{conda_path}"

from amplifinder.tools.isescan import run_isescan, parse_isescan_results, get_isescan_results_file
from tests.env import RUN_ISESCAN_TESTS
from amplifinder.steps.locate_tns.locate_tns import SEQ_COL, START_COL, END_COL, STRAND_COL, FAMILY_COL, CLUSTER_COL

skip_no_isescan = pytest.mark.skipif(not RUN_ISESCAN_TESTS, reason="ISEScan tests disabled")


@pytest.fixture
def isescan_output(tiny_ref_fasta, tmp_path):
    """Run ISEScan and return output directory."""
    from amplifinder.env import ISESCAN_ENV_NAME
    
    output_dir = tmp_path / "isescan_output"
    run_isescan(
        seqfile=tiny_ref_fasta,
        output_dir=output_dir,
        env_name=ISESCAN_ENV_NAME,
        threads=1,
    )
    return output_dir


@skip_no_isescan
def test_conda_available():
    """Test that conda is available in PATH."""
    import shutil
    import subprocess
    conda = shutil.which("conda")
    assert conda is not None, f"conda not found in PATH: {os.environ.get('PATH', '')}"
    # Test conda run works
    result = subprocess.run(["conda", "--version"], capture_output=True, text=True)
    assert result.returncode == 0, f"conda command failed: {result.stderr}"


@skip_no_isescan
def test_run_isescan(isescan_output, tiny_ref_fasta):
    """Should run ISEScan and create output directory."""
    assert isescan_output.exists()
    # ISEScan may not produce output files if no IS elements found (expected for synthetic data)
    # Just verify it ran without errors (fixture would have failed if ISEScan crashed)


@skip_no_isescan
def test_parse_isescan_results(tiny_ref_fasta, isescan_output):
    """Should parse ISEScan results into DataFrame if results exist."""
    results_file = get_isescan_results_file(tiny_ref_fasta, isescan_output)
    assert results_file is not None

    df = parse_isescan_results(tiny_ref_fasta, isescan_output)
    
    # Should return a DataFrame
    assert df is not None
    assert hasattr(df, 'columns')
    
    # Check basic structure if data exists
    if len(df) > 0:
        # ISEScan output typically has these columns
        expected_cols = [SEQ_COL, FAMILY_COL, CLUSTER_COL, START_COL, END_COL, STRAND_COL]
        for col in expected_cols:
            assert col in df.columns

    # Check the strand column is + or -
    for _, row in df.iterrows():
        assert row[STRAND_COL] in {"+", "-"}


@skip_no_isescan
def test_parse_isescan_no_results_raises(tmp_path, tiny_ref_fasta):
    """Should raise FileNotFoundError when no results exist."""
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    
    with pytest.raises(FileNotFoundError):
        parse_isescan_results(tiny_ref_fasta, empty_dir)


@pytest.fixture
def ecoli_with_is_fasta(tmp_path):
    """Download or use cached E. coli genome with known IS elements."""
    # Use U00096 (E. coli K-12 MG1655) which has known IS elements
    from amplifinder.steps import GetRefGenomeStep
    from tests.conftest import TEST_GENOMES_CACHE
    from amplifinder.utils.file_utils import ensure_dir
    
    cache_dir = ensure_dir(TEST_GENOMES_CACHE)
    
    step = GetRefGenomeStep(
        ref_name="U00096",
        ref_path=cache_dir,
        ncbi=True,
    )
    genome = step.run()
    return genome.fasta_path


@pytest.fixture
def ecoli_isescan_output(ecoli_with_is_fasta, tmp_path):
    """Run ISEScan on real E. coli genome."""
    from amplifinder.env import ISESCAN_ENV_NAME
    
    output_dir = tmp_path / "ecoli_isescan_output"
    print(f"Running ISEScan on {ecoli_with_is_fasta} ...", flush=True)
    run_isescan(
        seqfile=ecoli_with_is_fasta,
        output_dir=output_dir,
        env_name=ISESCAN_ENV_NAME,
        threads=2,
    )
    print(f"ISEScan output directory: {output_dir}", flush=True)
    return output_dir


@skip_no_isescan
@pytest.mark.slow
@pytest.mark.integration
def test_isescan_on_real_genome(ecoli_isescan_output, ecoli_with_is_fasta):
    """Run ISEScan on real E. coli genome and parse actual output."""
    assert ecoli_isescan_output.exists()
    
    # Should find results file
    results_file = get_isescan_results_file(ecoli_with_is_fasta, ecoli_isescan_output)
    assert results_file is not None, "ISEScan should produce results for E. coli K-12"
    assert results_file.exists()
    
    # Parse the actual ISEScan output
    df = parse_isescan_results(ecoli_with_is_fasta, ecoli_isescan_output)
    
    # E. coli K-12 MG1655 has known IS elements
    print(f"\n=== ISEScan Results for E. coli K-12 ===")
    print(f"Total IS elements found: {len(df)}")
    print(f"Columns: {list(df.columns)}")
    
    if len(df) > 0:
        print(f"\nFirst few results:")
        print(df.head())
        
        # Verify expected columns
        assert 'seqID' in df.columns
        assert 'family' in df.columns
        assert 'isBegin' in df.columns
        assert 'isEnd' in df.columns
        
        # E. coli K-12 should have IS elements (known to have IS1, IS5, etc.)
        assert len(df) > 0, "E. coli K-12 should have detectable IS elements"
        
        # Check that we got reasonable data
        for idx, row in df.head(3).iterrows():
            print(f"\nIS element {idx}:")
            print(f"  Family: {row.get('family', 'N/A')}")
            print(f"  Position: {row.get('isBegin', 'N/A')} - {row.get('isEnd', 'N/A')}")
            print(f"  Length: {row.get('isLen', 'N/A')}")
    else:
        pytest.fail("ISEScan found no IS elements in E. coli K-12, but it should have several")

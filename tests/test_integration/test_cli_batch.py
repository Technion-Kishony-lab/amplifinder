"""Integration tests for CLI batch processing."""

import pytest
import pandas as pd
from pathlib import Path
from amplifinder.utils.file_utils import ensure_dir


pytestmark = pytest.mark.integration


@pytest.fixture
def batch_csv_path(isolate_srr25242877, cleared_output_dir, tmp_path):
    """Create a batch CSV with two copies of the same isolate."""
    csv_path = tmp_path / "test_batch.csv"
    
    # Create batch CSV with same isolate duplicated under different names
    # Note: CSV column names must match Config field names exactly
    batch_data = [
        {
            "iso_fastq_path": str(isolate_srr25242877["fastq_path"]),
            "ref_name": "U00096",
            "iso_name": "sample1",
            "ref_path": str(cleared_output_dir.parent / "genomesDB"),
            "output_dir": str(cleared_output_dir / "sample1"),
            "iso_breseq_path": str(isolate_srr25242877["breseq_path"]),
            "ncbi": "true",
            "use_isfinder": "false",
            "create_plots": "false",
        },
        {
            "iso_fastq_path": str(isolate_srr25242877["fastq_path"]),
            "ref_name": "U00096",
            "iso_name": "sample2",
            "ref_path": str(cleared_output_dir.parent / "genomesDB"),
            "output_dir": str(cleared_output_dir / "sample2"),
            "iso_breseq_path": str(isolate_srr25242877["breseq_path"]),
            "ncbi": "true",
            "use_isfinder": "false",
            "create_plots": "false",
        },
    ]
    
    df = pd.DataFrame(batch_data)
    df.to_csv(csv_path, index=False)
    
    return csv_path


@pytest.fixture
def status_csv_path(tmp_path):
    """Path for batch output status CSV."""
    return tmp_path / "batch_status.csv"


@pytest.mark.slow
def test_cli_batch_processing(batch_csv_path, status_csv_path, cleared_output_dir):
    """Test CLI batch processing with two samples."""
    import subprocess
    import sys
    
    # Run CLI with batch CSV
    cmd = [
        sys.executable,
        "-m",
        "amplifinder",
        "--batch-input",
        str(batch_csv_path),
        "--batch-output",
        str(status_csv_path),
        "--max-parallel",
        "1",  # Run sequentially for clearer logs
    ]
    
    result = subprocess.run(cmd, capture_output=False, text=True)
    
    # Check exit code
    assert result.returncode == 0, f"CLI failed with exit code {result.returncode}"
    
    # Verify status CSV was created
    assert status_csv_path.exists(), "Status CSV not created"
    
    # Load and check status CSV
    status_df = pd.read_csv(status_csv_path)
    assert len(status_df) == 2, f"Expected 2 runs in status CSV, got {len(status_df)}"
    
    # Check both runs succeeded
    assert all(status_df["exit_code"] == 0), f"Some runs failed:\n{status_df}"
    
    # Check output directories were created
    sample1_dir = cleared_output_dir / "sample1"
    sample2_dir = cleared_output_dir / "sample2"
    assert sample1_dir.exists(), "sample1 output directory not created"
    assert sample2_dir.exists(), "sample2 output directory not created"
    
    # Check that each sample has expected output files
    for sample_name, sample_dir in [("sample1", sample1_dir), ("sample2", sample2_dir)]:
        # Directory structure: output_dir/ref_name/iso_name/iso_name/
        iso_run_dir = sample_dir / "U00096" / sample_name / sample_name
        assert iso_run_dir.exists(), f"{sample_name} run directory not created"
        
        # Check for key output files (CSVs that should exist)
        tnjc_csv = iso_run_dir / "tnjc.csv"
        tnjc2_raw_csv = iso_run_dir / "tnjc2_A_raw.csv"
        assert tnjc_csv.exists(), f"{sample_name} tnjc.csv not found"
        assert tnjc2_raw_csv.exists(), f"{sample_name} tnjc2_A_raw.csv not found"
        
        # Verify CSVs are readable and have content
        tnjc = pd.read_csv(tnjc_csv)
        assert len(tnjc) > 0, f"{sample_name} tnjc.csv is empty"
    
    print(f"\n✓ Successfully processed 2 samples in batch mode")
    print(f"  Sample1 results: {sample1_dir / 'U00096' / 'sample1' / 'sample1'}")
    print(f"  Sample2 results: {sample2_dir / 'U00096' / 'sample2' / 'sample2'}")

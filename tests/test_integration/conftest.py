"""Fixtures for integration tests using real AmpliFinder test data."""

import pytest
from pathlib import Path

# Paths to test data
TEST_DATA_ROOT = Path("/zdata/user-data/rkishony/amplifinder/AmpliFinder_test")
MATLAB_OUTPUT = TEST_DATA_ROOT / "AmpliFinderWorkspace" / "output"

# External data paths (from isolates.xlsx)
FASTQ_PATH = Path("/zdata/user-data/idan/small_projects/breseq_to_amplification_from_dropbox/morbidostat_sra")
BRESEQ_PATH = Path("/zdata/user-data/idan/small_projects/breseq_to_amplification_from_dropbox/Breseq1")


@pytest.fixture
def test_data_root():
    """Root directory for AmpliFinder test data."""
    return TEST_DATA_ROOT


@pytest.fixture
def matlab_output():
    """Directory containing MATLAB output for comparison."""
    return MATLAB_OUTPUT


@pytest.fixture
def isolate_srr25242877():
    """Test isolate SRR25242877 configuration."""
    return {
        "iso_name": "SRR25242877",
        "ref_name": "U00096",
        "fastq_path": FASTQ_PATH / "SRR25242877",
        "breseq_path": BRESEQ_PATH / "SRR25242877",
        "matlab_output": MATLAB_OUTPUT / "SRR25242877",
    }


@pytest.fixture
def isolate_srr25242906():
    """Test isolate SRR25242906 (ancestor)."""
    return {
        "iso_name": "SRR25242906",
        "ref_name": "U00096",
        "fastq_path": FASTQ_PATH / "SRR25242906",
        "breseq_path": BRESEQ_PATH / "SRR25242906",
        "matlab_output": MATLAB_OUTPUT / "SRR25242906",
    }


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as integration test (requires real data)"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow (e.g., runs breseq)"
    )

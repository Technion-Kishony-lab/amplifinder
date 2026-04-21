"""Tests for ISEScan tools."""

import os
import pytest
import pandas as pd

# Add conda to PATH before imports
conda_path = "/opt/anaconda3/bin"
if conda_path not in os.environ.get("PATH", ""):
    os.environ["PATH"] = f"{os.environ.get('PATH', '')}:{conda_path}"

from amplifinder.tools.isescan import run_isescan, parse_isescan_results, get_isescan_results_file  # noqa: E402
from amplifinder.steps.locate_tns.locate_tns import (  # noqa: E402
    LocateTNsUsingISEScanStep, SEQ_COL, START_COL, END_COL, STRAND_COL, FAMILY_COL, CLUSTER_COL,
)
from tests.env import RUN_ISESCAN_TESTS  # noqa: E402

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
    """Should run ISEScan and create output directory with results."""
    assert isescan_output.exists()
    results_file = get_isescan_results_file(tiny_ref_fasta, isescan_output)
    assert results_file is not None, "ISEScan should produce a results file for a genome with real IS sequences"
    assert results_file.exists()


@skip_no_isescan
def test_parse_isescan_results(tiny_ref_fasta, isescan_output):
    """Should detect both IS elements embedded in the tiny reference genome.

    tiny_ref layout (1-based):
        500 bp flank | IS1 (902 bp) @ 501-1402 | 400 bp flank | IS5-RC (1195 bp) @ 1803-2997 | 500 bp flank
    """
    df = parse_isescan_results(tiny_ref_fasta, isescan_output)

    expected_cols = [SEQ_COL, FAMILY_COL, CLUSTER_COL, START_COL, END_COL, STRAND_COL]
    for col in expected_cols:
        assert col in df.columns

    assert len(df) == 2

    df = df.sort_values(START_COL).reset_index(drop=True)
    TOL = 100  # ISEScan HMM boundaries may differ slightly from embedded coords

    is1 = df.iloc[0]
    assert is1[SEQ_COL] == "tiny"
    assert is1[FAMILY_COL] == "IS1"
    assert abs(is1[START_COL] - 501) <= TOL
    assert abs(is1[END_COL] - 1402) <= TOL
    assert is1[STRAND_COL] == "-"

    is5 = df.iloc[1]
    assert is5[SEQ_COL] == "tiny"
    assert is5[FAMILY_COL] == "IS5"
    assert abs(is5[START_COL] - 1803) <= TOL
    assert abs(is5[END_COL] - 2997) <= TOL
    assert is5[STRAND_COL] == "+"


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


# --- Unit tests for _calculate_output (no ISEScan execution needed) ---

_TSV_COLS = [SEQ_COL, START_COL, END_COL, STRAND_COL, FAMILY_COL, CLUSTER_COL]


def _write_fake_isescan_tsv(output_dir, seqfile, rows: list[dict]) -> None:
    """Write a fake ISEScan .tsv results file in the expected nested location."""
    nested_dir = output_dir / seqfile.parent.name
    nested_dir.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=_TSV_COLS)
    df.to_csv(nested_dir / f"{seqfile.name}.tsv", sep="\t", index=False)


def _make_step(tiny_genome, tmp_path, output_dir=None):
    output_dir = output_dir or (tmp_path / "isescan_step_output")
    return LocateTNsUsingISEScanStep(
        genome=tiny_genome,
        output_dir=output_dir,
        env_name=None,
        exec_name="isescan.py",
        threads=1,
        force=True,
    )


def test_calculate_output_nan_family(tiny_genome, tmp_path):
    """NaN family should produce tn_name from cluster, not 'nan'."""
    output_dir = tmp_path / "nan_test"
    _write_fake_isescan_tsv(output_dir, tiny_genome.fasta_path, [
        {SEQ_COL: "tiny", START_COL: 100, END_COL: 800, STRAND_COL: "+",
         FAMILY_COL: float("nan"), CLUSTER_COL: "IS1"},
    ])
    step = _make_step(tiny_genome, tmp_path, output_dir)
    result = step._calculate_output()
    tn = list(result.values())[0]
    assert tn.tn_name == "IS1", f"Expected 'IS1', got '{tn.tn_name}'"


def test_calculate_output_nan_both(tiny_genome, tmp_path):
    """NaN family and cluster should produce 'unknown'."""
    output_dir = tmp_path / "nan_both_test"
    _write_fake_isescan_tsv(output_dir, tiny_genome.fasta_path, [
        {SEQ_COL: "tiny", START_COL: 100, END_COL: 800, STRAND_COL: "+",
         FAMILY_COL: float("nan"), CLUSTER_COL: float("nan")},
    ])
    step = _make_step(tiny_genome, tmp_path, output_dir)
    result = step._calculate_output()
    tn = list(result.values())[0]
    assert tn.tn_name == "unknown", f"Expected 'unknown', got '{tn.tn_name}'"


def test_calculate_output_normal(tiny_genome, tmp_path):
    """Normal family+cluster should produce 'family:cluster'."""
    output_dir = tmp_path / "normal_test"
    _write_fake_isescan_tsv(output_dir, tiny_genome.fasta_path, [
        {SEQ_COL: "tiny", START_COL: 100, END_COL: 800, STRAND_COL: "+",
         FAMILY_COL: "IS1", CLUSTER_COL: "IS1A"},
    ])
    step = _make_step(tiny_genome, tmp_path, output_dir)
    result = step._calculate_output()
    tn = list(result.values())[0]
    assert tn.tn_name == "IS1:IS1A"


def test_calculate_output_reverse_strand(tiny_genome, tmp_path):
    """Minus strand should produce REVERSE orientation."""
    from amplifinder.data_types import Orientation
    output_dir = tmp_path / "strand_test"
    _write_fake_isescan_tsv(output_dir, tiny_genome.fasta_path, [
        {SEQ_COL: "tiny", START_COL: 100, END_COL: 800, STRAND_COL: "-",
         FAMILY_COL: "IS5", CLUSTER_COL: "IS5B"},
    ])
    step = _make_step(tiny_genome, tmp_path, output_dir)
    result = step._calculate_output()
    tn = list(result.values())[0]
    assert tn.orientation == Orientation.REVERSE


def test_calculate_output_empty(tiny_genome, tmp_path):
    """Empty results should return empty RecordTypedDf."""
    output_dir = tmp_path / "empty_test"
    _write_fake_isescan_tsv(output_dir, tiny_genome.fasta_path, [])
    step = _make_step(tiny_genome, tmp_path, output_dir)
    result = step._calculate_output()
    assert len(result) == 0

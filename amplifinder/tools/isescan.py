"""ISEScan runner and output parsing."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pandas as pd

from amplifinder.logger import logger
from amplifinder.utils.file_utils import ensure_dir
from amplifinder.utils.run_utils import run_command


def run_isescan(
    seqfile: Path,
    output_dir: Path,
    *,
    env_name: Optional[str],
    exec_name: str = "isescan.py",
    threads: int = 2,
) -> Path:
    """Run ISEScan in a conda environment or current environment.

    Args:
        seqfile: Reference FASTA file path.
        output_dir: Directory where ISEScan writes outputs.
        env_name: Conda environment name containing ISEScan, or None to use current environment.
        exec_name: ISEScan entrypoint (default: isescan.py).
        threads: Threads for ISEScan (--nthread).

    Returns:
        Output directory path.
    """
    seqfile = Path(seqfile)
    output_dir = Path(output_dir)
    ensure_dir(output_dir)

    cmd = [
        exec_name,
        "--seqfile", str(seqfile),
        "--output", str(output_dir),
        "--nthread", str(threads),
    ]
    
    if env_name:
        cmd = ["conda", "run", "-n", env_name] + cmd
        logger.info(f"Running ISEScan in conda env '{env_name}': {exec_name} ...")
    else:
        logger.info(f"Running ISEScan: {exec_name} ...")

    run_command(cmd, check=True, capture_output=True, text=True, error_msg="ISEScan failed", cwd=output_dir)
    return output_dir


def get_isescan_results_file(seqfile: Path, output_dir: Path) -> Optional[Path]:
    """Find ISEScan results file (.tsv preferred, else .csv)."""
    seqfile = Path(seqfile)
    output_dir = Path(output_dir)
    # ISEScan nests outputs under a directory named after the input folder
    nested_dir = output_dir / seqfile.parent.name
    tsv_path = nested_dir / f"{seqfile.name}.tsv"
    if tsv_path.exists():
        return tsv_path
    csv_path = nested_dir / f"{seqfile.name}.csv"
    if csv_path.exists():
        return csv_path
    # Fallback: search recursively for matching result file
    for suffix in (".tsv", ".csv"):
        matches = list(output_dir.rglob(f"{seqfile.name}{suffix}"))
        if matches:
            return matches[0]
    return None


def parse_isescan_results(seqfile: Path, output_dir: Path) -> pd.DataFrame:
    """Parse ISEScan results into a DataFrame.

    Returns:
        DataFrame with ISEScan result rows (may be empty).
    """
    results_file = get_isescan_results_file(seqfile, output_dir)
    if results_file is None:
        raise FileNotFoundError(
            f"ISEScan results not found in {output_dir} for {seqfile.name}"
        )

    sep = "\t" if results_file.suffix == ".tsv" else ","
    df = pd.read_csv(results_file, sep=sep)
    return df

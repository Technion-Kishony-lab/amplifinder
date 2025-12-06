"""breseq runner and output parsing."""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd

from amplifinder.logger import info, warning
from amplifinder.data_types.junction import Junction
from amplifinder.data import load_all_field_defs

if TYPE_CHECKING:
    from amplifinder.data_types.genome import Genome


BRESEQ_DOCKER_IMAGE = "ummidock/breseq:0.32.1"

# Load field definitions from CSV files
RECORD_TYPES = load_all_field_defs()


# =============================================================================
# Runner functions
# =============================================================================

def run_breseq(
    ref_paths: List[Path],
    fastq_path: Path,
    output_path: Path,
    docker: bool = True,
    threads: int = 4,
) -> Path:
    """Run breseq alignment pipeline.
    
    Args:
        ref_paths: List of paths to reference files (.gb or .fasta)
        fastq_path: Path to directory containing FASTQ files
        output_path: Path to output directory
        docker: Use Docker to run breseq (default: True)
        threads: Number of threads for breseq (default: 4)
        
    Returns:
        Path to breseq output directory
        
    Raises:
        RuntimeError: If breseq fails
    """
    output_path = Path(output_path)
    fastq_path = Path(fastq_path)
    
    # Check if already completed
    output_gd = output_path / "output" / "output.gd"
    if output_gd.exists():
        info(f"breseq output already exists: {output_gd}")
        return output_path
    
    # Find FASTQ files
    fastq_files = list(fastq_path.glob("*.fastq*"))
    if not fastq_files:
        raise FileNotFoundError(f"No FASTQ files found in {fastq_path}")
    
    info(f"Found {len(fastq_files)} FASTQ file(s)")
    
    # Clean output directory if partial run exists
    if output_path.exists():
        warning(f"Cleaning incomplete breseq output: {output_path}")
        subprocess.run(["rm", "-rf", str(output_path)], check=True)
    
    output_path.mkdir(parents=True, exist_ok=True)
    
    if docker:
        _run_breseq_docker(ref_paths, fastq_path, fastq_files, output_path, threads)
    else:
        _run_breseq_local(ref_paths, fastq_files, output_path, threads)
    
    # Verify output
    if not output_gd.exists():
        raise RuntimeError(f"breseq failed: {output_gd} not created")
    
    info("breseq completed successfully")
    return output_path


def _run_breseq_docker(
    ref_paths: List[Path],
    fastq_path: Path,
    fastq_files: List[Path],
    output_path: Path,
    threads: int,
) -> None:
    """Run breseq using Docker."""
    # Get absolute paths for Docker mounts
    ref_dir = ref_paths[0].parent.resolve()
    fastq_dir = fastq_path.resolve()
    out_dir = output_path.resolve()
    
    # Build Docker command
    # Mount reference directory
    mounts = [
        "-v", f"{ref_dir}:/ref:ro",
        "-v", f"{fastq_dir}:/fastq:ro",
        "-v", f"{out_dir}:/out",
    ]
    
    # Build reference arguments
    ref_args = []
    for ref_path in ref_paths:
        ref_args.extend(["-r", f"/ref/{ref_path.name}"])
    
    # Build FASTQ arguments
    fastq_args = [f"/fastq/{f.name}" for f in fastq_files]
    
    cmd = [
        "docker", "run", "--rm", "-u", "root", "-i",
        *mounts,
        BRESEQ_DOCKER_IMAGE,
        "breseq",
        "-j", str(threads),
        "-o", "/out",
        *ref_args,
        *fastq_args,
    ]
    
    info(f"Running breseq (Docker): {' '.join(cmd[:10])}...")
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        warning(f"breseq stderr: {result.stderr}")
        raise RuntimeError(f"breseq failed with exit code {result.returncode}")


def _run_breseq_local(
    ref_paths: List[Path],
    fastq_files: List[Path],
    output_path: Path,
    threads: int,
) -> None:
    """Run breseq locally (requires breseq in PATH)."""
    # Build reference arguments
    ref_args = []
    for ref_path in ref_paths:
        ref_args.extend(["-r", str(ref_path.resolve())])
    
    # Build FASTQ arguments
    fastq_args = [str(f.resolve()) for f in fastq_files]
    
    cmd = [
        "breseq",
        "-j", str(threads),
        "-o", str(output_path.resolve()),
        *ref_args,
        *fastq_args,
    ]
    
    info(f"Running breseq (local): {' '.join(cmd[:8])}...")
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
    )
    
    if result.returncode != 0:
        warning(f"breseq stderr: {result.stderr}")
        raise RuntimeError(f"breseq failed with exit code {result.returncode}")


def get_ref_file(genome: "Genome", use_annotations: bool = True) -> Path:
    """Get reference file for breseq.
    
    Args:
        genome: Genome object
        use_annotations: If True, prefer GenBank (has annotations);
                        If False, prefer FASTA (for ISfinder mode)
    
    Returns:
        Path to reference file (.gb or .fasta)
    """
    if use_annotations and genome.genbank_path:
        return genome.genbank_path
    elif genome.fasta_path:
        return genome.fasta_path
    elif genome.genbank_path:
        return genome.genbank_path
    else:
        raise ValueError(f"No reference file available for {genome.name}")


def get_breseq_version(breseq_path: Path) -> Optional[str]:
    """Extract breseq version from output directory.
    
    Args:
        breseq_path: Path to breseq output directory
        
    Returns:
        Version string or None if not found
    """
    log_file = breseq_path / "output" / "log.txt"
    if not log_file.exists():
        return None
    
    with open(log_file) as f:
        for line in f:
            if "breseq" in line.lower() and "version" in line.lower():
                return line.strip()
    
    return None


# =============================================================================
# Parser functions
# =============================================================================

def parse_breseq_output(breseq_path: Path) -> Dict[str, pd.DataFrame]:
    """Parse breseq output.gd file into DataFrames.
    
    Args:
        breseq_path: Path to breseq output directory
        
    Returns:
        Dictionary with keys: JC, SNP, MOB, DEL, UN
        Each value is a DataFrame with parsed records
    """
    output_gd = Path(breseq_path) / "output" / "output.gd"
    
    if not output_gd.exists():
        raise FileNotFoundError(f"breseq output not found: {output_gd}")
    
    # Initialize empty record lists
    records: Dict[str, List[Dict[str, Any]]] = {
        name: [] for name in RECORD_TYPES.keys()
    }
    
    info(f"Parsing breseq output: {output_gd}")
    
    with open(output_gd) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            
            parts = line.split("\t")
            record_type = parts[0]
            
            if record_type not in RECORD_TYPES:
                continue
            
            record = _parse_record(parts, RECORD_TYPES[record_type])
            records[record_type].append(record)
    
    # Convert to DataFrames
    results = {}
    for name, recs in records.items():
        df = pd.DataFrame(recs) if recs else pd.DataFrame()
        results[name] = df
        if not df.empty:
            info(f"  {name}: {len(df)} records")
    
    return results


def _parse_record(parts: List[str], field_defs: pd.DataFrame) -> Dict[str, Any]:
    """Parse a single breseq record line.
    
    Args:
        parts: Tab-separated parts of the line
        field_defs: DataFrame with columns: fields, types, optional
        
    Returns:
        Dictionary with parsed fields
    """
    record: Dict[str, Any] = {}
    
    # Get ordered non-optional fields
    non_optional = field_defs[~field_defs["optional"]]
    
    # Parse positional (non-optional) fields
    for i, (_, row) in enumerate(non_optional.iterrows()):
        if i < len(parts):
            record[row["fields"]] = _convert_value(parts[i], row["types"])
    
    # Parse optional key=value fields
    n_non_optional = len(non_optional)
    type_lookup = field_defs.set_index("fields")["types"]
    for part in parts[n_non_optional:]:
        if "=" in part:
            key, value = part.split("=", 1)
            if key in type_lookup.index:
                record[key] = _convert_value(value, type_lookup[key])
            else:
                record[key] = value
    
    return record


def _convert_value(value: str, dtype: str) -> Any:
    """Convert string value to appropriate type.
    
    Args:
        value: String value from breseq output
        dtype: Target type (string, int, float)
        
    Returns:
        Converted value
    """
    if value == "NA" or value == "":
        return None
    
    try:
        if dtype == "string":
            return value
        elif dtype == "int":
            return int(float(value))  # Handle "1.0" -> 1
        elif dtype == "float":
            return float(value)
    except (ValueError, TypeError):
        return value
    
    return value


def parse_coverage(breseq_path: Path, ref_names: List[str]) -> Dict[str, np.ndarray]:
    """Parse coverage files from breseq output.
    
    Args:
        breseq_path: Path to breseq output directory
        ref_names: List of reference scaffold names
        
    Returns:
        Dictionary mapping scaffold name to coverage array
    """
    breseq_path = Path(breseq_path)
    coverage_dir = breseq_path / "08_mutation_identification"
    
    coverage = {}
    
    for ref_name in ref_names:
        cov_file = coverage_dir / f"{ref_name}.coverage.tab"
        
        if not cov_file.exists():
            warning(f"Coverage file not found: {cov_file}")
            continue
        
        info(f"Parsing coverage for {ref_name}")
        
        # Read coverage table
        df = pd.read_csv(
            cov_file,
            sep="\t",
            usecols=[0, 1],  # unique_top_cov, unique_bot_cov
            skiprows=1,
            names=["top", "bot"],
            dtype={"top": np.int32, "bot": np.int32},
        )
        
        # Total coverage is sum of top and bottom strand
        coverage[ref_name] = df["top"].values + df["bot"].values
    
    return coverage


def get_breseq_summary(breseq_path: Path) -> Dict[str, Any]:
    """Extract summary statistics from breseq output.
    
    Args:
        breseq_path: Path to breseq output directory
        
    Returns:
        Dictionary with summary statistics including mapped_bases
    """
    summary_file = Path(breseq_path) / "05_alignment_correction" / "summary.json"
    
    if not summary_file.exists():
        warning(f"Summary file not found: {summary_file}")
        return {}
    
    with open(summary_file) as f:
        data = json.load(f)
    
    # Extract mapped bases from summary
    summary = {}
    
    if "references" in data:
        mapped_bases = []
        for ref in data["references"].values():
            if "bases_mapped_to_reference" in ref:
                mapped_bases.append(ref["bases_mapped_to_reference"])
        summary["mapped_bases"] = mapped_bases
        summary["total_mapped_bases"] = sum(mapped_bases)
    
    return summary


def get_read_length(breseq_path: Path, n_reads: int = 10) -> List[int]:
    """Extract read lengths from FASTQ files in breseq data directory.
    
    Args:
        breseq_path: Path to breseq output directory
        n_reads: Number of reads to sample (default: 10)
        
    Returns:
        List of read lengths
    """
    data_dir = Path(breseq_path) / "data"
    
    if not data_dir.exists():
        warning(f"Data directory not found: {data_dir}")
        return []
    
    fastq_files = list(data_dir.glob("*.fastq"))
    if not fastq_files:
        fastq_files = list(data_dir.glob("*.fq"))
    
    if not fastq_files:
        warning(f"No FASTQ files found in {data_dir}")
        return []
    
    # Read first few sequences from first file
    read_lengths = []
    with open(fastq_files[0]) as f:
        for i, line in enumerate(f):
            if i % 4 == 1:  # Sequence line
                read_lengths.append(len(line.strip()))
                if len(read_lengths) >= n_reads:
                    break
    
    return read_lengths


def get_junctions(breseq_path: Path) -> List[Junction]:
    """Parse and return Junction objects from breseq output.
    
    Args:
        breseq_path: Path to breseq output directory
        
    Returns:
        List of Junction objects
    """
    results = parse_breseq_output(breseq_path)
    jc_df = results.get("JC", pd.DataFrame())
    
    if jc_df.empty:
        return []
    
    junctions = []
    for _, row in jc_df.iterrows():
        junctions.append(Junction.from_dict(row.to_dict()))
    
    return junctions


def count_output_lines(breseq_path: Path) -> int:
    """Count lines in breseq output.gd file.
    
    Args:
        breseq_path: Path to breseq output directory
        
    Returns:
        Number of lines in output.gd
    """
    output_gd = Path(breseq_path) / "output" / "output.gd"
    
    if not output_gd.exists():
        return 0
    
    with open(output_gd) as f:
        return sum(1 for _ in f)


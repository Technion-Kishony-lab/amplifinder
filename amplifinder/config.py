"""Configuration management for AmpliFinder."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Tuple, ClassVar
import json
import yaml

from amplifinder.utils.file_utils import ensure_dir


# Default configuration values (from config.txt)
DEFAULT_CONFIG = {
    # External tools
    "breseq_docker": True,
    "blastn_path": None,  # Auto-detect from PATH
    "samtools_path": None,  # Auto-detect from PATH
    "isdb_path": None,  # Use bundled ISfinderDB

    # IS/TN detection parameters
    "max_dist_to_IS": 10,
    "trim_jc_flanking": 5,
    "length_seq_into_is": 200,
    "reference_IS_out_span": 100,
    "isfinder_evalue": 1e-4,
    "isfinder_critical_coverage": 0.9,

    # breseq parameters
    "breseq_threads": 4,

    # Junction filtering
    "min_jct_cov": 5,
    "min_amplicon_length": 30,
    "max_amplicon_length": 1_000_000,
    "filter_amplicon_length": 100,

    # Coverage parameters
    "ncp_limit1": -1,
    "ncp_limit2": 3,
    "ncp_n": 150,

    # Copy number thresholds
    "copy_number_threshold": 1.5,
    "del_copy_number_threshold": 0.3,

    # Alignment parameters
    "min_overlap": 12,
    "score_min": (0, -0.1),
    "mismatch_penalty": (5, 5),

    # File size thresholds
    "breseq_output_size_threshold": 10_000,
    "max_fastq_size": 500_000_000,
    "min_num_bases": 80_000_000,

    # Filtering options
    "shortest_amplicon": 1,
    "true_transposition_length": 0,
    "remove_jc_breseq_reject": False,
    "remove_isjc2_breseq_reject": True,

    # Logging
    "log_path": "amplifinder.log",
}


@dataclass
class Config:
    """AmpliFinder configuration."""

    RUN_CONFIG_FILENAME: ClassVar[str] = "run_config.yaml"

    # Required paths
    iso_path: Path
    ref_name: str

    # Optional paths
    anc_path: Optional[Path] = None
    output_dir: Path = field(default_factory=lambda: Path("output"))
    ref_path: Path = field(default_factory=lambda: Path("genomesDB"))
    iso_breseq_path: Optional[Path] = None
    anc_breseq_path: Optional[Path] = None

    # Read length (if None, calculated from FASTQ)
    iso_read_length: Optional[int] = None
    anc_read_length: Optional[int] = None

    # Names (derived from paths if not provided)
    iso_name: Optional[str] = None
    anc_name: Optional[str] = None

    # Reference options
    ncbi: bool = True
    use_isfinder: bool = False

    # External tools
    breseq_docker: bool = True
    blastn_path: Optional[Path] = None
    samtools_path: Optional[Path] = None
    isdb_path: Optional[Path] = None

    # IS/TN detection parameters
    max_dist_to_IS: int = 10
    trim_jc_flanking: int = 5
    length_seq_into_is: int = 200
    reference_IS_out_span: int = 100
    isfinder_evalue: float = 1e-4
    isfinder_critical_coverage: float = 0.9

    # breseq parameters
    breseq_threads: int = 4

    # Junction filtering
    min_jct_cov: int = 5
    min_amplicon_length: int = 30
    max_amplicon_length: int = 1_000_000
    filter_amplicon_length: int = 100

    # Coverage parameters
    ncp_limit1: int = -1
    ncp_limit2: int = 3
    ncp_n: int = 150

    # Copy number thresholds
    copy_number_threshold: float = 1.5
    del_copy_number_threshold: float = 0.3

    # Alignment parameters
    min_overlap: int = 12
    score_min: Tuple[float, float] = (0, -0.1)
    mismatch_penalty: Tuple[int, int] = (5, 5)

    # File size thresholds
    breseq_output_size_threshold: int = 10_000
    max_fastq_size: int = 500_000_000
    min_num_bases: int = 80_000_000

    # Filtering options
    shortest_amplicon: int = 1
    true_transposition_length: int = 0
    remove_jc_breseq_reject: bool = False
    remove_isjc2_breseq_reject: bool = True

    # Logging
    log_path: str = "amplifinder.log"

    @property
    def has_ancestor(self) -> bool:
        """True if ancestor comparison should be performed."""
        return self.anc_path is not None

    @property
    def iso_run_dir(self) -> Path:
        """Get the isolate run directory path.
        
        Structure: {output_dir}/{ref_name}/{anc_name}/{iso_name}/
        """
        return self.output_dir / self.ref_name / self.anc_name / self.iso_name

    @property
    def anc_run_dir(self) -> Path:
        """Get the ancestor run directory path.
        
        Structure: {output_dir}/{ref_name}/{anc_name}/{anc_name}/
        """
        return self.output_dir / self.ref_name / self.anc_name / self.anc_name

    def get_anc_breseq_path(self) -> Path:
        """Return ancestor breseq path (provided or default run dir)."""
        if not self.has_ancestor:
            raise ValueError("Ancestor breseq path requested but no ancestor configured.")
        return self.anc_breseq_path or (self.anc_run_dir / "breseq")

    def get_iso_breseq_path(self) -> Path:
        """Return isolate breseq path (provided or default run dir)."""
        return self.iso_breseq_path or (self.iso_run_dir / "breseq")

    def save(self, run_dir: Path) -> None:
        """Save config to run_config.yaml in run directory."""
        run_dir = ensure_dir(run_dir)
        config_path = run_dir / self.RUN_CONFIG_FILENAME
        
        # Convert Config to dict, handling Path objects and tuples
        config_dict = {}
        for key, value in self.__dict__.items():
            if isinstance(value, Path):
                config_dict[key] = str(value)
            elif isinstance(value, tuple):
                config_dict[key] = list(value)
            elif value is not None:
                config_dict[key] = value
        
        with open(config_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
        
        from amplifinder.logger import info
        info(f"Saved config to: {config_path}")

    @classmethod
    def load_from_run(cls, run_dir: Path) -> "Config":
        """Load config from run_config.yaml in run directory.
        
        Raises:
            FileNotFoundError: If run_config.yaml doesn't exist
        """
        run_dir = Path(run_dir)
        config_path = run_dir / cls.RUN_CONFIG_FILENAME
        
        if not config_path.exists():
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        config_dict = load_config(config_path)
        
        # Convert string paths back to Path objects
        path_fields = ['iso_path', 'anc_path', 'output_dir', 'ref_path', 
                       'iso_breseq_path', 'anc_breseq_path', 'blastn_path', 
                       'samtools_path', 'isdb_path']
        for field in path_fields:
            if field in config_dict and config_dict[field] is not None:
                config_dict[field] = Path(config_dict[field])
        
        # Convert lists back to tuples for tuple fields
        tuple_fields = ['score_min', 'mismatch_penalty']
        for field in tuple_fields:
            if field in config_dict and isinstance(config_dict[field], list):
                config_dict[field] = tuple(config_dict[field])
        
        return cls(**config_dict)

    def __post_init__(self):
        """Convert paths and set derived values."""
        # Convert string paths to Path objects
        self.iso_path = Path(self.iso_path)
        self.output_dir = Path(self.output_dir)
        self.ref_path = Path(self.ref_path)

        if self.anc_path is not None:
            self.anc_path = Path(self.anc_path)
        if self.iso_breseq_path is not None:
            self.iso_breseq_path = Path(self.iso_breseq_path)
        if self.anc_breseq_path is not None:
            self.anc_breseq_path = Path(self.anc_breseq_path)
        if self.blastn_path is not None:
            self.blastn_path = Path(self.blastn_path)
        if self.samtools_path is not None:
            self.samtools_path = Path(self.samtools_path)
        if self.isdb_path is not None:
            self.isdb_path = Path(self.isdb_path)

        # Derive names from paths if not provided
        if self.iso_name is None:
            self.iso_name = self.iso_path.stem or "isolate"
        if self.anc_name is None:
            if self.anc_path is not None:
                self.anc_name = self.anc_path.stem or "ancestor"
            else:
                # No ancestor: isolate is its own "ancestor" for folder structure
                self.anc_name = self.iso_name

        # Convert tuples if passed as lists
        if isinstance(self.score_min, list):
            self.score_min = tuple(self.score_min)
        if isinstance(self.mismatch_penalty, list):
            self.mismatch_penalty = tuple(self.mismatch_penalty)

        # Validate: ISfinder required for non-NCBI genomes
        if not self.ncbi and not self.use_isfinder:
            raise ValueError(
                "ISfinder is required for local (non-NCBI) genomes. "
                "Use --use-isfinder flag."
            )

        # Validate: genomesDB is reserved for NCBI references
        if str(self.ref_path) == "genomesDB" and not self.ncbi:
            raise ValueError(
                "genomesDB directory is reserved for NCBI references. "
                "Choose a different --ref-path for local genomes."
            )


def load_config(config_path: Path) -> dict[str, Any]:
    """Load configuration from YAML or JSON file.

    Args:
        config_path: Path to config file (.yaml, .yml, or .json)

    Returns:
        Dictionary of configuration values
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    suffix = config_path.suffix.lower()

    with open(config_path) as f:
        if suffix in (".yaml", ".yml"):
            return yaml.safe_load(f) or {}
        elif suffix == ".json":
            return json.load(f)
        else:
            raise ValueError(
                f"Unsupported config format: {suffix}. "
                "Use .yaml, .yml, or .json"
            )


def merge_config(
    cli_args: dict[str, Any],
    config_file: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Merge configuration from defaults, config file, and CLI args.

    Priority: CLI args > config file > defaults

    Args:
        cli_args: Arguments from command line
        config_file: Loaded config file (optional)

    Returns:
        Merged configuration dictionary
    """
    merged = DEFAULT_CONFIG.copy()

    # Apply config file values
    if config_file:
        for key, value in config_file.items():
            if value is not None:
                merged[key] = value

    # Apply CLI args (highest priority)
    for key, value in cli_args.items():
        if value is not None:
            merged[key] = value

    return merged



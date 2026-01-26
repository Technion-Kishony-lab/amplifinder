"""Configuration management for AmpliFinder."""

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Optional, Tuple, ClassVar, Union
import json
import yaml
from pydantic import BaseModel, ConfigDict

from amplifinder.data_types import AverageMethod
from amplifinder.utils.file_utils import ensure_dir


class FrozenParams(BaseModel):
    """Base class for frozen parameter models."""
    model_config = ConfigDict(frozen=True)


class AlignmentFilterParams(FrozenParams):
    """Parameters for read filtering."""

    max_nm_score: Optional[int] = 3  # 3,
    min_as_score: Optional[int] = -25  # -25
    length_tolerance: float = 0.1


class AlignmentClassifyParams(FrozenParams):
    """Parameters for junction alignment analysis."""

    min_overlap_len: int = 13

    # Max distance of the far end of the hit from the junction
    # None to allow any distance
    # Positive value to allow reads within the specified distance from the junction
    # Negative value to specify as distance from the jc arm end
    max_dist_from_junction: Optional[int] = -10

    def get_max_dist_from_junction(self, arm_len: int) -> int:
        max_dist = self.max_dist_from_junction
        if max_dist is None or max_dist > arm_len:
            return arm_len
        if max_dist < 0:
            return arm_len + max_dist
        return max_dist


class BowtieParams(FrozenParams):
    """Parameters for bowtie2 alignment."""

    score_min: Union[str, Tuple[float, float]] = (0, -0.2)
    mismatch_penalty: Union[str, Tuple[int, int]] = (5, 5)
    local: bool = False
    num_alignments: int = 100
    min_qlen: Optional[int] = 136


class JcCallParams(FrozenParams):
    """Parameters for junction coverage calling (exists, not exists, ambiguous)."""

    # for a jct to be negative, the number of spanning reads should be less than:
    neg_threshold_abs: int = 5  # absolute number of spanning reads, or
    neg_threshold_rel: float = 0.01  # fraction of the expected number of spanning reads

    # for a jct to be positive, the number of spanning reads should be
    # greater than the expected number minus x standard deviations
    pos_threshold_in_num_std_below_expected: int = 3
    pos_threshold_rel: float = 0.4  # non-stochastic copy number fluctuations


# Default configuration values (from config.txt and added alignment params)
DEFAULT_CONFIG = {
    # External tools
    "breseq_docker": True,
    "blastn_path": None,  # Auto-detect from PATH
    "samtools_path": None,  # Auto-detect from PATH
    "isdb_path": None,  # Use bundled ISfinderDB

    # IS/TN detection parameters
    "max_dist_to_IS": 10,
    "trim_jc_flanking": 5,
    "reference_IS_in_span": None,
    "reference_IS_out_span": 100,
    "isfinder_evalue": 1e-4,
    "isfinder_critical_coverage": 0.9,

    # Alignment threads (breseq, bowtie2)
    "threads": 4,

    # Junction filtering
    "min_amplicon_length": 30,
    "max_amplicon_length": 1_000_000,
    "filter_amplicon_length": 100,

    # Coverage parameters
    "ncp_min": 0.1,
    "ncp_max": 1000.0,
    "ncp_n": 150,
    "average_method": "median",  # 'median', 'mode', or 'mean' (converted to AverageMethod enum)

    # Copy number thresholds
    "copy_number_threshold": 1.5,
    "del_copy_number_threshold": 0.3,

    # File size thresholds
    "breseq_output_size_threshold": 10_000,  # Terminates if breseq output.gd exceeds this line count
    # TODO: These are not used yet
    "max_fastq_size": 500_000_000,
    "min_num_bases": 80_000_000,

    # Filtering options
    # TODO: These are not used yet
    "shortest_amplicon": 1,
    "true_transposition_length": 0,
    "remove_jc_breseq_reject": False,
    "remove_isjc2_breseq_reject": True,

    # Logging
    "log_path": "amplifinder.log",

    # Plotting
    "create_plots": True,
}


def _normalize_alignment_params(cfg: dict[str, Any]) -> None:
    """Move legacy scalar alignment params into param objects."""
    def _collect(keys: tuple[str, ...]) -> dict[str, Any]:
        return {k: cfg.pop(k) for k in keys if k in cfg}

    param_specs = [
        (AlignmentFilterParams, "alignment_filter_params"),
        (AlignmentClassifyParams, "alignment_analysis_params"),
        (BowtieParams, "bowtie_params"),
        (JcCallParams, "jc_call_params"),
    ]
    for param_class, param_name in param_specs:
        keys = tuple(param_class.__annotations__.keys())
        vals = _collect(keys)
        cfg[param_name] = vals or cfg.get(param_name) or {}


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
    reference_IS_in_span: Optional[int] = None  # None = use the entire IS element
    reference_IS_out_span: int = 100
    isfinder_evalue: float = 1e-4
    isfinder_critical_coverage: float = 0.9

    # Alignment threads (breseq, bowtie2)
    threads: int = 4

    # Junction filtering
    min_amplicon_length: int = 30
    max_amplicon_length: int = 1_000_000
    filter_amplicon_length: int = 100

    # Coverage parameters
    ncp_min: float = 0.1
    ncp_max: float = 1000.0
    ncp_n: int = 150
    average_method: AverageMethod = AverageMethod.MEDIAN

    # Copy number thresholds
    copy_number_threshold: float = 1.5
    del_copy_number_threshold: float = 0.3

    # Alignment parameters
    alignment_filter_params: Optional[AlignmentFilterParams] = None
    alignment_analysis_params: Optional[AlignmentClassifyParams] = None
    bowtie_params: Optional[BowtieParams] = None
    jc_call_params: Optional[JcCallParams] = None

    # File size thresholds
    breseq_output_size_threshold: int = 10_000  # Terminates if breseq output.gd exceeds this line count
    # TODO: These are not used yet
    max_fastq_size: int = 500_000_000
    min_num_bases: int = 80_000_000

    # Filtering options
    # TODO: These are not used yet
    shortest_amplicon: int = 1
    true_transposition_length: int = 0
    remove_jc_breseq_reject: bool = False
    remove_isjc2_breseq_reject: bool = True

    # Logging
    log_path: str = "amplifinder.log"

    # Plotting
    create_plots: bool = True

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

    def get_anc_breseq_path(self) -> Optional[Path]:
        """Return ancestor breseq path (provided or default run dir).

        Returns:
            Path to ancestor breseq directory, or None if no ancestor configured.
        """
        if not self.has_ancestor:
            return None
        return self.anc_breseq_path or (self.anc_run_dir / "breseq")

    def get_iso_breseq_path(self) -> Path:
        """Return isolate breseq path (provided or default run dir)."""
        return self.iso_breseq_path or (self.iso_run_dir / "breseq")

    def get_anc_name(self) -> Optional[str]:
        """Return ancestor name, or None if no ancestor configured.

        Note: anc_name property is set to iso_name when no ancestor (for folder structure),
        but this method returns None for semantic correctness when passing to steps.

        Returns:
            Ancestor name, or None if no ancestor.
        """
        if not self.has_ancestor:
            return None
        return self.anc_name

    def get_breseq_paths(self) -> tuple[Path, Optional[Path]]:
        """Get isolate and ancestor breseq paths.
        """
        return self.get_iso_breseq_path(), self.get_anc_breseq_path()

    def to_yaml_dict(self) -> dict[str, Any]:
        """Convert config to YAML-serializable dict."""
        def convert_value(val):
            """Recursively convert values for YAML serialization."""
            if isinstance(val, Path):
                return str(val)
            elif isinstance(val, Enum):
                return val.value
            elif isinstance(val, tuple):
                return [convert_value(v) for v in val]
            elif isinstance(val, dict):
                return {k: convert_value(v) for k, v in val.items()}
            elif isinstance(val, list):
                return [convert_value(v) for v in val]
            elif isinstance(val, BaseModel):
                # pydantic v1/v2 compatibility
                dump_fn = getattr(val, "model_dump", None)
                dumped = dump_fn() if dump_fn else val.dict()
                return convert_value(dumped)
            else:
                return val

        # Convert Config to dict, handling Path objects and tuples
        config_dict = {}
        for key, value in self.__dict__.items():
            if value is not None:
                config_dict[key] = convert_value(value)
        
        return config_dict

    def save_to_file(self, config_path: Path, log: bool = True) -> None:
        """Save config to specified path.
        
        Args:
            config_path: Path to save config file
            log: Whether to log the save location (default: True)
        """
        config_path = Path(config_path)
        config_path.parent.mkdir(parents=True, exist_ok=True)
        
        config_dict = self.to_yaml_dict()
        
        with open(config_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False, sort_keys=False)
        
        if log:
            from amplifinder.logger import info
            info(f"Saved config to:\n{config_path}")
    
    def save(self, run_dir: Path) -> None:
        """Save config to run_config.yaml in run directory."""
        run_dir = ensure_dir(run_dir)
        config_path = run_dir / self.RUN_CONFIG_FILENAME
        self.save_to_file(config_path, log=True)

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
        _normalize_alignment_params(config_dict)

        # Convert string paths back to Path objects
        path_fields = ['iso_path', 'anc_path', 'output_dir', 'ref_path',
                       'iso_breseq_path', 'anc_breseq_path', 'blastn_path',
                       'samtools_path', 'isdb_path']
        for path_field in path_fields:
            if path_field in config_dict and config_dict[path_field] is not None:
                config_dict[path_field] = Path(config_dict[path_field])

        # Convert lists back to tuples for tuple fields
        tuple_fields = ['score_min', 'mismatch_penalty']
        for tuple_field in tuple_fields:
            if tuple_field in config_dict and isinstance(config_dict[tuple_field], list):
                config_dict[tuple_field] = tuple(config_dict[tuple_field])

        return cls(**config_dict)

    def __post_init__(self):
        """Convert paths and set derived values."""
        # Convert string paths to absolute Path objects
        self.iso_path = Path(self.iso_path).resolve()
        self.output_dir = Path(self.output_dir).resolve()
        self.ref_path = Path(self.ref_path).resolve()

        if self.anc_path is not None:
            self.anc_path = Path(self.anc_path).resolve()
        if self.iso_breseq_path is not None:
            self.iso_breseq_path = Path(self.iso_breseq_path).resolve()
        if self.anc_breseq_path is not None:
            self.anc_breseq_path = Path(self.anc_breseq_path).resolve()
        if self.blastn_path is not None:
            self.blastn_path = Path(self.blastn_path).resolve()
        if self.samtools_path is not None:
            self.samtools_path = Path(self.samtools_path).resolve()
        if self.isdb_path is not None:
            self.isdb_path = Path(self.isdb_path).resolve()

        # Derive names from paths if not provided
        if self.iso_name is None:
            self.iso_name = self.iso_path.stem or "isolate"
        if self.anc_name is None:
            if self.anc_path is not None:
                self.anc_name = self.anc_path.stem or "ancestor"
            else:
                # No ancestor: isolate is its own "ancestor" for folder structure
                self.anc_name = self.iso_name

        # Normalize params
        def validate(cls, obj):
            return (cls.model_validate(obj or {}) if hasattr(cls, "model_validate")
                    else cls.parse_obj(obj or {}))
        self.alignment_filter_params = validate(AlignmentFilterParams, self.alignment_filter_params)
        self.alignment_analysis_params = validate(AlignmentClassifyParams, self.alignment_analysis_params)
        self.bowtie_params = validate(BowtieParams, self.bowtie_params)
        self.jc_call_params = validate(JcCallParams, self.jc_call_params)

        # Validate: ISfinder required for non-NCBI genomes
        if not self.ncbi and not self.use_isfinder:
            raise ValueError(
                "ISfinder is required for local (non-NCBI) genomes. "
                "Use --use-isfinder flag."
            )

        # Validate: genomesDB is reserved for NCBI references
        if self.ref_path.name == "genomesDB" and not self.ncbi:
            raise ValueError(
                "genomesDB directory is reserved for NCBI references. "
                "Choose a different --ref-path for local genomes."
            )

        # Validate: average_method (convert string to enum if needed)
        if isinstance(self.average_method, str):
            try:
                self.average_method = AverageMethod(self.average_method)
            except ValueError:
                valid_options = ", ".join(f"'{m.value}'" for m in AverageMethod)
                raise ValueError(
                    f"average_method must be one of {valid_options}, got: {self.average_method}"
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

    _normalize_alignment_params(merged)
    return merged

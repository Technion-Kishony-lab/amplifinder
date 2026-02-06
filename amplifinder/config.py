"""Configuration management for AmpliFinder."""

from dataclasses import dataclass, field, fields, is_dataclass, MISSING
from enum import Enum
from pathlib import Path
from typing import Any, Optional, Tuple, ClassVar, Union

from amplifinder.data_types import AverageMethod
from amplifinder.utils.file_utils import ensure_dir
from amplifinder.utils.yaml_utils import load_config, save_annotated_yaml


@dataclass(frozen=True)
class AlignmentFilterParams:
    """Parameters for read filtering."""
    max_nm_score: Optional[int] = field(default=3, metadata={
        "comment": "max edit distance (NM tag)"})
    min_as_score: Optional[int] = field(default=-25, metadata={
        "comment": "min alignment score (AS tag)"})
    filter_len_tolerance: float = 0.1


@dataclass(frozen=True)
class AlignmentClassifyParams:
    """Parameters for junction alignment analysis."""
    min_overlap_len: int = field(default=13, metadata={
        "comment": "min overlap with junction"})
    max_dist_from_junction: Optional[int] = field(default=-10, metadata={
        "comment": ">0: distance from junction, <0: distance from arm end"})

    def get_max_dist_from_junction(self, arm_len: int) -> int:
        max_dist = self.max_dist_from_junction
        if max_dist is None or max_dist > arm_len:
            return arm_len
        if max_dist < 0:
            return arm_len + max_dist
        return max_dist


@dataclass(frozen=True)
class BowtieParams:
    """Parameters for bowtie2 alignment."""
    score_min: Union[str, Tuple[float, float]] = field(default=(0, -0.2), metadata={
        "comment": "bowtie2 --score-min (default: 0, -0.2)"})
    mismatch_penalty: Union[str, Tuple[int, int]] = field(
        default=(5, 5), metadata={"comment": "bowtie2 --mp"})
    local: bool = field(default=False, metadata={
        "comment": "true: local alignment mode, false: end-to-end alignment mode"})
    num_alignments: int = field(default=100, metadata={
        "comment": "bowtie2 -k (default: 100)"})
    bowtie_len_tolerance: Optional[float] = field(default=0.1, metadata={
        "comment": "bowtie2 len filter relative to read length (default: 0.1)"})


@dataclass(frozen=True)
class JcCallParams:
    """Parameters for junction coverage calling (exists, not exists, ambiguous)."""

    # for a jct to be negative, the number of spanning reads should be less than:
    neg_threshold_abs: int = field(default=5, metadata={
        "comment": "absolute threshold for negative call"})
    neg_threshold_rel: float = field(default=0.01, metadata={
        "comment": "fraction of expected spanning reads"})

    # for a jct to be positive, the number of spanning reads should be
    # greater than the expected number minus x standard deviations
    pos_threshold_in_num_std_below_expected: int = field(default=3, metadata={
        "comment": "std devs below expected for positive call"})
    pos_threshold_rel: float = field(default=0.4, metadata={
        "comment": "non-stochastic relative fluctuation margin"})


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


@dataclass(kw_only=True)
class Config:
    """AmpliFinder configuration."""

    RUN_CONFIG_FILENAME: ClassVar[str] = "run_config.yaml"

    # --- Reference ---
    ref_name: Optional[str] = field(default=None, metadata={
        "section": "Reference",
        "comment": "NCBI accession for reference genome"})
    ref_path: Path = field(default_factory=lambda: Path("genomesDB"), metadata={
        "comment": "path to reference genome folder (default: genomesDB)"})
    ncbi: bool = field(default=True, metadata={
        "comment": "fetch reference from NCBI if not found in ref_path"})

    # --- Isolate ---
    iso_fastq_path: Optional[Path] = field(default=None, metadata={
        "section": "Isolate",
        "comment": "path to isolate FASTQ file(s)"})
    iso_name: Optional[str] = field(default=None, metadata={
        "comment": "isolate name (default: derived from iso_fastq_path)"})
    iso_read_length: Optional[int] = field(default=None, metadata={
        "comment": "read length (if not provided, auto-detected from FASTQ)"})
    iso_breseq_path: Optional[Path] = field(default=None, metadata={
        "comment": "path to create/load breseq (default: <iso_run_folder>/breseq)"})

    # --- Ancestor (optional) ---
    anc_fastq_path: Optional[Path] = field(default=None, metadata={
        "section": "Ancestor (optional)",
        "comment": "path to ancestor FASTQ file(s)"})
    anc_name: Optional[str] = field(default=None, metadata={
        "comment": "ancestor name (default: derived from anc_fastq_path)"})
    anc_read_length: Optional[int] = field(default=None, metadata={
        "comment": "read length (if not provided, auto-detected from FASTQ)"})
    anc_breseq_path: Optional[Path] = field(default=None, metadata={
        "comment": "path to create/load breseq (default: <anc_run_folder>/breseq)"})

    # --- Output ---
    output_dir: Path = field(default_factory=lambda: Path("output"), metadata={
        "comment": "output directory (default: output)"})

    # --- IS-finder ---
    use_isfinder: bool = field(default=False, metadata={
        "section": "IS-finder",
        "comment": "use ISfinder database for IS detection (default: False)"})
    isfinder_evalue: float = 1e-4
    isfinder_critical_coverage: float = 0.9

    # --- Detection parameters ---
    reference_IS_in_span: Optional[int] = field(default=None, metadata={
        "section": "Detect IS-junctions",
        "comment": "null = use entire IS element"})
    reference_IS_out_span: int = 100
    max_dist_to_IS: int = 10
    trim_jc_flanking: int = 5

    # --- Coverage ---
    average_method: AverageMethod = field(
        default=AverageMethod.MEDIAN, metadata={"section": "Coverage", "comment": "mean, median, or mode"})

    # --- Filtering ---
    min_amplicon_length: int = field(default=100, metadata={
        "section": "Filtering"})
    max_amplicon_length: int = 1_000_000
    replication_copy_number_threshold: float = 1.5
    deletion_copy_number_threshold: float = 0.3

    # --- Synthetic junction alignment and detection ---
    alignment_filter_params: Optional[AlignmentFilterParams] = field(default=None, metadata={
        "section": "Synthetic junction alignment and detection"})
    alignment_analysis_params: Optional[AlignmentClassifyParams] = None
    bowtie_params: Optional[BowtieParams] = None
    jc_call_params: Optional[JcCallParams] = None

    # --- Misc ---
    threads: int = field(default=4, metadata={
        "section": "Misc",
        "comment": "alignment threads"})
    breseq_output_size_threshold: int = field(default=10_000, metadata={
        "comment": "num lines above which breseq output is considered too large"})
    min_num_bases: int = field(default=80_000_000, metadata={
        "comment": "warn if total fastq sequencing depth is below this threshold"})
    remove_jc_breseq_reject: bool = field(default=False, metadata={
        "comment": "whether to remove JC records labeled as breseq-rejected"})
    create_plots: bool = field(default=True, metadata={
        "comment": "whether to create coverage plots"})

    @property
    def has_ancestor(self) -> bool:
        """True if ancestor comparison should be performed."""
        return self.anc_fastq_path is not None

    @property
    def _anc_folder_name(self) -> str:
        """Folder name for organizing runs (anc_name if exists, else iso_name)."""
        return self.anc_name if self.has_ancestor else self.iso_name

    @property
    def iso_run_dir(self) -> Path:
        """Get the isolate run directory path.

        Structure: {output_dir}/{ref_name}/{anc_folder}/{iso_name}/
        """
        return self.output_dir / self.ref_name / self._anc_folder_name / self.iso_name

    @property
    def anc_run_dir(self) -> Path:
        """Get the ancestor run directory path.

        Structure: {output_dir}/{ref_name}/{anc_folder}/{anc_folder}/
        """
        return self.output_dir / self.ref_name / self._anc_folder_name / self._anc_folder_name

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

    def get_breseq_paths(self) -> tuple[Path, Optional[Path]]:
        """Get isolate and ancestor breseq paths.
        """
        return self.get_iso_breseq_path(), self.get_anc_breseq_path()

    def validate_args(self) -> list[str]:
        """Validate required arguments for this config."""
        errors: list[str] = []

        if self.ref_name is None:
            errors.append("ref_name is required")
        if self.iso_fastq_path is None:
            errors.append("iso_fastq_path is required")

        return errors

    def validate_paths(self) -> list[str]:
        """Validate input paths for this config."""
        errors: list[str] = []

        def check(name: str, path: Optional[Path]) -> None:
            if path is not None and not path.exists():
                errors.append(f"{name} does not exist: {path}")

        check("iso_fastq_path", self.iso_fastq_path)
        check("anc_fastq_path", self.anc_fastq_path)
        check("iso_breseq_path", self.iso_breseq_path)
        check("anc_breseq_path", self.anc_breseq_path)
        if not self.ncbi:
            check("ref_path", self.ref_path)

        return errors

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
            elif is_dataclass(val) and not isinstance(val, type):
                return convert_value({f.name: getattr(val, f.name) for f in fields(val)})
            else:
                return val

        # Convert Config to dict, handling Path objects and tuples
        config_dict = {}
        for key, value in self.__dict__.items():
            config_dict[key] = convert_value(value)

        return config_dict

    def save_to_file(self, config_path: Path, log: bool = True,
                     header: Optional[list[str]] = None) -> None:
        """Save config to specified path with annotated YAML.

        Args:
            config_path: Path to save config file
            log: Whether to log the save location (default: True)
            header: Optional header comment lines for the YAML file
        """
        config_path = Path(config_path)
        save_annotated_yaml(self.to_yaml_dict(), fields(self), config_path, header)

        if log:
            from amplifinder.logger import logger
            logger.info(f"Saved config to:\n{config_path}")

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
        path_fields = ['iso_fastq_path', 'anc_fastq_path', 'output_dir', 'ref_path',
                       'iso_breseq_path', 'anc_breseq_path']
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
        if self.iso_fastq_path is not None:
            self.iso_fastq_path = Path(self.iso_fastq_path).resolve()
        self.output_dir = Path(self.output_dir).resolve()
        self.ref_path = Path(self.ref_path).resolve()

        if self.anc_fastq_path is not None:
            self.anc_fastq_path = Path(self.anc_fastq_path).resolve()
        if self.iso_breseq_path is not None:
            self.iso_breseq_path = Path(self.iso_breseq_path).resolve()
        if self.anc_breseq_path is not None:
            self.anc_breseq_path = Path(self.anc_breseq_path).resolve()

        # Derive names from paths if not provided
        if self.iso_name is None and self.iso_fastq_path is not None:
            self.iso_name = self.iso_fastq_path.stem
        if self.anc_name is None and self.anc_fastq_path is not None:
            self.anc_name = self.anc_fastq_path.stem

        # Normalize params
        def init_params(cls, obj):
            if obj is None:
                return cls()
            if isinstance(obj, cls):
                return obj
            return cls(**obj)
        self.alignment_filter_params = init_params(AlignmentFilterParams, self.alignment_filter_params)
        self.alignment_analysis_params = init_params(AlignmentClassifyParams, self.alignment_analysis_params)
        self.bowtie_params = init_params(BowtieParams, self.bowtie_params)
        self.jc_call_params = init_params(JcCallParams, self.jc_call_params)

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


def _get_config_defaults() -> dict[str, Any]:
    """Extract default values from Config dataclass fields.

    Returns dictionary of defaults for optional Config fields,
    converting Enum defaults to their values.
    """
    defaults = {}
    for f in fields(Config):
        # Skip fields without defaults
        if f.default is MISSING and f.default_factory is MISSING:
            continue

        # Get default value
        if f.default is not MISSING:
            value = f.default
        elif f.default_factory is not MISSING:
            value = f.default_factory()
        else:
            continue

        # Convert Enum to value for serialization
        if isinstance(value, Enum):
            value = value.value

        defaults[f.name] = value

    return defaults


def merge_config(*configs: Optional[dict[str, Any]]) -> dict[str, Any]:
    """Merge configuration dictionaries with later configs taking priority.

    Priority: rightmost config > ... > leftmost config > defaults

    Args:
        *configs: Configuration dictionaries to merge (None values are skipped)

    Returns:
        Merged configuration dictionary
    """
    merged = _get_config_defaults()

    for config in configs:
        if config:
            for key, value in config.items():
                if value is not None:
                    merged[key] = value

    _normalize_alignment_params(merged)
    return merged

# Configuration

AmpliFinder has two independent configuration systems that control different things.

## Environment Config (`amplifinder.yaml`)

Server-specific tool paths and environment settings. Loaded automatically at import time.

**Search locations** (first found wins):
1. User home: `~/.amplifinder/amplifinder.yaml`
2. System-wide: `/etc/amplifinder/amplifinder.yaml`
3. Project root: `amplifinder.yaml` (bundled defaults)

To customize, copy the bundled [`amplifinder.yaml`](../amplifinder.yaml) to `~/.amplifinder/amplifinder.yaml` and edit it there.

**Fields:**

| Field | Description | Default |
|-------|-------------|---------|
| `blast_path` | Directory containing BLAST+ binaries | system PATH |
| `samtools_path` | Directory containing samtools | system PATH |
| `bowtie2_path` | Directory containing bowtie2 | system PATH |
| `breseq_docker` | Run breseq via Docker | `true` |
| `isdb_path` | ISfinder database path | bundled |
| `isfinder_cache` | Cache for bundled ISfinder DB | `~/.amplifinder/ISfinderDB` |
| `isescan_env_name` | Conda env name for ISEScan (null = use current environment) | `isescan` |

**Builtin Configuration:** see [`amplifinder.yaml`](../amplifinder.yaml) in the project root.

## Run Config (`--config`)

Per-run parameters: input files, thresholds, and detection settings.

**Priority:** CLI arguments > `--config` file > built-in defaults.

A complete template with all fields and their built-in defaults is provided in [`examples/base_config.yaml`](../examples/base_config.yaml).

### Using config files

```bash
# Run with a config file (CLI args override file values)
amplifinder -i isolate.fastq -r U00096 --config params.yaml

# Generate a complete config template with all defaults
amplifinder -i isolate.fastq -r U00096 --create-config my_run.yaml

# Rerun from a saved run config
amplifinder --config output/U00096/ancestor/isolate/run_config.yaml
```

Every run saves its full configuration to `run_config.yaml` in the output directory. This file can be passed back to `--config` to reproduce the run.

### Key run parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `iso_fastq_path` | Isolate FASTQ path | Required |
| `ref_name` | NCBI accession | Required |
| `anc_fastq_path` | Ancestor FASTQ path | None |
| `output_dir` | Output directory | `output` |
| `threads` | Alignment threads | 4 |
| `is_detection_method` | IS detection method to USE: `genbank` (default), `isfinder`, or `isescan` | `genbank` |
| `run_comparison_methods` | Additional methods to RUN for comparison (e.g., `[isfinder]`) | `[]` |
| `average_method` | Coverage statistic (mean/median/mode) | median |
| `min_amplicon_length` | Minimum amplicon length (bp) | 100 |
| `max_amplicon_length` | Maximum amplicon length (bp) | 1,000,000 |
| `replication_copy_number_threshold` | Min copy number for amplification | 1.5 |
| `deletion_copy_number_threshold` | Max copy number for deletion | 0.3 |
| `create_plots` | Generate coverage plots | True |

Advanced alignment parameters (`alignment_filter_params`, `alignment_analysis_params`, `bowtie_params`, `jc_call_params`) can be set as nested objects. Use `--create-config` to see all available fields with their defaults.

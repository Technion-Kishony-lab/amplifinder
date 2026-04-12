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

A complete template with all fields and their built-in defaults is provided in [`base_config.yaml`](../base_config.yaml) at the repository root (use `amplifinder --create-config` to regenerate one with your environment).

### Using config files

```bash
# Run with a config file (CLI args override file values)
amplifinder -i path/to/isolate_fastq/ -r U00096 --config params.yaml

# Generate a complete config template with all defaults
amplifinder -i path/to/isolate_fastq/ -r U00096 --create-config my_run.yaml

# Rerun from a saved run config
amplifinder --config output/U00096/ancestor/isolate/run_config.yaml
```

Every run saves its full configuration to `run_config.yaml` in the output directory. This file can be passed back to `--config` to reproduce the run.

**FASTQ inputs:** `iso_fastq_path` and `anc_fastq_path` must be **directories** containing the sample’s `*.fastq*` files (paired-end mates and/or extra lanes in the same folder; single-end can be one file in that folder). All matching files are passed to breseq and used for read-length estimation and downstream alignment.

### Reference genome

AmpliFinder does **not** take a path to a reference FASTA on the command line. The run always uses **`ref_name`** (identifier, usually an NCBI accession) plus **`ref_path`** (a directory) and **`ncbi`** (whether to download from NCBI when the genome is not already cached).

**Default (`ncbi: true`):** If the genome is not yet under `ref_path`, it is fetched from NCBI: GenBank is saved and a FASTA is derived. Cached layout:

```text
ref_path/
├── {ref_name}.json          # optional: maps accession → GenBank locus name
├── genbank/
│   └── {locus}.gb
└── fasta/
    └── {locus}.fasta
```

**`--no-ncbi` / `ncbi: false` (local / offline only):** No download is attempted. The genome **must** already exist under `ref_path` with **at least one** of `genbank/{locus}.gb` or `fasta/{locus}.fasta`. If your files use a different stem than `ref_name`, add `ref_path/{ref_name}.json` with `{"accession": "<ref_name>", "name": "<locus>"}` so the pipeline resolves the correct `locus` name.

Additional requirements when **`ncbi` is false**:

1. **`ref_path` directory name:** The default folder name `genomesDB` is reserved for NCBI-backed references. For purely local references, set `ref_path` to a **different** directory name (or path whose final component is not `genomesDB`).
2. **IS detection:** You **must** set `is_detection_method` to **`isfinder`** or **`isescan`**. The configuration rejects **`genbank`** when `ncbi` is false (even if you provide a local `.gb` file). The default **`genbank`** mode is intended for NCBI-fetched references whose GenBank features are parsed for insertion-sequence features.

Breseq still receives GenBank when present (preferred), otherwise FASTA.

### Key run parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `iso_fastq_path` | Directory of isolate FASTQ files | Required |
| `ref_name` | Reference identifier (usually NCBI accession); not a path to a FASTA | Required |
| `ref_path` | Directory where reference GenBank/FASTA files are stored or downloaded | `genomesDB` |
| `ncbi` | Fetch `ref_name` from NCBI if missing from `ref_path` | `true` |
| `anc_fastq_path` | Directory of ancestor FASTQ files | None |
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

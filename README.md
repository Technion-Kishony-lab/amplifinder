# AmpliFinder

A Python tool for detecting gene amplifications and deletions via Insertion Sequence (IS) junction analysis from whole-genome sequencing data.

## Overview

AmpliFinder analyzes paired-end sequencing data to identify IS-mediated genomic rearrangements by:
- Aligning reads using breseq
- Detecting junction reads spanning IS elements
- Classifying candidate amplicons/deletions based on coverage patterns
- Comparing isolate vs ancestor samples for normalized coverage analysis

## Installation

### Requirements

- Python 3.9+
- External tools:
  - **breseq** (v0.37+) - for read alignment and variant calling
  - **bowtie2** (v2.3.5.1+) - for junction alignment
  - **BLAST+** - for ISfinder database searches (optional)

### Install from Source

```bash
git clone https://github.com/kishonylab/amplifinder.git
cd amplifinder
pip install -e .
```

### Dependencies

Python dependencies are automatically installed:
- numpy, pandas, pydantic
- biopython, pysam
- openpyxl (for Excel export)
- matplotlib (for visualization)
- click, pyyaml

## Usage

### Basic Usage (Isolate Only)

Run analysis on a single isolate without ancestor comparison (raw coverage):

```bash
amplifinder -i isolate.fastq --ref-name U00096 --iso-name my_isolate
```

### With Ancestor Comparison

Run analysis comparing isolate to ancestor (normalized coverage):

```bash
amplifinder -i isolate.fastq --ref-name U00096 --iso-name my_isolate \
            -a ancestor.fastq --anc-name my_ancestor
```

The pipeline will automatically:
1. Run the ancestor analysis first (if not already done)
2. Reuse ancestor results for multiple isolates
3. Normalize isolate coverage by ancestor coverage

### Multiple Isolates with Same Ancestor

```bash
# First isolate
amplifinder -i isolate1.fastq --ref-name U00096 --iso-name iso1 \
            -a ancestor.fastq --anc-name anc

# Second isolate (reuses ancestor run)
amplifinder -i isolate2.fastq --ref-name U00096 --iso-name iso2 \
            -a ancestor.fastq --anc-name anc
```

### Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `-i, --iso-path` | Path to isolate FASTQ file(s) or directory | Required |
| `-r, --ref-name` | Reference genome name (NCBI accession, e.g., U00096) | Required |
| `-a, --anc-path` | Path to ancestor FASTQ file(s) or directory | None |
| `--iso-name` | Isolate name (default: derived from path) | Auto |
| `--anc-name` | Ancestor name (default: derived from path) | Auto |
| `--ref-path` | Path to reference genome cache | `genomesDB` |
| `-o, --output-dir` | Output directory | `output` |
| `--iso-breseq-path` | Path to existing breseq output for isolate | None |
| `--anc-breseq-path` | Path to existing breseq output for ancestor | None |
| `--ncbi/--no-ncbi` | Fetch reference from NCBI | True |
| `--use-isfinder/--no-use-isfinder` | Use ISfinder database for IS detection | False |
| `--config` | Path to YAML/JSON config file | None |
| `--log-level` | Logging level (DEBUG/INFO/WARNING/ERROR) | INFO |
| `--breseq-only` | Only run through breseq step, then exit | False |
| `--visualize` | Visualize results from a completed run directory | None |
| `--save-plots` | Save plots to PNG files instead of displaying | False |

### Visualization

Visualize coverage plots for completed runs:

```bash
# Interactive display
amplifinder --visualize output/U00096/my_ancestor/my_isolate

# Save plots to PNG files
amplifinder --visualize output/U00096/my_ancestor/my_isolate --save-plots
```

### Batch Mode

Run multiple jobs from CSV with concurrent execution:

```bash
amplifinder --batch-input runs.csv --max-parallel 2
```

With base config:

```bash
amplifinder --batch-input runs.csv --config base_config.yaml
```

Outputs status to `run_status.csv` (use `--batch-output` to change). Uses thread pool (add `--use-processes` for real parallelism).

**CSV format:** Each row = one run. Columns = any `Config` field:
- Required: `iso_fastq_path`, `ref_name`
- Optional: `run_id`, `anc_fastq_path`, `iso_name`, `anc_name`, `output_dir`, `iso_breseq_path`, `anc_breseq_path`
- Booleans: `ncbi`, `use_isfinder`, `create_plots` (`true`/`false`, `1`/`0`, `yes`/`no`)

Row values override `--config` base values.

## Output

### Directory Structure

Output is organized by reference → ancestor → isolate:

```
output/
└── {ref_name}/                    # e.g., U00096
    └── {anc_name}/                # ancestor group
        ├── {anc_name}/            # ancestor run
        │   ├── breseq/
        │   ├── tnjc.csv
        │   ├── tnjc2_*.csv
        │   └── jc_*/              # per-candidate directories
        │
        └── {iso_name}/            # isolate run
            ├── breseq/
            ├── tnjc.csv
            ├── tnjc2_*.csv
            ├── ISJC2.csv
            ├── candidate_amplifications.csv
            ├── run_config.yaml
            └── jc_*/              # per-candidate directories
                ├── junctions.fasta
                ├── sorted.bam
                └── coverage_plot.png  # (if --save-plots)
```

### CSV Output Files

**`ISJC2.csv`**: All analyzed candidates with:
- Genomic positions and directions
- IS element identifiers
- Coverage statistics (mean/median/mode copy number)
- Event classification (`raw_event`, `isolate_architecture`)
- Junction coverage patterns

**`candidate_amplifications.csv`**: Filtered high-confidence candidates:
- Copy number ≥ 1.5x (amplifications) or ≤ 0.3x (deletions)
- Amplicon length ≥ 100 bp
- Based on `isolate_architecture` classification

### Coverage Analysis

- **Without ancestor**: Raw coverage depth
- **With ancestor**: Normalized coverage (isolate/ancestor ratio)

## Configuration

### Global Configuration (`amplifinder.yaml`)

Server-specific tool paths and environment settings. Auto-loaded at import time.

**Location:** Project root (or `~/.amplifinder/` or `/etc/amplifinder/`)

```yaml
# amplifinder.yaml
blast_path: /path/to/blast/bin
samtools_path: /path/to/samtools
breseq_docker: true
```

### Run Configuration (`--config`)

Override default run parameters with a custom config file:

```yaml
# params.yaml
iso_path: "path/to/isolate.fastq"
ref_name: "U00096"
anc_path: "path/to/ancestor.fastq"
output_dir: "output"
threads: 8
copy_number_threshold: 1.5
```

**Usage:**
```bash
# With CLI args (CLI overrides config)
amplifinder -i isolate.fastq -r U00096 --config params.yaml

# Config only (all params from file)
amplifinder --config params.yaml

# Rerun using saved run_config.yaml
amplifinder --config output/U00096/ancestor/isolate/run_config.yaml

# Create a complete config file with all defaults (doesn't run)
amplifinder -i isolate.fastq -r U00096 --create-config my_run.yaml
# Then run with: amplifinder --config my_run.yaml
```

## Examples

### Example 1: Single Isolate Analysis

```bash
amplifinder \
    -i data/isolate_R1.fastq,data/isolate_R2.fastq \
    --ref-name U00096 \
    --iso-name E_coli_isolate_1 \
    -o results
```

### Example 2: Isolate vs Ancestor

```bash
amplifinder \
    -i data/isolate.fastq \
    --ref-name U00096 \
    --iso-name evolved_strain \
    -a data/ancestor.fastq \
    --anc-name ancestor_strain \
    -o results
```

### Example 3: Using Existing breseq Output

If you've already run breseq separately:

```bash
amplifinder \
    -i data/isolate.fastq \
    --ref-name U00096 \
    --iso-name my_isolate \
    --iso-breseq-path existing_breseq_output \
    -a data/ancestor.fastq \
    --anc-name my_ancestor \
    --anc-breseq-path existing_anc_breseq_output
```

### Example 4: Create and Use Config Files

```bash
# Create a complete config file from CLI args (shows all defaults)
amplifinder -i data/isolate.fastq -r U00096 -a data/ancestor.fastq \
    --create-config my_run.yaml

# Edit my_run.yaml to customize parameters, then run:
amplifinder --config my_run.yaml

# Rerun from saved run_config.yaml:
amplifinder --config output/U00096/ancestor/isolate/run_config.yaml
```

## Troubleshooting

### Missing External Tools

Ensure breseq, bowtie2, and BLAST+ are installed and in your PATH:

```bash
which breseq
which bowtie2
which blastn
```

### Reference Genome Download Issues

If NCBI download fails, check:
- Internet connectivity
- NCBI accession format (e.g., U00096, CP000000)
- Disk space in `--ref-path` directory

### Memory Issues

For large genomes or high-coverage data:
- Ensure sufficient RAM (8GB+ recommended)
- Consider running breseq separately with `--breseq-only`, then continuing

## Citation

If you use AmpliFinder in your research, please cite:

[Citation information to be added]

## License

MIT License

## Support

For issues and questions, please open an issue on GitHub: https://github.com/kishonylab/amplifinder

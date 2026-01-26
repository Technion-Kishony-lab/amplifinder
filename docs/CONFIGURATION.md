# AmpliFinder Configuration System

## Overview

AmpliFinder uses a three-level configuration system:

```
DEFAULT_CONFIG → amplifinder.yaml → --config params.yaml → CLI args
(hardcoded)      (global env)        (run params)         (highest priority)
```

---

## 1. `amplifinder.yaml` - Global Environment Config

**Purpose:** Server-specific tool paths and environment settings

**Location:** Auto-loaded from:
- Project root: `/path/to/amplifinder/amplifinder.yaml`
- User home: `~/.amplifinder/amplifinder.yaml`
- System-wide: `/etc/amplifinder/amplifinder.yaml`

**Status:** READ at import time by `amplifinder/env.py`

**Example:**
```yaml
# amplifinder.yaml
blast_path: /opt/anaconda3/envs/shared_env/bin
samtools_path: /zdata/user-data/Apps/samtools-1.14
bowtie2_path: null  # use system PATH
breseq_docker: true
isdb_path: null     # use bundled ISfinderDB
```

---

## 2. `--config` - Run Parameter Overrides

**Purpose:** Override default run parameters (optional)

**Usage:** `amplifinder -i isolate.fastq -r U00096 --config params.yaml`

**Status:** READ when provided via CLI

**Example:**
```yaml
# params.yaml
threads: 8
copy_number_threshold: 1.5
output_dir: /custom/output
use_isfinder: true
```

---

## 3. `run_config.yaml` - Saved Run Configuration

**Purpose:** Complete snapshot of config used for a specific run

**Location:** Written to each run directory:
```
{output_dir}/{ref_name}/{anc_name}/{iso_name}/run_config.yaml
```

**Status:** WRITTEN by pipeline after merging all configs

**Rerun with saved config:**
```bash
# Use the saved config to rerun or reproduce results
amplifinder --config output/U00096/ancestor/isolate/run_config.yaml
```

All parameters come from the config file - no CLI args needed!

---

## Configuration Priority

When the same parameter is specified in multiple places:

**Highest → Lowest:**
1. CLI arguments (e.g., `-i isolate.fastq`)
2. `--config` file (e.g., `params.yaml`)
3. `amplifinder.yaml` (global environment)
4. `DEFAULT_CONFIG` (hardcoded defaults in `config.py`)

---

## Creating Config Files

### Generate complete config with `--create-config`

Create a config file showing all parameters with their values (merged defaults + CLI overrides):

```bash
# Create from CLI args
amplifinder -i isolate.fastq -r U00096 -a ancestor.fastq \
    --create-config my_run.yaml

# Creates a file with ALL parameters (including nested alignment params)
# Edit the file, then run with: amplifinder --config my_run.yaml
```

**This is useful for:**
- Seeing all available parameters and their defaults
- Creating a template to customize
- Documenting exactly what parameters were used

---

## Batch Processing

For processing many isolates, use a shared params file:

**`shared_params.yaml`:**
```yaml
threads: 8
ref_name: U00096
output_dir: /project/output
copy_number_threshold: 1.5
```

**`run_batch.sh`:**
```bash
for iso in isolates/*.fastq; do
  amplifinder -i $iso -a ancestor.fastq --config shared_params.yaml
done
```

---

## Key Files Summary

| File | Purpose | Read/Write | Location |
|------|---------|------------|----------|
| `amplifinder.yaml` | Global environment | READ (auto) | Project root |
| `--config params.yaml` | Run overrides | READ (optional) | User-specified |
| `run_config.yaml` | Run snapshot | WRITTEN | Output run dir |

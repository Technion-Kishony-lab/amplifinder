# AmpliFinder

Detect Insertion Sequence (IS)-mediated gene amplifications and deletions from whole-genome sequencing data. AmpliFinder analyzes paired-end reads by aligning them with breseq, detecting junction reads at IS element boundaries, pairing junctions into candidate amplicons, and classifying events based on coverage and read-alignment patterns. An optional ancestor comparison normalizes coverage and compare IS-junctions to distinguish de novo from ancestral events.

## Requirements

- **Python 3.11+**
- **breseq** (read alignment and junction detection)
- **bowtie2** (synthetic junction alignment)
- **BLAST+** (optional, for ISfinder-based IS detection)
- **ISEScan** (optional, for ISEScan-based IS detection; run via separate conda env)

## Installation

```bash
git clone https://github.com/kishonylab/amplifinder.git
cd amplifinder
pip install -e .
```

Python dependencies (numpy, pandas, biopython, pysam, click, pyyaml, matplotlib, rich, etc.) are installed automatically.

## Quick Start

### Isolate only (raw coverage)

```bash
amplifinder -i isolate.fastq -r U00096 --iso-name my_isolate
```

### Isolate vs ancestor (normalized coverage and IS-junction comparison)

```bash
amplifinder -i isolate.fastq -r U00096 --iso-name my_isolate \
            -a ancestor.fastq --anc-name my_ancestor
```

Ancestor results are cached and reused across isolates sharing the same ancestor.

### Batch mode

```bash
amplifinder --batch-input runs.csv --config base_params.yaml --max-parallel 4
```

`runs.csv` should be a CSV files with columns matching any subset of a run configuration fields, see: [`examples/base_config.yaml`](examples/base_config.yaml). Row values override the base config. 

The batch run status is written to `run_status.csv`.

## CLI Reference

| Option | Description | Default | Config field |
|--------|-------------|---------|--------------|
| `-i, --iso-path` | Isolate FASTQ file(s) or directory | Required | `iso_fastq_path` |
| `-r, --ref-name` | Reference genome NCBI accession (e.g. U00096) | Required | `ref_name` |
| `-a, --anc-path` | Ancestor FASTQ file(s) or directory | None | `anc_fastq_path` |
| `--iso-name` | Isolate name | derived from path | `iso_name` |
| `--anc-name` | Ancestor name | derived from path | `anc_name` |
| `--ref-path` | Reference genome cache directory | `genomesDB` | `ref_path` |
| `-o, --output-dir` | Output directory | `output` | `output_dir` |
| `--iso-breseq-path` | Existing breseq output for isolate | None | `iso_breseq_path` |
| `--anc-breseq-path` | Existing breseq output for ancestor | None | `anc_breseq_path` |
| `--ncbi/--no-ncbi` | Fetch reference from NCBI | True | `ncbi` |
| `--is-detection-method` | IS detection method: `genbank` (default), `isfinder`, or `isescan` | `genbank` | `is_detection_method` |
| `--config` | YAML/JSON config file for run parameters | None | |
| `--create-config` | Save merged config to file and exit | None | |
| `--create-plots/--no-create-plots` | Generate coverage plots | True | `create_plots` |
| `--breseq-only` | Run only breseq, then exit | False | |
| `--verbose` | Report step progress | False | |
| `--debug` | Enable debug logging | False | |
| `--batch-input` | Batch CSV input file | None | |
| `--batch-output` | Batch status CSV output | `run_status.csv` | |
| `--max-parallel` | Max concurrent batch runs | 2 | |
| `--use-processes` | Use process pool (vs thread pool) | False | |
| `--version` | Show version | | |

## Output

```
output/
└── {ref_name}/
    └── {anc_name}/
        ├── {anc_name}/                  # ancestor (breseq + junction alignments)
        │   └── junctions/
        │       └── jc_.../
        │           ├── junctions.fasta
        │           ├── sorted.bam
        │           └── jc_read_counts.csv
        └── {iso_name}/                  # isolate run folder
            ├── run_config.yaml          # full config snapshot (rerunnable)
            ├── run_log.txt              # pipeline log
            ├── warnings.txt             # warnings only
            ├── summary.yaml             # final results (yaml)
            ├── summary.json             # final results (json)
            ├── genbank_ref_tnjc.csv     # reference IS junctions
            ├── tnjc.csv                 # junctions matched to IS elements
            ├── tnjc2_A_raw.csv          # pairs of IS-junctions matching the
                                           two sides of same IS
            ├── tnjc2_B_linked.csv       # linking each IS-pair to a single-locus
                                           IS pairs with a matching IS-junction
            ├── tnjc2_C_covered.csv      # adding amplicon coverage
            ├── tnjc2_D_syn_jcts.csv     # adding synthetic junction paths
            ├── tnjc2_E_analyzed.csv     # adding junction read counts
            ├── tnjc2_F_classified.csv   # adding classification
            └── junctions/
                └── jc_.../
                    ├── junctions.fasta
                    ├── sorted.bam
                    ├── jc_read_counts.csv
                    ├── jct_coverages.png
                    └── amplicon_coverage.png
```

**`summary.yaml`** / **`summary.json`** — classified amplification/deletion candidates with genomic positions, IS identifiers, copy number, and event classification.

**`tnjc2_F_classified.csv`** — full table of all analyzed junction pairs across all pipeline stages.

## Configuration

AmpliFinder uses two independent configuration systems:

- **`amplifinder.yaml`** — server-specific tool paths (blast, samtools, bowtie2, breseq docker). Loaded automatically at import time.
- **`--config`** — per-run parameters (input files, thresholds, detection settings). CLI args take priority over the config file.

See [docs/configuration.md](docs/configuration.md) for details.

## Pipeline

The pipeline runs 16 steps from FASTQ to classified amplification candidates. See [docs/pipeline.md](docs/pipeline.md) for the full flow diagram and step descriptions.

## Credit

*AmpliFinder* was developed by Idan Yelin and Roy Kishony.

## License

MIT License. See [LICENSE](LICENSE) for details.

## Support

Issues and questions: https://github.com/kishonylab/amplifinder/issues

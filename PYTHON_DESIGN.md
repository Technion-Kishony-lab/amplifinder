# AmpliFinder Python Translation Design

## Overview

AmpliFinder is a bioinformatics pipeline for detecting **IS-mediated gene amplifications and deletions** from whole-genome sequencing data. It analyzes paired-end reads to identify junction sequences at insertion sequence (IS) element boundaries.

---

## Package Structure

```
amplifinder/
├── __init__.py
├── __main__.py              # CLI entry point
├── cli.py                   # Argument parsing (click)
├── config.py                # Configuration management
├── logger.py                # Logging utilities
│
├── steps/                   # Pipeline step classes
│   ├── __init__.py
│   ├── base.py              # Step base class (caching, I/O tracking)
│   ├── initialize.py        # InitStep - directory setup
│   ├── get_reference.py     # GetReferenceStep - fetch genome
│   ├── isfinder.py          # ISfinderStep - BLAST vs ISfinder DB
│   └── run_breseq.py        # RunBreseqStep - breseq execution
│
├── tools/                   # External tool runners + parsers
│   ├── __init__.py
│   ├── blast.py             # run_blastn + parse_blast_csv
│   └── breseq.py            # run_breseq (Docker) + output parsing
│
├── utils/                   # General utilities
│   ├── __init__.py
│   └── fasta.py             # FASTA utilities (read_fasta_lengths)
│
├── data_types/              # Data structures
│   ├── __init__.py
│   ├── genome.py            # Genome dataclass
│   ├── junction.py          # Junction dataclass
│   ├── schema.py            # DataFrame schemas
│   └── tabularable.py       # Tabularable mixin
│
├── data/                    # Bundled data
│   ├── __init__.py          # Data loaders
│   ├── fields/              # breseq field definitions
│   │   ├── JC_fields.csv
│   │   ├── SNP_fields.csv
│   │   └── ...
│   └── ISfinderDB/
│       └── IS.fna
│
└── tests/                   # Test suite
    ├── __init__.py
    └── fixtures/
```

---

## Key Modules

### 1. `steps/base.py` — Step Base Class
```python
class Step(ABC):
    """Base class for pipeline steps with input/output file tracking.
    
    Handles:
    - Skip if all outputs exist (unless force=True)
    - Clean partial outputs before re-run
    - Global and step-specific force control
    """
    global_force: bool = False
    
    def __init__(self, inputs=None, outputs=None, force=None): ...
    def run(self) -> bool: ...  # Returns True if ran, False if skipped
    @abstractmethod
    def _run(self) -> None: ...  # Override in subclass
```

### 2. `data_types/genome.py` — Genome Data
```python
@dataclass
class Genome:
    name: str
    fasta_path: Optional[Path]
    genbank_path: Optional[Path]
    length: int
    scaffolds: List[str]
```

### 3. `data_types/junction.py` — Junction Data (from breseq)
```python
@dataclass
class Junction(Tabularable):
    """breseq junction record"""
    num: int
    side_1_seq_id: str
    side_1_position: int
    side_1_strand: int
    side_2_seq_id: str
    side_2_position: int
    side_2_strand: int
    ...
```

### 4. `data_types/schema.py` — DataFrame Schemas
```python
class Schema:
    """DataFrame schema with column names and dtypes."""
    columns: Dict[str, type]
    def empty(self) -> pd.DataFrame: ...
    def validate(self, df: pd.DataFrame) -> bool: ...
```

### 5. IS Element Classes (TBD)
```python
@dataclass
class ISElement:
    """Base IS element with ID and orientation"""
    id: int
    orientation: int  # -1, 0, +1
    
    @property
    def name(self) -> str:
        postfix = {-1: "R", 0: "u", 1: "F"}
        return f"{self.id}{postfix[self.orientation]}"
```

### 6. Junction Pairs (TBD) — `ISJC2`
```python
@dataclass
class ISJC2:
    """Pair of IS-associated junctions (candidate amplification)"""
    jc_num: Tuple[int, int]
    pos_chr: Tuple[int, int]
    pos_IS: Tuple[int, int]
    dir_chr: Tuple[int, int]
    dir_IS: Tuple[int, int]
    ref_is: Tuple[int, int]
    scaf_chr: int
    IS_elements: List[DirectedIS]
    span_origin: bool
    amplicon_length: int
    amplicon_coverage: float
    amplicon_coverage_mode: float
    event: str  # "flanked", "unflanked", "transposition", etc.
```

---

## Key Functions

### Pipeline Flow

| Step | MATLAB Function | Python Module |
|------|-----------------|---------------|
| 1. Parse args | `parse_AmpliFinder_arguments` | `cli.main()` |
| 2. Init dirs | - | `steps.InitStep` |
| 3. Get reference | `curate_reference` | `steps.GetReferenceStep` |
| 4. Find IS in ref | `findISinRef` | `steps.ISfinderStep` |
| 5. Run breseq | `run_breseq` | `steps.RunBreseqStep` → `tools.breseq.run_breseq()` |
| 6. Parse breseq | `breseq2mat` | `tools.breseq.parse_breseq_output()` |
| 7. Assign ISs | `assign_potential_ISs` | TBD |
| 8. Combine pairs | `combine_ISJC_pairs` | TBD |
| 9. Calc coverage | `calc_coverage_ISJC2` | `tools.breseq.parse_coverage()` |
| 10. Classify | `classify_ISJC2` | TBD |
| 11. Export | `export_ISJC2` | TBD |

---

## Main Data Structures

### 1. `Isolate` Table (DataFrame)
```
| Column         | Type   | Description                    |
|----------------|--------|--------------------------------|
| iso            | str    | Isolate name                   |
| anc            | str    | Ancestor name                  |
| fastq_path     | Path   | Path to FASTQ files            |
| breseq_path    | Path   | Path to breseq output          |
| iso_outpath    | Path   | Output directory               |
| ref            | str    | Reference genome name(s)       |
| ref_path       | Path   | Path to reference files        |
| use_ISfinder   | bool   | Use ISfinder DB                |
| read_length    | int    | Mean read length               |
```

### 2. `JC` Table (Junction Candidates from breseq)
```
| Column              | Type   | Description                |
|---------------------|--------|----------------------------|
| num                 | int    | Junction ID                |
| jscaf1, jscaf2      | int    | Scaffold indices           |
| pos1, pos2          | int    | Positions                  |
| dir1, dir2          | int    | Directions (-1/+1)         |
| flanking_left/right | int    | Flanking sequence lengths  |
| reject              | str    | breseq rejection reason    |
| refIS               | int    | Reference IS index (0=none)|
| IS_side             | list   | Matching IS sides          |
| IS_dis              | list   | Distances to IS            |
```

### 3. `ISJC2` Table (Junction Pairs)
```
| Column                  | Type   | Description                    |
|-------------------------|--------|--------------------------------|
| JC_num                  | (2,)   | Junction IDs [left, right]     |
| pos_Chr                 | (2,)   | Chromosome positions           |
| dir_Chr                 | (2,)   | Chromosome directions          |
| pos_IS                  | (2,)   | IS positions                   |
| dir_IS                  | (2,)   | IS directions                  |
| refIS                   | (2,)   | Reference IS indices           |
| scaf_Chr                | int    | Chromosome scaffold            |
| IS                      | list   | DirectedIS objects             |
| span_origin             | bool   | Spans replication origin       |
| amplicon_length         | int    | Amplified region length        |
| amplicon_coverage       | float  | Median copy number             |
| amplicon_coverage_mode  | float  | Mode copy number               |
| event                   | str    | Classification label           |
| chosen_IS               | IS     | Selected IS for analysis       |
```

### 4. `IS_loc` Table (IS Element Locations)
```
| Column      | Type   | Description              |
|-------------|--------|--------------------------|
| ID          | int    | IS identifier            |
| IS_Name     | str    | IS element name          |
| IS_scaf     | str    | Scaffold name            |
| IS_jscaf    | int    | Scaffold index           |
| LocLeft     | int    | Left position            |
| LocRight    | int    | Right position           |
| Complement  | bool   | On complement strand     |
```

---

## Configuration (`config.txt`)

```ini
# External tools
BRESEQ_DOCKER = 1
BLASTN_PATH = /path/to/blastn
SAMTOOLS_PATH = /path/to/samtools

# IS detection
MAX_DIST_TO_IS = 10
LENGTH_SEQ_INTO_IS = 200
REFERENCE_IS_OUT_SPAN = 100

# Junction filtering
MIN_JCT_COV = 5
MIN_AMPLICON_LENGTH = 30
MAX_AMPLICON_LENGTH = 1000000
FILTER_AMPLICON_LENGTH = 100

# Copy number thresholds
COPY_NMBR_THRS = 1.5        # Amplification threshold
DEL_COPY_NMBR_THRS = 0.3    # Deletion threshold

# Alignment parameters
SCORE_MIN1 = 0
SCORE_MIN2 = -0.1
MISMATCH_PENALTY = 5,5
REQ_OVERLAP = 12
```

---

## Installation & Distribution

### Package Setup (`pyproject.toml`)
```toml
[build-system]
requires = ["setuptools>=61.0", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "amplifinder"
version = "1.0.0"
description = "Detect IS-mediated gene amplifications from WGS data"
readme = "README.md"
requires-python = ">=3.9"
license = {text = "MIT"}
authors = [{name = "Your Name", email = "you@example.com"}]
keywords = ["bioinformatics", "genomics", "amplification", "IS-elements"]
classifiers = [
    "Development Status :: 4 - Beta",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
    "License :: OSI Approved :: MIT License",
    "Programming Language :: Python :: 3.9",
    "Programming Language :: Python :: 3.10",
    "Programming Language :: Python :: 3.11",
]

dependencies = [
    "numpy>=1.21",
    "pandas>=1.3",
    "biopython>=1.79",
    "pysam>=0.19",
    "openpyxl>=3.0",
    "click>=8.0",
]

[project.optional-dependencies]
dev = ["pytest", "pytest-cov", "black", "ruff", "mypy"]

[project.scripts]
amplifinder = "amplifinder.cli:main"

[project.urls]
Homepage = "https://github.com/yourusername/amplifinder"
Documentation = "https://amplifinder.readthedocs.io"
```

### Requirements
```
# requirements.txt
numpy>=1.21
pandas>=1.3
biopython>=1.79
pysam>=0.19
openpyxl>=3.0
click>=8.0
```

### External Dependencies
- **breseq** (Docker or local)
- **bowtie2**
- **samtools**
- **BLAST+**

### Bioconda Distribution (Recommended)

Bioconda is the standard distribution channel for bioinformatics tools. It handles external dependencies (bowtie2, samtools, blast) automatically.

```bash
# Install via conda (once recipe is published)
conda install -c bioconda amplifinder
```

**Bioconda recipe** (`recipes/amplifinder/meta.yaml`):
```yaml
package:
  name: amplifinder
  version: "1.0.0"

source:
  url: https://github.com/yourusername/amplifinder/archive/v1.0.0.tar.gz

build:
  number: 0
  script: python -m pip install . --no-deps -vv

requirements:
  host:
    - python >=3.9
    - pip
  run:
    - python >=3.9
    - numpy >=1.21
    - pandas >=1.3
    - biopython >=1.79
    - pysam >=0.19
    - openpyxl >=3.0
    - click >=8.0
    - bowtie2
    - samtools
    - blast

test:
  commands:
    - amplifinder --help

about:
  home: https://github.com/yourusername/amplifinder
  license: MIT
  summary: Detect IS-mediated gene amplifications from WGS data
```

---

## Usage Examples

### Command Line
```bash
# Basic usage
amplifinder -i /path/to/isolate.fastq -a /path/to/ancestor.fastq -r U00096

# With ISfinder
amplifinder -i isolate/ -a ancestor/ -r U00096 --use-isfinder

# Local reference (non-NCBI)
amplifinder -i isolate/ -r mygenome --ref-path /genomes/ --ncbi false --use-isfinder

# Using config file (recommended for reproducibility)
amplifinder --config params.yaml
amplifinder --config params.json
```

### Config File Format

Supports **YAML** (recommended) or **JSON** for reproducible analyses.

**Example `params.yaml`:**
```yaml
# AmpliFinder configuration
# -------------------------

# Input/Output
isolate_path: /data/isolate/
ancestor_path: /data/ancestor/
output_dir: ./results/
ref_name: U00096

# Reference options
ncbi: true
use_ISfinder: false
ref_path: genomesDB/

# IS detection parameters
max_dist_to_IS: 10
length_seq_into_is: 200
reference_IS_out_span: 100

# Junction filtering
min_jct_cov: 5
min_amplicon_length: 30
max_amplicon_length: 1000000
filter_amplicon_length: 100

# Copy number thresholds
copy_number_threshold: 1.5      # amplification
del_copy_number_threshold: 0.3  # deletion

# Alignment parameters
score_min: [0, -0.1]
mismatch_penalty: [5, 5]
req_overlap: 12

# External tools (optional, auto-detected if in PATH)
# breseq_docker: true
# samtools_path: /usr/bin/samtools
# blastn_path: /usr/bin/blastn
```

**Equivalent `params.json`:**
```json
{
  "isolate_path": "/data/isolate/",
  "ancestor_path": "/data/ancestor/",
  "output_dir": "./results/",
  "ref_name": "U00096",
  "ncbi": true,
  "use_ISfinder": false,
  "copy_number_threshold": 1.5,
  "del_copy_number_threshold": 0.3,
  "min_amplicon_length": 30
}
```

**Priority:** CLI args > config file > defaults

### Python API
```python
from amplifinder import AmpliFinder

# Run pipeline
af = AmpliFinder(
    iso_path="path/to/isolate.fastq",
    anc_path="path/to/ancestor.fastq", 
    ref_name="U00096"
)
results = af.run()

# Access results
print(results.isjc2)  # DataFrame of candidate amplifications
results.to_excel("output.xlsx")
```

---

## Event Classification Legend

```
>>> IS element
=== Amplified region (cassette)
~~~ Flanking chromosome
--  IS-associated junction (ISJC)

TRANSPOSITION:
~~~AB~~~                    reference
~~~A>>>B~~~                 isolate

UNFLANKED:
~~~Aa======bB~~~            reference  
~~~Aa======b>>>a======bB~~~ isolate (internal IS only)

HEMI-FLANKED:
~~~Aa======bB~~~            reference
~~~A>>>a======b>>>B======~~ isolate (one flanking IS)

FLANKED:
~~~Aa======b>>>B~~~         reference
~~~A>>>a======b>>>B~~~      isolate (both flanking ISs)
```

---

## Key Python Libraries Mapping

| MATLAB Function | Python Equivalent |
|-----------------|-------------------|
| `fastaread/fastawrite` | `Bio.SeqIO` |
| `table` | `pandas.DataFrame` |
| `struct` | `dataclass` / `dict` |
| `save/load .mat` | `pickle` / `pandas.to_parquet` |
| `seqrcomplement` | `Bio.Seq.reverse_complement()` |
| `BioMap` | `pysam.AlignmentFile` |
| `readtable` | `pandas.read_csv/read_excel` |
| `writetable` | `pandas.to_excel` |
| `system` | `subprocess.run` |
| `fullfile` | `pathlib.Path` |



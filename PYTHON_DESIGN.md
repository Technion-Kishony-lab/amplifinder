# AmpliFinder Python Translation Design

## Overview

AmpliFinder is a bioinformatics pipeline for detecting **IS-mediated gene amplifications and deletions** from whole-genome sequencing data. It analyzes paired-end reads to identify junction sequences at insertion sequence (IS) element boundaries.

---

## Package Structure

```
amplifinder/
├── __init__.py
├── __main__.py              # CLI entry point
├── cli.py                   # Argument parsing (argparse)
├── config.py                # Configuration management
├── logger.py                # Logging utilities
│
├── core/                    # Core pipeline logic
│   ├── __init__.py
│   ├── pipeline.py          # Main AmpliFinder orchestration
│   ├── reference.py         # Reference genome handling
│   ├── breseq.py            # breseq wrapper & parsing
│   ├── alignment.py         # bowtie2/samtools wrappers
│   └── coverage.py          # Coverage calculation
│
├── analysis/                # Analysis modules
│   ├── __init__.py
│   ├── junctions.py         # Junction (JC) processing
│   ├── is_elements.py       # IS element detection & assignment
│   ├── isjc.py              # IS-Junction pairs (ISJC)
│   ├── isjc2.py             # Junction pair combinations (ISJC2)
│   ├── classification.py    # Event classification
│   └── amplicons.py         # Amplicon candidate curation
│
├── io/                      # Input/Output
│   ├── __init__.py
│   ├── fastq.py             # FASTQ handling
│   ├── fasta.py             # FASTA read/write
│   ├── bam.py               # BAM parsing
│   ├── genbank.py           # GenBank parsing
│   └── export.py            # Excel/CSV export
│
├── models/                  # Data structures
│   ├── __init__.py
│   ├── isolate.py           # Isolate dataclass
│   ├── junction.py          # JC dataclass
│   ├── is_element.py        # IS element classes
│   └── reference.py         # Reference genome properties
│
├── external/                # External tool wrappers
│   ├── __init__.py
│   ├── breseq.py            # breseq runner (Docker)
│   ├── bowtie2.py           # bowtie2 runner
│   ├── samtools.py          # samtools runner
│   └── blast.py             # BLAST runner
│
├── data/                    # Bundled data
│   └── ISfinderDB/
│       └── IS.fna
│
└── tests/                   # Test suite
    ├── __init__.py
    ├── test_pipeline.py
    └── fixtures/
```

---

## Key Modules

### 1. `core/pipeline.py` — Main Orchestrator
```python
class AmpliFinder:
    """Main pipeline controller"""
    def __init__(self, iso_path, anc_path, ref_name, **kwargs): ...
    def run(self) -> pd.DataFrame: ...
```

### 2. `models/isolate.py` — Isolate Data
```python
@dataclass
class Isolate:
    name: str
    fastq_path: Path
    breseq_path: Path
    output_path: Path
    ref_name: str
    ref_path: Path
    isfinder: bool = False
    read_length: int = None
    ancestor: Optional["Isolate"] = None
```

### 3. `models/is_element.py` — IS Element Classes
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

class DirectedIS(ISElement):
    """IS with direction (Forward/Reverse/Unknown)"""
    POSTFIX = {-1: "R", 0: "u", 1: "F"}

class ISSide(ISElement):
    """IS side (Left/Both/Right)"""
    POSTFIX = {-1: "L", 0: "B", 1: "R"}
```

### 4. `models/junction.py` — Junction Data
```python
@dataclass
class Junction:
    """breseq junction record"""
    num: int
    jscaf1: int
    pos1: int
    dir1: int  # -1 or 1
    jscaf2: int
    pos2: int
    dir2: int
    flanking_left: int
    flanking_right: int
    reject: Optional[str] = None
    # ... other breseq fields
```

### 5. `analysis/isjc2.py` — Junction Pairs
```python
@dataclass
class ISJC2:
    """Pair of IS-associated junctions (candidate amplification)"""
    jc_num: Tuple[int, int]
    pos_chr: Tuple[int, int]
    pos_is: Tuple[int, int]
    dir_chr: Tuple[int, int]
    dir_is: Tuple[int, int]
    ref_is: Tuple[int, int]
    scaf_chr: int
    is_elements: List[DirectedIS]
    span_origin: bool
    amplicon_length: int
    amplicon_coverage: float
    amplicon_coverage_mode: float
    event: str  # "flanked", "unflanked", "transposition", etc.
```

---

## Key Functions

### Pipeline Flow

| Step | MATLAB Function | Python Function |
|------|-----------------|-----------------|
| 1. Parse args | `parse_AmpliFinder_arguments` | `cli.parse_args()` |
| 2. Analyze FASTQ | `analyze_fastq_files` | `io.fastq.analyze()` |
| 3. Get reference | `curate_reference` | `core.reference.curate()` |
| 4. Run breseq | `run_breseq` | `external.breseq.run()` |
| 5. Parse breseq | `breseq2mat` | `core.breseq.parse_output()` |
| 6. Find IS in ref | `findISinRef` | `analysis.is_elements.find_in_ref()` |
| 7. Assign ISs | `assign_potential_ISs` | `analysis.is_elements.assign_to_junctions()` |
| 8. Combine pairs | `combine_ISJC_pairs` | `analysis.isjc2.combine_pairs()` |
| 9. Calc coverage | `calc_coverage_ISJC2` | `core.coverage.calculate()` |
| 10. Classify | `classify_ISJC2` | `analysis.classification.classify_isjc2()` |
| 11. Write FASTA | `junction2fasta` | `io.fasta.write_junctions()` |
| 12. Align | `bowtie2_alignment` | `external.bowtie2.align()` |
| 13. BAM analysis | `bamAnalysis` | `io.bam.analyze()` |
| 14. Final classify | `classify_candidates` | `analysis.classification.classify_candidates()` |
| 15. Export | `export_ISJC2` | `io.export.to_excel()` |

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
| isfinder       | bool   | Use ISfinder DB                |
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
amplifinder -i isolate/ -a ancestor/ -r U00096 --isfinder

# Local reference (non-NCBI)
amplifinder -i isolate/ -r mygenome --ref-path /genomes/ --ncbi false --isfinder

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
isfinder: false
ref_path: genomesDB/

# IS detection parameters
max_dist_to_is: 10
length_seq_into_is: 200
reference_is_out_span: 100

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
  "isfinder": false,
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



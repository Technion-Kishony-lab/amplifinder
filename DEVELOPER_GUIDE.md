# AmpliFinder Developer Documentation

## Architecture Overview

AmpliFinder is a Python pipeline for detecting Insertion Sequence (IS)-mediated gene amplifications and deletions from whole-genome sequencing data. The pipeline processes paired-end sequencing reads through a series of steps that identify, classify, and analyze candidate amplicons.

### Pipeline Flow

```
{xxx}  input
xxx[]  array
- - -  optional
(?)    optional

                                      INPUTS
   ┌────────────────────────────────────┬────────────────────────────────────┐
{FASTQ}                             {ref_name}                          {anc_path}
   │                                    │                               (optional)
   │                                    ▼                                    │
   │                          ┌───────────────────┐                          │
   │                          │ 1. GetRefGenome   │                          │
   │                          └─────────┬─────────┘                          │
   │                                    ▼                                    │
   │           ┌──────────────────    Genome    ──────────────┐              │
   │           │                  (FASTA + GBK)               │              │
   │           │                        │                     │              │
   │           │                        ▼                     │              │
   │           │          ┌────────────────────────────┐      │              │
   │           │          │ 2a. LocateTNsUsingGenbank  │      │              │
   │           │          │ 2b. LocateTNsUsingISfinder │      │              │
   │           │          └─────────────┬──────────────┘      │              │
   │           │                        ▼                     │              │
   │           │                     TnLoc[]                  │              │
   │           │                        │                     │              │
   │           │                        ▼                     │              │
   │           │           ┌────────────────────────┐         │              │
   │           │           │    3a.CreateRefTnJc    │         │              │
   │           │           └────────────┬───────────┘         │              │
   │           │                        ▼                     │              │
   │           │             ┌── RefTnJunction[]              │              │
   │           │             │          │                     │              │
   │           ▼             │          ▼                     ▼              │
   │     ┌───────────┐       │       ┌───────────────────────────┐           │
   │────►│ 4. Breseq │       │       │ 3b. CreateRefTnEndSeqs    │           │
   │     └─────┬─────┘       │       └──────────────┬────────────┘           │
   │           ▼             │                      ▼                        │
   │        breseq JC        │                  SeqRefTnSide[]               │
   │      + coverage         ▼                      │                        │
   │            └─────┬──────┘                      ▼                        │
   │                  ▼            ┌───────────────────┐                     │
   │               all Jc ────────►│  5. CreateTnJc    │                     │
   │                               └─────────┬─────────┘                     │
   │                                         ▼                               │
   │                                    TnJunction[]                         │
   │                                         │                               │
   │                                         ▼                               │
   │                               ┌───────────────────┐                     │
   │                               │  6. CreateTnJc2   │                     │
   │                               └─────────┬─────────┘                     │
   │                                         ▼                               │
   │                                      TnJc2[]                            │
   │                                         │                               │
   │                                         ▼                               │
   │                           ┌─────────────────────────┐                   │
   │                           │ 7. CalcAmpliconCoverage │◄ - - - - - - - - -┤
   │                           └────────────┬────────────┘  anc              │ 
   │                                        ▼               coverage         │
   │                                   CoveredTnJc2[]                        │
   │                                        │                                │
   │                                        ▼                                │
   │                             ┌─────────────────────┐                     │
   │                             │ 8. ClassifyStructure│                     │
   │                             └──────────┬──────────┘                     │
   │                                        ▼                                │
   │                                 ClassifiedTnJc2[] (`RawEvent`)          │
   │                                        │                                │
   │                                        ▼                                │
   │                             ┌─────────────────────┐                     │
   │                             │ 9. FilterCandidates │                     │
   │                             └──────────┬──────────┘                     │
   │                                        ▼                                │
   │                                    CandidateTnJc2[]                     │
   │                                        │                                │
   │                                        ▼                                │
   │                         ┌──────────────────────────────┐                │
   │                         │ 10. CreateSyntheticJunctions │                │
   │                         └──────────────┬───────────────┘                │
   │                                        ▼                                │
   │                                 junctions.fasta                         │
   │                          (7 `JunctionType` per candidate)               │
   │                                        │                                │
   │                                        ▼                                │
   │                     ┌─────────────────────────────────────┐             │
   └────────────────────►│ 11. AlignReadsToSyntheticJunctions  │◄- - - - - - ┘
            {FASTQ}      └──────────────────┬──────────────────┘ {anc_FASTQ}
                                            ▼             
                                         iso.bam           
                                         anc.bam(?)          
                                            │             
                                            ▼             
                                ┌───────────────────────┐   
                                │ 12. AnalyzeAlignments │
                                └────────────┬──────────┘
                                             ▼
                                  AnalyzedCandidateTnJc2[] (left/right/spanning, `JunctionCoverage`)
                                  c_cov, anc_jc_cov(?)
                                             │
                                             ▼
                                ┌────────────────────────┐
                                │ 13. ClassifyCandidates │
                                └────────────┬───────────┘
                                             ▼
                                    AnalyzedCandidateTnJc2[] (RawEvent + EventModifiers)
                                             │
                                             ▼
                                   ┌───────────────────┐
                                   │ 14. Export        │
                                   └─────────┬─────────┘
                                             ▼
                                        ISJC2.csv
                                        candidate_amplifications.csv
                                             │
                                             ▼
                                ┌──────────────────────────┐
                                │ 15. VisualizeCandidates  │◄─── interactive
                                └────────────┬─────────────┘
                                             ▼
                                   coverage plots (matplotlib)
```

## Data Structures

### Base Abstractions

**`Step[T]`** (`steps/base.py`): Abstract pipeline step with:
- Input/output file tracking
- Caching: skip if outputs exist (unless `force=True`)
- Methods: `_calculate_output() → T`, `_save_output()`, `load_outputs()`
- To implement steps without file outputs (memory only), set `output_files=None`

**`Record`** (`data_types/records.py`): Pydantic model with:
- Auto-generated `schema()` from fields
- `from_other()` to convert between record types (copy shared fields)
- Extra fields support

**`RecordTypedDF[T]`** (`data_types/typed_df.py`): Typed DataFrame wrapper:
- Wraps pandas DataFrame with type information
- `from_records()`, `from_csv()`, `to_csv()` methods
- Type-safe access to records

### Record Type Hierarchy

**TN Location Records:**
```python
TnLoc(Record)
├── ID: int
├── TN_Name: str
├── TN_scaf: str
├── LocLeft: int
├── LocRight: int
├── Complement: bool
└── Join: bool

RefTnSide(Record)
├── tn_id: int
├── side: Side
└── distance: Optional[int]

SeqRefTnSide(RefTnSide)
├── seq_fwd: str
└── seq_rc: str

BlastHit(Record)
├── query: str
├── subject: str
├── percent_identical: float
├── length: int
├── mismatch: int
├── gapopen: int
├── qstart: int
├── qend: int
├── sstart: int
├── send: int
├── evalue: float
└── bitscore: float
```

**Junction Records:**
```python
Junction(Record)  # Base junction
├── num: int
├── scaf1: str
├── pos1: int
├── dir1: Orientation
├── scaf2: str
├── pos2: int
├── dir2: Orientation
├── flanking_left: int
└── flanking_right: int

RefTnJunction(Junction)  # Reference TN junction
└── ref_tn_side: RefTnSide

TnJunction(Junction)  # Junction matched to TN element(s)
├── ref_tn_sides: List[RefTnSide]
└── switched: bool
```

**Amplicon Candidate Records (Pipeline Flow):**
```python
TnJc2(Record)  # Base paired junction record
├── jc_num_L: int
├── jc_num_R: int
├── scaf_chr: str
├── pos_chr_L: int
├── pos_chr_R: int
├── pos_tn_L: int
├── pos_tn_R: int
├── dir_chr_L: Orientation
├── dir_chr_R: Orientation
├── dir_tn_L: Orientation
├── dir_tn_R: Orientation
├── tn_ids: List[int]
├── tn_orientations: List[Orientation]
├── span_origin: bool
├── amplicon_length: int
└── complementary_length: int

CoveredTnJc2(TnJc2)  # Step 7: Coverage added
├── ref_name: str
├── iso_name: str
├── anc_name: Optional[str]
├── amplicon_coverage: float
├── genome_coverage: float
├── copy_number: float
├── amplicon_coverage_mode: float
└── copy_number_ratio: Optional[float]

ClassifiedTnJc2(CoveredTnJc2)  # Step 8: Structure classified
├── raw_event: RawEvent
├── shared_tn_ids: List[int]
└── chosen_tn_id: Optional[int]

CandidateTnJc2(ClassifiedTnJc2)  # Step 9: Filtered candidates
└── analysis_dir: str

AnalyzedTnJc2(CandidateTnJc2)  # Step 12: Junction coverage analyzed
├── jc_cov_left: List[int]
├── jc_cov_right: List[int]
├── jc_cov_spanning: List[int]
├── anc_jc_cov_left: Optional[List[int]]
├── anc_jc_cov_right: Optional[List[int]]
├── anc_jc_cov_spanning: Optional[List[int]]
├── isolate_architecture: RawEvent
├── ancestor_architecture: Optional[RawEvent]
├── event: str
└── event_modifiers: List[EventModifier]
```

### Helper Types

```python
class Coverage(NamedTuple):
    """Coverage statistics for a genomic region."""
    mean: float
    median: float
    mode: float

class JunctionCoverage(NamedTuple):
    """Read coverage at a synthetic junction."""
    spanning: int  # reads crossing junction
    left: int      # reads ending at junction
    right: int     # reads starting at junction

```

### Key Enums

```python
class Side(int, Enum):
    """Side of a TN element (left or right)."""
    LEFT = -1
    RIGHT = 1

class Orientation(int, Enum):
    """Orientation relative to reference."""
    FORWARD = 1
    REVERSE = -1
    BOTH = 0

class RawEvent(str, Enum):
    """Structural classification based on junction pair relationships."""
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"
    UNFLANKED = "unflanked"
    HEMI_FLANKED_LEFT = "hemi-flanked left"
    HEMI_FLANKED_RIGHT = "hemi-flanked right"
    FLANKED = "flanked"
    MULTIPLE_SINGLE_LOCUS = "multiple single locus"
    UNRESOLVED = "unresolved"

class JunctionType(int, Enum):
    """The 7 synthetic junction types."""
    LEFT_REF = 1
    LEFT_IS_TRANS = 2
    LEFT_MID_IS = 3
    LOST_IS = 4
    RIGHT_MID_IS = 5
    RIGHT_IS_TRANS = 6
    RIGHT_REF = 7

class EventModifier(str, Enum):
    """Modifiers for classified events."""
    ANCESTRAL = "ancestral"
    DE_NOVO = "de novo"
    LOW_COVERAGE = "low coverage near junction"
```

## Folder Structure

### Bundled Data
```
amplifinder/data/
├── ISfinderDB/IS.fna               # ISfinder sequences
└── breseq_fields/*.csv             # breseq output schemas
```

### Reference Cache (`--ref-path`, default: `genomesDB/`)

Naming convention:
- `{accession}` = user-provided NCBI accession (e.g., CP000000)
- `{locus}` = actual genome locus name from NCBI GenBank file (e.g., NC_000000)
- `{accession}.json` contains: `{"accession": "CP000000", "name": "NC_000000"}`

```
{ref_path}/
│
│ # NCBI downloads (GetRefGenomeStep)
├── {accession}.json                # mapping: accession → locus name
├── fasta/{locus}.fasta             # genome sequence
├── genbank/{locus}.gb              # GenBank file
│
│ # TN location processing
└── tn_loc/{locus}/
    ├── genbank_tn_loc.csv        # TN parsed from GenBank annotations (always)
    ├── isfinder_blast.txt        # BLAST output (always)
    ├── isfinder_tn_loc.csv       # TN parsed from BLAST results (always)
    ├── genbank_ref_tn_jc.csv     # reference TN junctions (selected source only)
    ├── genbank_tn_end_seqs.csv   # TN boundary sequences (selected source only)
    ├── isfinder_ref_tn_jc.csv    # reference TN junctions (selected source only)
    └── isfinder_tn_end_seqs.csv  # TN boundary sequences (selected source only)
```

### Run Output (`--output-dir`, default: `output/`)

Organized by: reference → ancestor → isolate

```
{output}/
└── {ref_name}/                     # reference genome (e.g., U00096)
    └── {anc_name}/                 # ancestor defines the group
        │
        │ # Ancestor folder (stores ancestor breseq and alignments for sharing)
        ├── {anc_name}/                 # ancestor folder
        │   ├── breseq/
        │   │   └── output/output.gd    # ancestor breseq output (Step 7)
        │   └── jc_{start}_{end}_{tn_id}_L{read_len}/  # per-candidate (from isolate runs)
        │       ├── junctions.fasta     # copied from isolate run
        │       └── iso.sorted.bam      # ancestor reads aligned to junctions
        │
        │ # Isolate runs (normalized coverage: iso/anc ratio)
        ├── {iso_name_1}/               # isolate 1 vs this ancestor
        │   ├── breseq/
        │   ├── tn_jc.csv
        │   ├── tn_jc2.csv
        │   ├── tn_jc2_*.csv
        │   ├── run_config.yaml
        │   ├── ISJC2.csv
        │   ├── candidate_amplifications.csv
        │   └── jc_{start}_{end}_{tn_id}_L{read_len}/  # per-candidate
        │       ├── junctions.fasta     # created from isolate candidates
        │       ├── iso.sorted.bam      # isolate reads aligned to junctions
        │       ├── anc.sorted.bam      # copied from {anc_name}/jc_.../iso.sorted.bam
        │       └── coverage_plot.png   # (if --save-plots)
        │
        └── {iso_name_2}/               # isolate 2 vs this ancestor
            └── jc_{start}_{end}_{tn_id}_L{read_len}/
                ├── junctions.fasta     # created from isolate candidates
                ├── iso.sorted.bam      # isolate reads aligned to junctions
                └── anc.sorted.bam      # copied from {anc_name}/jc_.../iso.sorted.bam (shared)
```

**Per-candidate folder naming:** `jc_{start}_{end}_{tn_id}_L{read_len}`
- `start`, `end`: amplicon coordinates (genomic positions)
- `tn_id`: TN index from TnLoc table (e.g., 001, 002)
- `read_len`: read length (affects junction flank size)

**Junction and alignment workflow:**
1. Isolate run creates junction files (in `{iso_name}/jc_.../`):
   - `junctions.fasta` (created from isolate candidates)
2. Junction files are copied to ancestor folder (in `{anc_name}/jc_.../`):
   - `junctions.fasta` (copied from isolate, allows ancestor alignments to be shared)
3. Ancestor reads are aligned in ancestor folder:
   - `{anc_name}/jc_.../iso.sorted.bam` (ancestor reads aligned to junctions)
4. Ancestor alignments are copied back to isolate folder:
   - `{anc_name}/jc_.../iso.sorted.bam` → `{iso_name}/jc_.../anc.sorted.bam`
5. Isolate reads are aligned in isolate folder:
   - `{iso_name}/jc_.../iso.sorted.bam` (isolate reads aligned to junctions)

**Key design:**
- Junction files are created from isolate candidates (not from ancestor)
- Ancestor alignments are stored in ancestor folder and can be reused for multiple isolate runs
- Only ancestor breseq output (Step 7) and ancestor reads (Step 11) are needed from ancestor

**Coverage:**
- Isolate run without ancestor (`anc_path=None`): RAW coverage (no normalization)
- Isolate run with ancestor (`anc_path=set`): iso/anc RATIO (uses ancestor breseq coverage from Step 7)

**Ancestor requirements:**
- Only ancestor breseq output (Step 7) is needed for coverage calculation
- Ancestor reads (Step 11) are aligned to isolate candidate junctions
- Ancestor alignments are stored in ancestor folder and can be shared across multiple isolate runs

## Implemented Steps

| Step | Python Class | Description | MATLAB Source |
|------|--------------|-------------|---------------|
| 0 | `InitializingStep` | Creates output directory | `makeDirs.m` |
| 1 | `GetRefGenomeStep` | Downloads genome FASTA/GBK | `get_reference.m`, `efetch_genbank.m` |
| 2 | `LocateTNsUsingGenbankStep` | Finds TN from GenBank annotations | `findISinRef.m` |
| 2 | `LocateTNsUsingISfinderStep` | Finds TN via BLAST to ISfinder DB | `ISfinder.m` |
| 3 | `CreateRefTnJcStep` | Creates junctions for ref TN | `create_JC_of_reference_IS.m` |
| 3 | `CreateRefTnEndSeqsStep` | TN flanking sequences for matching | `create_IS_end_seqs.m` |
| 4 | `BreseqStep` | Runs breseq on FASTQ + reference | `run_breseq.m` |
| 5 | `CreateTnJcStep` | Matches breseq JC to TN elements | `assign_potential_ISs.m` |
| 6 | `CreateTnJc2Step` | Pairs junctions into TnJc2 | `combine_ISJC_pairs.m`, `calculate_amplicon_length.m` |
| 7 | `CalcAmpliconCoverageStep` | Calculates amplicon coverage | `calc_coverage_ISJC2.m` |
| 8 | `ClassifyStructureStep` | Classifies junction pairs | `classify_ISJC2.m`, `directed_IS.m` |
| 9 | `FilterCandidatesStep` | Filters by amplicon length | `curate_candidate_amplicons.m` |
| 10 | `CreateSyntheticJunctionsStep` | Creates synthetic junction FASTA | `junction2fasta.m` |
| 11 | `AlignReadsToJunctionsStep` | Aligns reads to junctions | `bowtie2_alignment.m` |
| 12 | `AnalyzeAlignmentsStep` | Parses BAM, counts reads | `bamAnalysis.m`, `count_reads.m` |
| 13 | `ClassifyCandidatesStep` | Classifies based on read patterns | `classify_candidates.m` |
| 14 | `ExportStep` | Exports to Excel | `export_ISJC2.m` |

## Convergence Points (anc_path=None vs anc_path=set)

**Step 7 (CalcAmpliconCoverage):**
- `anc_path=None` → raw coverage only
- `anc_path=set` → load anc coverage, compute amplicon_coverage ratio

**Steps 10-11 (CreateSyntheticJunctions, AlignReads):**
- `anc_path=None` → create junctions, align isolate reads only
- `anc_path=set` → create junctions (from isolate candidates), copy to ancestor folder, align ancestor reads in ancestor folder, copy ancestor alignments back to isolate folder

**Step 12 (AnalyzeAlignments):**
- `anc_path=None` → parse own BAM → `jc_cov_*` fields
- `anc_path=set` → parse BOTH BAMs → `iso_jc_cov_*` + `anc_jc_cov_*` fields

**Step 13 (ClassifyCandidates):**
- `anc_path=None` → limited classification
- `anc_path=set` → full iso vs anc pattern comparison

## Utilities

### Coverage Parser (`utils/coverage.py`)
- `load_breseq_coverage()`: Load coverage from breseq output
- `get_coverage_in_range()`: Extract coverage in genomic range
- `calc_coverage_stats()`: Calculate mean/median/mode
- `calc_genome_stats()`: Genome-wide stats for normalization

### BAM Parser (`utils/bam.py`)
- `parse_bam_reads()`: Parse BAM, extract read positions and CIGAR
- `count_junction_reads()`: Count spanning/left/right reads per junction

### Bowtie2 Wrapper (`tools/bowtie2.py`)
- `run_bowtie2_build()`: Build bowtie2 index from FASTA
- `run_bowtie2_align()`: Align reads to index, output sorted BAM

## Development Guidelines

### Column Naming Consistency

The original MATLAB code has inconsistent column naming. **Best practice:**
1. Use **record inheritance** to ensure column names are consistent across all downstream tables
2. Define column names **once** in the base record and inherit them - never rename inherited fields
3. Only add **new** fields in derived records, never duplicate existing fields with different names
4. Apply renaming only at **export time** (Step 14) if user-friendly names are needed for Excel output

### Reference and Ancestor Context

Every output table should be self-documenting:
- `ref_name: str` in `TnJc2` (base record) - all junctions are relative to a reference
- `iso_name: str` in `TnJc2` - identifies which sample the junctions came from
- `anc_name: Optional[str]` in `CoveredTnJc2` - present when ancestor comparison is performed

### `raw_event` vs `isolate_architecture`

Both fields use similar vocabulary but represent **different analyses**:
- `raw_event` (Step 8): **Structural prediction** based on junction pair relationships
- `isolate_architecture` (Step 13): **Empirical verification** based on read alignment patterns

**Keep both fields** - comparing them can reveal discrepancies between predicted and observed architectures.

## Testing Patterns

### Test Structure
```
tests/
├── conftest.py              # Shared fixtures
├── env.py                   # Environment flags (e.g., RUN_BOWTIE2_TESTS)
├── fixtures/                # Test data files (breseq output, etc.)
├── test_data_types/         # Record/DataFrame tests
├── test_steps/              # Step unit tests
├── test_tools/              # Tool tests (blast, breseq_parser)
├── test_utils/              # Utility tests (tn_loc, coverage, bam)
└── test_integration/        # MATLAB comparison tests
    ├── conftest.py          # Integration test fixtures
    ├── matlab_compare.py    # MATLAB comparison utilities
    └── test_pipeline.py     # End-to-end pipeline tests
```

### Integration Tests (MATLAB Comparison)

Integration tests compare Python outputs against MATLAB reference outputs:

```python
# tests/test_integration/test_tnjc2.py
import scipy.io as sio

def test_tnjc2_matches_matlab(tmp_path):
    """Compare CreateTnJc2Step output to MATLAB reference."""
    # Load MATLAB reference
    matlab_data = sio.loadmat("AmpliFinder_test/ISJC2.mat")
    
    # Run Python pipeline up to this step
    # ...
    
    # Compare key fields
    assert python_df["amplicon_start"].tolist() == matlab_expected["amplicon_start"]
```

Reference data location: `AmpliFinder_test/` (MATLAB outputs for known test cases)

### Skip Markers for External Tools

```python
# tests/env.py
RUN_BOWTIE2_TESTS = shutil.which("bowtie2") is not None

# tests/test_tools/test_bowtie2.py
skip_no_bowtie = pytest.mark.skipif(not RUN_BOWTIE2_TESTS, reason="bowtie2 not found")
```

## Quick Reference

### System Tools & Libraries
- **bowtie2** - use `shutil.which("bowtie2")` for discovery
- **pysam** (v0.23.3) - replaces MATLAB BioMap + samtools
- **PyYAML** - for run_config.yaml
- **matplotlib** (v3.7+) - for coverage visualization
- **scipy** - for reading MATLAB files (`scipy.io.loadmat()`)

### MATLAB → Python Mapping

| MATLAB | Python |
|--------|--------|
| `bowtie2-build`, `bowtie2` | `subprocess.run([shutil.which("bowtie2"), ...])` |
| `samtools view/sort/index` | `pysam.AlignmentFile(...).write()` |
| `BioMap(bamfile)` | `pysam.AlignmentFile(bam_path, "rb")` |
| `load(COV.mat)` | `pd.read_csv(...coverage.tab)` |
| config file I/O | `yaml.safe_load()` / `yaml.dump()` |

### Key File Paths

```python
# Directory structure
anc_group_dir = output_dir / ref_name / anc_name      # ancestor group level
anc_run_dir = anc_group_dir / anc_name                # ancestor's own run
iso_run_dir = anc_group_dir / iso_name                # isolate run

# breseq coverage file
breseq_path / "08_mutation_identification" / f"{ref_name}.coverage.tab"

# run config
iso_run_dir / "run_config.yaml"

# per-candidate directory naming
jc_name = f"jc_{start}_{end}_{tn_id:03d}_L{read_len}"

# Ancestor run's candidate dir (stores ancestor alignments for sharing)
anc_jc_dir = anc_run_dir / jc_name
anc_jc_dir / "junctions.fasta"           # copied from isolate run
anc_jc_dir / "iso.sorted.bam"            # ancestor reads aligned to junctions

# Isolate run's candidate dir
iso_jc_dir = iso_run_dir / jc_name
iso_jc_dir / "junctions.fasta"           # created from isolate candidates
iso_jc_dir / "anc.sorted.bam"            # copied from anc_jc_dir/iso.sorted.bam
iso_jc_dir / "iso.sorted.bam"            # isolate reads aligned to junctions
iso_jc_dir / "coverage_plot.png"         # (if --save-plots)
```

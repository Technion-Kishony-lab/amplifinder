# AmpliFinder Remaining Steps Implementation Plan

## Overall Plan

```
{xxx}  input
xxx[]  array
- - -  optional
(?)    optional

                                 INPUTS
                 ┌────────────────┬────────────────┐
                 │                │                │
            {ref_name}         {FASTQ}         {anc_path}
                 │                │            (optional)
                 │                │                │
                 ▼                │                │
       ┌───────────────────┐      │                │
       │ 1. GetReference   │      │                │
       └─────────┬─────────┘      │                │
                 ▼                │                │
              Genome              │                │
          (FASTA + GBK)           │                │
                 │                │                │
        ┌────────┼────────┐       │                │
        │        │        │       │                │
        ▼        │        ▼       ▼                │
 ┌────────────┐  │  ┌───────────────────┐          │
 │ 2. LocateTN│  │  │ 4. BreseqStep     │          │
 └──────┬─────┘  │  └─────────┬─────────┘          │
        ▼        │            ▼                    │
     TnLoc[]     │       breseq JC                 │
        │        │       coverage                  │
        ▼        │            │                    │
 ┌────────────┐  │            │                    │
 │ 3. RefTnJC │  │            │                    │
 │ + EndSeqs  │  │            │                    │
 └──────┬─────┘  │            │                    │
        ▼        │            │                    │
 RefTnJunction[] │            │                    │
 TnEndSeq[]      │            │                    │
        │        │            │                    │
        └────────┴─────┬──────┘                    │
                       ▼                           │
             ┌───────────────────┐                 │
             │ 5. CreateTNJCStep │                 │
             └─────────┬─────────┘                 │
                       ▼                           │
                 TnJunction[]                      │
                       │                           │
                       ▼                           │
             ┌───────────────────┐                 │
             │ 6. CreateTNJC2Step│                 │
             └─────────┬─────────┘                 │
                       ▼                           │
                TnJunctionPair[]                   │
                       │                           │
                       ▼                           │
        ┌─────────────────────────────┐            │
        │ 7. CalcAmpliconCoverageStep │◄- - - - - -┤
        └──────────────┬──────────────┘    anc     │  
                       ▼                coverage   │
                 CoveredPair[]                     │
                       │                           │ 
                       ▼                           │
             ┌───────────────────┐                 │
             │ 8. ClassifyStruct │                 │
             └─────────┬─────────┘                 │
                       ▼                           │
                ClassifiedPair[] (with `RawEvent`) │
                       │                           │
                       ▼                           │
            ┌─────────────────────┐                │
            │ 9. FilterCandidates │                │
            └──────────┬──────────┘                │
                       ▼                           │
                  Candidate[]                      │
                       │                           │
                       ▼                           │
             ┌───────────────────┐                 │
             │10. SyntheticJuncs │                 │
             └─────────┬─────────┘                 │
                       ▼                           │
               junctions.fasta                     │
        (7 `JunctionType` per candidate)           │
                       │                           │
                       ▼                           │
             ┌───────────────────┐                 │
 {FASTQ} ───►│ 11. AlignReads    │◄- - - - - - - - ┘
             └─────────┬─────────┘    {anc_FASTQ}
                       ▼             
                   iso.bam           
                   anc.bam(?)          
                       │             
                       ▼             
             ┌───────────────────┐   
             │12. AnalyzeAligns  │
             └─────────┬─────────┘
                       ▼
              AnalyzedCandidate[]
                c_cov, anc_jc_cov(?)
                       │
                       ▼
             ┌───────────────────┐
             │13. ClassifyCands  │
             └─────────┬─────────┘
                       ▼
              AnalyzedCandidate[] (with event: Tuple[RawEvent, List[EventModifier]])
                       │
                       ▼
             ┌───────────────────┐
             │ 14. Export        │
             └─────────┬─────────┘
                       ▼
                  ISJC2.xlsx
          candidate_amplifications.xlsx
```

## Folder Structure

```
# Bundled data
amplifinder/data/
├── ISfinderDB/IS.fna               # ISfinder sequences
└── breseq_fields/*.csv             # breseq output schemas

# Reference cache (--ref-path, default: output/reference/)
{ref_path}/
├── {ref_name}.fasta                # genome sequence
├── {ref_name}.gbk                  # GenBank annotations
├── {ref_name}.json                 # genome metadata
├── {ref_name}_TN_end_seqs.csv      # TN boundary sequences
├── genbank/{ref_name}_tn_loc.csv   # TN from GenBank
└── isfinder/{ref_name}_tn_loc.csv  # TN from ISfinder

# Run output (--output)
{output}/
└── {iso_name}/
    ├── breseq/                     # breseq output
    │   └── output/output.gd
    ├── ref_tn_jc.csv               # reference TN junctions
    ├── TNJC.csv                    # TN-associated junctions
    ├── TNJC2.csv                   # junction pairs (candidates)
    └── TNJC2_*/                    # per-candidate analysis (future)
```


## What We Already Have

### Base Abstractions

**`Step[T]`** (`steps/base.py`): Abstract pipeline step with:
- Input/output file tracking
- Caching: skip if outputs exist (unless `force=True`)
- Methods: `_calculate_output() → T`, `_save_output()`, `load_outputs()`

**`Record`** (`data_types/records.py`): Pydantic model with:
- Auto-generated `schema()` from fields
- `from_other()` to convert between record types (copy shared fields)
- Extra fields support

### Implemented Steps (0-6)

0. **InitializingStep** — creates output directory (`makeDirs.m`)
1. **GetReferenceStep** — downloads genome FASTA/GBK (`get_reference.m`, `efetch_genbank.m`)
2. **LocateTNsUsingGenbankStep** — finds TN from GenBank annotations (`findISinRef.m`)
   **LocateTNsUsingISfinderStep** — finds TN via BLAST to ISfinder DB (`ISfinder.m`)
3. **CreateReferenceTnJunctionsStep** — synthetic junctions for ref TN (`create_JC_of_reference_IS.m`)
   **CreateRefTnEndSeqsStep** — TN boundary sequences for matching (`create_IS_end_seqs.m`)
4. **BreseqStep** — runs breseq on FASTQ + reference (`run_breseq.m`)
5. **CreateTNJCStep** — matches breseq JC to TN elements (`assign_potential_ISs.m`)
6. **CreateTNJC2Step** — pairs junctions into candidates (`combine_ISJC_pairs.m`, `calculate_amplicon_length.m`)



## Isolate/Ancestor Pipeline Design

### Core Pattern
Run the **same pipeline** on both ancestor and isolate:
- `anc_path: Optional[Path]` - path to completed ancestor run (None for ancestor run)
- `anc_fastq_path` is derived from ancestor's saved config (not a separate CLI arg)

```bash
# Option A: Two separate runs (useful when sharing ancestor across isolates)
# 1. Run ancestor first (standalone)
amplifinder --sample ancestor.fastq --ref U00096 --output output/ancestor/
# Saves config to output/ancestor/run_config.yaml including sample_path

# 2. Run isolate with ancestor reference (only need --anc-path)
amplifinder --sample isolate.fastq --ref U00096 --output output/isolate/ \
            --anc-path output/ancestor/
# Reads anc_fastq_path from output/ancestor/run_config.yaml

# Option B: Single command for both (runs ancestor first, then isolate)
amplifinder --sample isolate.fastq --anc-sample ancestor.fastq --ref U00096 \
            --output output/isolate/
# Automatically: 1) runs ancestor to output/isolate/ancestor/, 2) runs isolate with --anc-path
```

### Config Saving
`InitializingStep` saves config to `{output}/run_config.yaml`:

```yaml
sample_path: /path/to/sample.fastq
ref_name: U00096
breseq_path: /path/to/breseq/output
# ... other params
```

```python
# steps/initialize.py - add config param, call save_run_config()
# utils/run_config.py - save_run_config(), load_run_config()
```

### Convergence Points

```
Step 7:  anc_path=None → raw coverage only
         anc_path=set  → load anc coverage, compute amplicon_coverage ratio

Steps 10-11: anc_path=None → create junctions, align THIS sample only
             anc_path=set  → create junctions, align BOTH iso reads AND anc reads (anc_fastq from saved config)

Step 12: anc_path=None → parse own BAM → jc_cov_* fields  
         anc_path=set  → parse BOTH BAMs → iso_jc_cov_* + anc_jc_cov_* fields

Step 13: anc_path=None → limited classification
         anc_path=set  → full iso vs anc pattern comparison
```

---

## Record Type Progression

```python
# record_types.py

class TnJunctionPair(Record):       # Step 6 output (existing)
    # ...existing fields...
    # NOTE: Remove these optional fields from TnJunctionPair:
    #   amplicon_coverage: Optional[float] = None
    #   copy_number_mode: Optional[float] = None
    # They belong in CoveredPair

class CoveredPair(TnJunctionPair):  # Step 7 output
    # This sample's coverage stats (normalized by genome median coverage)
    amplicon_coverage_median: float   # median coverage in amplicon region
    amplicon_coverage_mode: float     # mode coverage in amplicon region
    copy_number_median: float         # amplicon_median / genome_median
    copy_number_mode: float           # amplicon_mode / genome_mode
    # Ratio vs ancestor (only when iso run, None for anc run)
    copy_number_ratio_median: Optional[float] = None  # iso_median / anc_median
    copy_number_ratio_mode: Optional[float] = None    # iso_mode / anc_mode

class RawEvent(str, Enum):
    """Structure classification from classify_ISJC2.m"""
    REFERENCE = "reference"
    TRANSPOSITION = "transposition"
    UNFLANKED = "unflanked"
    HEMI_FLANKED_LEFT = "hemi-flanked left"
    HEMI_FLANKED_RIGHT = "hemi-flanked right"
    FLANKED = "flanked"
    MULTIPLE_SINGLE_LOCUS = "multiple single locus"
    UNRESOLVED = "unresolved"

class EventModifier(str, Enum):
    """Modifiers for final event classification."""
    ANCESTRAL = "ancestral"
    DE_NOVO = "de novo"
    LOW_COVERAGE = "low coverage near junction"

class ClassifiedPair(CoveredPair):  # Step 8 output
    raw_event: RawEvent
    shared_IS: List[int]
    chosen_IS: int
    pos_of_paired_single_locus_junction: Tuple[int, int]

class Candidate(ClassifiedPair):    # Step 9 output
    analysis_directory: str         # "TNJC2_001"

class JunctionType(int, Enum):
    """The 7 synthetic junction types from junction2fasta.m
    
    Amplicon structure: ~~~>>>======>>>======>>>~~~
    """
    LEFT_REF = 1        # ~~==  left reference (chromosome-cassette)
    LEFT_IS_TRANS = 2   # ~~>>  left IS transposition (chromosome-IS)
    LEFT_MID_IS = 3     # ==>>  left of mid IS (cassette-IS)
    LOST_IS = 4         # ====  lost IS (cassette-cassette, no IS)
    RIGHT_MID_IS = 5    # >>==  right of mid IS (IS-cassette)
    RIGHT_IS_TRANS = 6  # >>~~  right IS transposition (IS-chromosome)
    RIGHT_REF = 7       # ==~~  right reference (cassette-chromosome)

class JunctionCoverage(NamedTuple):
    """Read counts for a single junction type."""
    spanning: int  # reads crossing junction (equiv to MATLAB 'green')
    left: int      # reads ending at junction
    right: int     # reads starting at junction

class AnalyzedCandidate(Candidate): # Step 12-13 output
    # Each list has 7 elements pf JunctionCoverage, one per JunctionType (1-7)
    jc_cov: List[JunctionCoverage]                    # THIS sample's coverage (always)
    anc_jc_cov: Optional[List[JunctionCoverage]]      # ancestor's coverage (only when iso run)
    # Final classification: (base_event, modifiers)
    event: Tuple[RawEvent, List[EventModifier]]
    isolate_architecture: RawEvent
    ancestor_architecture: Optional[RawEvent]
```

---

## Implementation Steps

Each step: implement → write tests → run/debug → flake8 → git commit

### 0. Record Types & Config Utils (prerequisite)
- Module: `data_types/record_types.py`
- Changes:
  - Remove `amplicon_coverage` and `copy_number_mode` from `TnJunctionPair`
  - Add: `CoveredPair`, `ClassifiedPair`, `Candidate`, `AnalyzedCandidate`
  - Add: `RawEvent`, `EventModifier`, `JunctionType`, `JunctionCoverage` enums
- Module: `utils/run_config.py`
  - Add: `save_run_config()`, `load_run_config()`
- Module: `steps/initialize.py`
  - Update: accept `config: Config`, call `save_run_config()` after creating dir
- Tests: existing tests + `test_utils/test_run_config.py`
- Commit: "Add remaining Record types and run config utils"

### 1. Coverage Parser (utility)
- Module: `utils/coverage.py`
- MATLAB: `parseBreseqCOV.m`, `get_coverage_in_range.m`, `calc_average_coverage.m`
- Functions:
  ```python
  def load_breseq_coverage(breseq_path: Path, ref_name: str) -> np.ndarray:
      """Load coverage from breseq output.
      
      Returns 1D array of coverage values (top + bottom strand).
      """
      cov_file = breseq_path / "08_mutation_identification" / f"{ref_name}.coverage.tab"
      df = pd.read_csv(cov_file, sep="\t", usecols=["unique_top_cov", "unique_bot_cov"])
      return (df["unique_top_cov"] + df["unique_bot_cov"]).values

  def get_coverage_in_range(cov: np.ndarray, start: int, end: int) -> np.ndarray:
      """Get coverage values in range [start, end] (1-based, inclusive)."""
      return cov[start - 1 : end]

  def calc_coverage_stats(cov: np.ndarray) -> Tuple[float, float]:
      """Calculate median and mode of coverage array.
      
      Returns (median, mode).
      """

  def calc_genome_median(cov: np.ndarray) -> float:
      """Calculate genome-wide median coverage for normalization."""
  ```
- Tests: `test_utils/test_coverage.py`
- Commit: "Add breseq coverage parser"

### 2. CalcAmpliconCoverageStep (Step 7)
- Module: `steps/calc_coverage.py`
- MATLAB: `calc_coverage_ISJC2.m`
- Input: `TnJunctionPair` → Output: `CoveredPair`
- File output: None (memory only)
- Logic:
  - Get coverage in amplicon region from breseq coverage files
  - If `anc_path=None`: raw coverage only
  - If `anc_path=set`: compute `amplicon_coverage = iso / anc`
- Tests: `test_steps/test_calc_coverage.py`
- Commit: "Implement CalcAmpliconCoverageStep"

### 3. ClassifyStructureStep (Step 8)
- Module: `steps/classify_structure.py`
- MATLAB: `classify_ISJC2.m`, `directed_IS.m`
- Input: `CoveredPair` → Output: `ClassifiedPair`
- File output: None (memory only)
- Logic: Classify as transposition/reference/unflanked/hemi-flanked/flanked
- Tests: `test_steps/test_classify_structure.py`
- Commit: "Implement ClassifyStructureStep"

### 4. FilterCandidatesStep (Step 9)
- Module: `steps/filter_candidates.py`
- MATLAB: filtering logic in `curate_candidate_amplicons.m`
- Input: `ClassifiedPair` → Output: `Candidate`
- File output: None (memory only)
- Logic: Filter by `MIN_AMPLICON_LENGTH < amplicon_length < MAX_AMPLICON_LENGTH`
- Tests: `test_steps/test_filter_candidates.py`
- Commit: "Implement FilterCandidatesStep"

### 5. CreateSyntheticJunctionsStep (Step 10)
- Module: `steps/synthetic_junctions.py`
- MATLAB: `junction2fasta.m`
- Input: `List[Candidate]`, genome sequences, TnLoc
- Output: `List[Candidate]` (unchanged, side effect: creates files)
- File output per candidate: `{analysis_dir}/junctions.fasta` (7 synthetic junctions)
- Logic: Create junction sequences for all expected architectural variations
- Tests: `test_steps/test_synthetic_junctions.py`
- Commit: "Implement CreateSyntheticJunctionsStep"

### 6. Bowtie2 Wrapper (tool)
- Module: `tools/bowtie2.py`
- MATLAB: `run_bowtie_line.m`
- System: bowtie2 (v2.3.5.1+), bowtie2-build
- Functions:
  ```python
  import shutil
  
  BOWTIE2 = shutil.which("bowtie2")
  BOWTIE2_BUILD = shutil.which("bowtie2-build")

  def run_bowtie2_build(ref_fasta: Path, index_prefix: Path) -> None:
      """Build bowtie2 index from FASTA."""

  def run_bowtie2_align(
      index: Path, 
      fastq_path: Path, 
      output_bam: Path,  # Direct BAM output via pysam
      score_min: str = "L,0,-0.6",
      num_alignments: int = 10,
      threads: int = 1,
  ) -> None:
      """Align reads to index, output sorted BAM."""
  ```
- Tests: `test_tools/test_bowtie2.py` (skip if bowtie2 not found)
- Commit: "Add bowtie2 wrapper"

### 7. AlignReadsStep (Step 11)
- Module: `steps/align_reads.py`
- MATLAB: `bowtie2_alignment.m`, `run_samtools_line.m`
- Input: `Candidate` + FASTQ → Output: `Candidate`
- File output: `{analysis_dir}/alignment/alignment.sorted.bam`
- Logic:
  - Align this sample's reads to synthetic junctions
  - If `anc_path` set: also align ancestor reads (from saved config)
- Tests: `test_steps/test_align_reads.py`
- Commit: "Implement AlignReadsStep"

### 8. BAM Parser (utility)
- Module: `utils/bam.py`
- MATLAB: `parse_reads_from_bam.m` (uses BioMap)
- Library: `pysam` (v0.23.3)
- Functions:
  ```python
  def parse_bam_reads(bam_path: Path) -> List[dict]:
      """Parse BAM, extract read positions and CIGAR info."""
      # Returns list of {ref, start, length, cigar, is_full}

  def count_junction_reads(reads, seq_length, read_length, req_overlap=12) -> List[JunctionCoverage]:
      """Count spanning/left/right reads per junction (1-7)."""
      # Implements MATLAB count_reads.m logic
  ```
- Tests: `test_utils/test_bam.py`
- Commit: "Add BAM parser using pysam"

### 9. AnalyzeAlignmentsStep (Step 12)
- Module: `steps/analyze_alignments.py`
- MATLAB: `bamAnalysis.m`, `count_reads.m`
- Input: `Candidate` with BAM → Output: `AnalyzedCandidate`
- File output: None (memory only)
- Logic: Parse BAM, count green/left/right reads per junction (1-7)
- Tests: `test_steps/test_analyze_alignments.py`
- Commit: "Implement AnalyzeAlignmentsStep"

### 10. ClassifyCandidatesStep (Step 13)
- Module: `steps/classify_candidates.py`
- MATLAB: `classify_candidates.m`
- Input: `AnalyzedCandidate` → Output: `AnalyzedCandidate`
- File output: None (memory only)
- Logic: Match junction read patterns to expected architectures
- Tests: `test_steps/test_classify_candidates.py`
- Commit: "Implement ClassifyCandidatesStep"

### 11. ExportStep (Step 14)
- Module: `steps/export.py`
- MATLAB: `export_ISJC2.m`
- Input: `AnalyzedCandidate`
- File output: `ISJC2.xlsx`, `candidate_amplifications.xlsx`
- Logic: Format and export, filter by copy number thresholds
- Tests: `test_steps/test_export.py`
- Commit: "Implement ExportStep"

### 12. Pipeline Integration
- Module: `pipeline.py` + `cli.py`
- MATLAB: `AmpliFinder.m` (main orchestration)
- Add: `--anc-path`, `--anc-sample` CLI options
- Tests: `test_integration/test_pipeline.py`
- Commit: "Integrate remaining steps"

### 13. Full Integration Tests
- Tests: `test_integration/test_full_run.py`
- Commit: "Add full integration tests"

---

## Development Workflow

```
1. Implement module
     ↓
2. Write unit tests
     ↓
3. Run tests, debug until passing
     ↓
4. Run flake8, fix errors
     ↓
5. Git commit
     ↓
6. Continue to next step
```

---

## Test Patterns

Follow existing patterns in `tests/`:

```
tests/
├── conftest.py              # Shared fixtures
├── test_steps/              # Step tests
├── test_tools/              # Tool tests (bowtie2)
└── test_utils/              # Utility tests (coverage, bam)
```

### Skip Markers for External Tools

```python
# tests/env.py
RUN_BOWTIE2_TESTS = shutil.which("bowtie2") is not None

# tests/test_tools/test_bowtie2.py
skip_no_bowtie = pytest.mark.skipif(not RUN_BOWTIE2_TESTS, reason="bowtie2 not found")
```

---

## Quick Reference

### System Tools & Libraries
- **bowtie2** - use `shutil.which("bowtie2")` for discovery
- **pysam** (v0.23.3) - replaces MATLAB BioMap + samtools
- **PyYAML** - for run_config.yaml

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
# breseq coverage file
breseq_path / "08_mutation_identification" / f"{ref_name}.coverage.tab"

# run config
output_dir / "run_config.yaml"

# synthetic junctions per candidate
output_dir / f"TNJC2_{idx:03d}" / "junctions.fasta"

# alignment output per candidate
output_dir / f"TNJC2_{idx:03d}" / "alignment.sorted.bam"
```


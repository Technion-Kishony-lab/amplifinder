# AmpliFinder Remaining Steps Implementation Plan

## Overall Plan

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
   │        breseq JC        │                  TnEndSeq[]                   │
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
                                        ISJC2.xlsx
                                        candidate_amplifications.xlsx
```

## Folder Structure (partially implemented)

```
# Bundled data
amplifinder/data/
├── ISfinderDB/IS.fna               # ISfinder sequences
└── breseq_fields/*.csv             # breseq output schemas

# Reference cache (--ref-path, default: genomesDB/)
#
# Naming convention:
#   {accession} = user-provided NCBI accession (e.g., CP000000)
#   {locus}     = actual genome locus name from NCBI GenBank file (e.g., NC_000000)
#
#   {accession}.json contains: {"accession": "CP000000", "name": "NC_000000"}
#   This mapping allows lookup by accession to find locus-named files
#
{ref_path}/
│
│ # NCBI downloads (GetRefGenomeStep)
├── {accession}.json                # mapping: accession → locus name
├── fasta/{locus}.fasta             # genome sequence
├── genbank/{locus}.gb              # GenBank file
│
│ # TN location processing (LocateTNsUsingGenbankStep, LocateTNsUsingISfinderStep)
│ # and TN end sequences (CreateRefTnJcStep, CreateRefTnEndSeqsStep)
└── tn_loc/{locus}/
    ├── genbank_tn_loc.csv        # TN parsed from GenBank annotations (always)
    ├── isfinder_blast.txt        # BLAST output (always)
    ├── isfinder_tn_loc.csv       # TN parsed from BLAST results (always)
    ├── genbank_ref_tn_jc.csv     # reference TN junctions (selected source only)
    ├── genbank_tn_end_seqs.csv   # TN boundary sequences (selected source only)
    ├── isfinder_ref_tn_jc.csv    # reference TN junctions (selected source only)
    └── isfinder_tn_end_seqs.csv  # TN boundary sequences (selected source only)


# Run output (--output, default: output/)
# Organized by ancestor: all isolates sharing an ancestor live under that ancestor's folder
{output}/                           # --output CLI argument (default: output/)
└── {anc_name}/                     # ancestor defines the group
    ├── _self/                      # ancestor's own standalone run (no anc comparison)
    │   ├── breseq/
    │   │   └── output/output.gd
    │   ├── tn_jc.csv
    │   ├── tn_jc2.csv
    │   └── run_config.yaml
    ├── {iso_name_1}/               # isolate 1 vs this ancestor
    │   ├── breseq/
    │   │   └── output/output.gd
    │   ├── tn_jc.csv
    │   ├── tn_jc2.csv
    │   └── run_config.yaml
    └── {iso_name_2}/               # isolate 2 vs this ancestor
        └── ...
```


## What We Already Have

### Base Abstractions

**`Step[T]`** (`steps/base.py`): Abstract pipeline step with:
- Input/output file tracking
- Caching: skip if outputs exist (unless `force=True`)
- Methods: `_calculate_output() → T`, `_save_output()`, `load_outputs()`
- To implement steps without file outputs (memory only), just set a class with `output_files=None`.

**`Record`** (`data_types/records.py`): Pydantic model with:
- Auto-generated `schema()` from fields
- `from_other()` to convert between record types (copy shared fields)
- Extra fields support

### Implemented Steps (0-6)

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



## New plan (not yet implemened)

### New Pattern - run the **same pipeline** on both ancestor and isolate:
- Output organized by ancestor: `{output}/{anc_name}/{iso_name}/`
- `iso_name` - sample name (required), the sample being analyzed
- `anc_name` - ancestor name (`None` → no ancestor comparison)

```bash
# Run ancestor standalone (no isolate comparison) - ancestor runs as an iso
amplifinder --sample my_ancestor.fastq --ref U00096 --iso-name my_ancestor
# Output: output/my_ancestor/_self   (no anc_name, so run as '_self')

# Run isolate with ancestor
amplifinder --sample my_isolate.fastq --ref U00096 --iso-name my_isolate --anc-name my_ancestor
# Output: output/my_ancestor/my_isolate/
# Auto-runs ancestor to output/my_ancestor/_self if not already done

# Run another isolate with same ancestor (reuses cached ancestor run)
amplifinder --sample my_isolate_2.fastq --ref U00096 --iso-name my_isolate_2 --anc-name my_ancestor
# Output: output/my_ancestor/my_isolate_2/
# Skips ancestor run (already exists at output/my_ancestor/_self)
```

### Convergence Points (anc_name=`None` vs anc_name=`{ancestor}`)

```
Step 7:  anc_name=None → raw coverage only
         anc_name=set  → load anc coverage, compute amplicon_coverage ratio

Steps 10-11: anc_name=None → create junctions, align THIS sample only
             anc_name=set  → create junctions, align BOTH iso reads AND anc reads

Step 12: anc_name=None → parse own BAM → jc_cov_* fields  
         anc_name=set  → parse BOTH BAMs → iso_jc_cov_* + anc_jc_cov_* fields

Step 13: anc_name=None → limited classification
         anc_name=set  → full iso vs anc pattern comparison
```

---

## Implementation Steps

IMPORTANT: For each implementation step, follow this scheme:
A. implement
B. write unit tests - run/debug 
C. write integration test with Matlab comparison (test_integration, AmpliFinder_test) - run/debug
D. flake8 - run and fix
E. git commit

### 0. Config Utils
- Module: `config.py` (existing `Config` class already has `iso_name`, `anc_name`, `iso_path`, `ref_name`, `output_dir`)
  - Add: `get_run_dir(config)` → path resolution
  - Add: `save_config(config, run_dir)`, `load_config_from_run(run_dir)`
- Module: `steps/initialize.py` - update to call `save_config()` after creating dir
- Output: `{run_dir}/run_config.yaml`
- Tests: `test_utils/test_run_config.py`
- Commit: "Add config save/load utils"

### 1. Coverage Parser (utility)
- Module: `utils/coverage.py`
- MATLAB: `parseBreseqCOV.m`, `get_coverage_in_range.m`, `calc_average_coverage.m`
- **Types to add:**
  ```python
  class Coverage(NamedTuple):
      mean: float
      median: float
      mode: float
  ```
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

  def calc_coverage_stats(cov: np.ndarray) -> Coverage:
      """Calculate mean, median and mode of coverage array."""

  def calc_genome_stats(cov: np.ndarray) -> Coverage:
      """Calculate genome-wide coverage stats for normalization."""
  ```
- Tests: `test_utils/test_coverage.py`
- Commit: "Add breseq coverage parser"

### 2. CalcAmpliconCoverageStep (Step 7)
- Module: `steps/calc_amplicon_coverage.py`
- MATLAB: `calc_coverage_ISJC2.m`
- **Types:**
  ```python
  class CoveredTnJc2(TnJc2):  # Record
      amplicon_coverage: Coverage
      genome_coverage: Coverage
      copy_number: Coverage             # amplicon / genome
      copy_number_ratio: Optional[Coverage] = None  # iso / anc (only when iso-anc run)
  ```
- **input_files:**
  - `breseq_path`: Path to breseq output (contains `08_mutation_identification/{ref_name}.coverage.tab`)
  - `anc_breseq_path`: Optional[Path] - ancestor breseq output (when `anc_path` is set)
- **data_inputs:**
  - `tnjc2_list`: RecordTypedDF[TnJc2]
- **Output:** RecordTypedDF[CoveredTnJc2]
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
- **Record types:**
  ```python
  class RawEvent(str, Enum):
      REFERENCE = "reference"
      TRANSPOSITION = "transposition"
      UNFLANKED = "unflanked"
      HEMI_FLANKED_LEFT = "hemi-flanked left"
      HEMI_FLANKED_RIGHT = "hemi-flanked right"
      FLANKED = "flanked"
      MULTIPLE_SINGLE_LOCUS = "multiple single locus"
      UNRESOLVED = "unresolved"

  class ClassifiedTnJc2(CoveredTnJc2):
      raw_event: RawEvent
      shared_IS: List[int]
      chosen_IS: int
      pos_of_paired_single_locus_junction: Tuple[int, int]
  ```
- **input_files:** None
- **data_inputs:**
  - `covered_tnjc2_list`: RecordTypedDF[CoveredTnJc2]
  - `tn_locs`: RecordTypedDF[TnLoc] - for IS direction lookup
- **Output:** RecordTypedDF[ClassifiedTnJc2]
- File output: None (memory only)
- Logic: Classify as transposition/reference/unflanked/hemi-flanked/flanked
- Tests: `test_steps/test_classify_structure.py`
- Commit: "Implement ClassifyStructureStep"

### 4. FilterCandidatesStep (Step 9)
- Module: `steps/filter_candidates.py`
- MATLAB: filtering logic in `curate_candidate_amplicons.m`
- **Record types:**
  ```python
  class CandidateTnJc2(ClassifiedTnJc2):
      analysis_directory: str  # "tn_jc2_001"
  ```
- **input_files:** None
- **data_inputs:**
  - `classified_tnjc2_list`: RecordTypedDF[ClassifiedTnJc2]
- **Output:** RecordTypedDF[CandidateTnJc2]
- File output: None (memory only)
- Logic: Filter by `MIN_AMPLICON_LENGTH < amplicon_length < MAX_AMPLICON_LENGTH`
- Tests: `test_steps/test_filter_candidates.py`
- Commit: "Implement FilterCandidatesStep"

### 5. CreateSyntheticJunctionsStep (Step 10)
- Module: `steps/synthetic_junctions.py`
- MATLAB: `junction2fasta.m`
- **Record types:**
  ```python
  class JunctionType(int, Enum):
      """The 7 synthetic junction types. Amplicon: ~~~>>>======>>>======>>>~~~"""
      LEFT_REF = '~~=='       # left reference (chromosome-cassette)
      LEFT_IS_TRANS = '~~>>'  # left IS transposition (chromosome-IS)
      LEFT_MID_IS = '==>>'    # left of mid IS (cassette-IS)
      LOST_IS = '===='        # lost IS (cassette-cassette, no IS)
      RIGHT_MID_IS = '>>=='   # right of mid IS (IS-cassette)
      RIGHT_IS_TRANS = '>>~~' # right IS transposition (IS-chromosome)
      RIGHT_REF = '==~~'      # right reference (cassette-chromosome)
  ```
- **input_files:**
  - `genome_fasta`: Path - reference genome FASTA
- **data_inputs:**
  - `candidates`: RecordTypedDF[CandidateTnJc2]
  - `tn_locs`: RecordTypedDF[TnLoc] - for IS sequences
- **Output:** RecordTypedDF[CandidateTnJc2] (unchanged, side effect: creates files)
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

### 7. AlignReadsToSyntheticJunctionsStep (Step 11)
- Module: `steps/align_reads.py`
- MATLAB: `bowtie2_alignment.m`, `run_samtools_line.m`
- **input_files:**
  - `fastq_path`: Path - this sample's FASTQ
  - `anc_fastq_path`: Optional[Path] - ancestor FASTQ (from saved config when `anc_path` is set)
  - `junctions_fasta`: Path - per-candidate `{analysis_dir}/junctions.fasta`
- **data_inputs:**
  - `candidates`: RecordTypedDF[CandidateTnJc2]
- **Output:** RecordTypedDF[CandidateTnJc2] (unchanged, side effect: creates BAM files)
- File output: `{analysis_dir}/iso.sorted.bam`, `{analysis_dir}/anc.sorted.bam` (optional)
- Logic:
  - Align this sample's reads to synthetic junctions
  - If `anc_path` set: also align ancestor reads (from saved config)
- Tests: `test_steps/test_align_reads.py`
- Commit: "Implement AlignReadsToSyntheticJunctionsStep"

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
- **Record types:**
  ```python
  class JunctionCoverage(NamedTuple):
      spanning: int  # reads crossing junction
      left: int      # reads ending at junction
      right: int     # reads starting at junction

  class AnalyzedCandidateTnJc2(CandidateTnJc2):
      jc_cov: List[JunctionCoverage]                    # 7 elements, one per JunctionType
      anc_jc_cov: Optional[List[JunctionCoverage]]      # ancestor's (only when iso run)
      event: Tuple[RawEvent, List[EventModifier]]
      isolate_architecture: RawEvent
      ancestor_architecture: Optional[RawEvent]
  ```
- **input_files:**
  - `iso_bam`: Path - `{analysis_dir}/iso.sorted.bam`
  - `anc_bam`: Optional[Path] - `{analysis_dir}/anc.sorted.bam` (when `anc_path` is set)
- **data_inputs:**
  - `candidates`: RecordTypedDF[CandidateTnJc2]
- **Output:** RecordTypedDF[AnalyzedCandidateTnJc2]
- File output: None (memory only)
- Logic: Parse BAM, count spanning/left/right reads per junction (1-7)
- Tests: `test_steps/test_analyze_alignments.py`
- Commit: "Implement AnalyzeAlignmentsStep"

### 10. ClassifyCandidatesStep (Step 13)
- Module: `steps/classify_candidates.py`
- MATLAB: `classify_candidates.m`
- **Record types:**
  ```python
  class EventModifier(str, Enum):
      ANCESTRAL = "ancestral"
      DE_NOVO = "de novo"
      LOW_COVERAGE = "low coverage near junction"
  ```
- **input_files:** None
- **data_inputs:**
  - `analyzed_candidates`: RecordTypedDF[AnalyzedCandidateTnJc2]
- **Output:** RecordTypedDF[AnalyzedCandidateTnJc2] (with `event`, `isolate_architecture`, `ancestor_architecture` populated)
- File output: None (memory only)
- Logic: Match junction read patterns to expected architectures
- Tests: `test_steps/test_classify_candidates.py`
- Commit: "Implement ClassifyCandidatesStep"

### 11. ExportStep (Step 14)
- Module: `steps/export.py`
- MATLAB: `export_ISJC2.m`
- **input_files:** None
- **data_inputs:**
  - `analyzed_candidates`: RecordTypedDF[AnalyzedCandidateTnJc2]
- **Output:** None (side effect: creates files)
- File output: `ISJC2.xlsx`, `candidate_amplifications.xlsx`
- Logic: Format and export, filter by copy number thresholds
- Tests: `test_steps/test_export.py`
- Commit: "Implement ExportStep"

### 12. Pipeline Integration
- Module: `pipeline.py` + `cli.py`
- MATLAB: `AmpliFinder.m` (main orchestration)
- **CLI args:** `--sample`, `--ref`, `--output`, `--iso-name`, `--anc-name` → `Config`
- **Logic:**
  - If `anc_name` set:
    - Check if ancestor run exists at `{output}/{anc_name}/_self/`
    - If missing → recursively run full pipeline with `iso_name=anc_name, anc_name=None`
    - Load ancestor config via `load_config_from_run()`
  - Run full pipeline for isolate
  - Pass ancestor results to convergence points (Steps 7, 10-12)
- Tests: `test_integration/test_pipeline.py`
- Commit: "Integrate remaining steps with ancestor support"

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
output_dir / f"tn_jc2_{idx:03d}" / "junctions.fasta"

# alignment output per candidate
output_dir / f"tn_jc2_{idx:03d}" / "alignment.sorted.bam"
```


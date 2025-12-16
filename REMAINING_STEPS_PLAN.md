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
                                             │
                                             ▼
                                ┌──────────────────────────┐
                                │ 15. VisualizeCandidates  │◄─── interactive
                                └────────────┬─────────────┘
                                             ▼
                                   coverage plots (matplotlib)
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
# Organized by: reference → ancestor → isolate
# Note: When no ancestor is provided, iso_name becomes its own anc_name
{output}/                           # --output CLI argument (default: output/)
└── {ref_name}/                     # reference genome (e.g., U00096)
    └── {anc_name}/                 # ancestor defines the group
        │
        │ # Ancestor run (same structure as isolates, RAW coverage)
        ├── {anc_name}/                 # ancestor run
        │   ├── breseq/
        │   │   └── output/output.gd
        │   ├── tn_jc.csv
        │   ├── tn_jc2.csv
        │   ├── tn_jc2_*.csv            # Output of different Steps
        │   ├── run_config.yaml
        │   └── jc_{start}_{end}_{tn_id}_L{read_len}/  # CANONICAL SOURCE
        │       ├── junctions.fasta     # created during ancestor run
        │       └── iso.sorted.bam      # ancestor's own alignment
        │
        │ # Isolate runs (normalized coverage: iso/anc ratio)
        ├── {iso_name_1}/               # isolate 1 vs this ancestor
        │   ├── breseq/
        │   │   └── output/output.gd
        │   ├── tn_jc.csv
        │   ├── tn_jc2.csv
        │   ├── tn_jc2_*.csv            # Output of different Steps
        │   ├── run_config.yaml
        │   └── jc_{start}_{end}_{tn_id}_L{read_len}/  # per-candidate
        │       ├── junctions.fasta     # copied from {anc_name}/{anc_name}/jc_.../
        │       ├── iso.sorted.bam      # isolate-specific alignment
        │       ├── anc.sorted.bam      # copied from {anc_name}/{anc_name}/jc_.../iso.sorted.bam
        │       └── coverage_plot.png   # (if --save-plots)
        │
        └── {iso_name_2}/               # isolate 2 vs this ancestor
            └── jc_{start}_{end}_{tn_id}_L{read_len}/
                ├── junctions.fasta     # copied from ancestor run
                ├── iso.sorted.bam      # isolate-specific
                └── anc.sorted.bam      # copied from ancestor run

# Per-candidate folder naming convention:
#   jc_{start}_{end}_{tn_id}_L{read_len}
#   - start, end: amplicon coordinates (genomic positions)
#   - tn_id: TN index from TnLoc table (e.g., 001, 002)
#   - read_len: read length (affects junction flank size)
#   Example: jc_123456_234567_001_L150
#
# Caching logic:
#   1. Ancestor run creates (in {anc_name}/{anc_name}/jc_.../):
#      - junctions.fasta
#      - iso.sorted.bam (ancestor reads aligned to junctions)
#   2. Isolate runs copy from ancestor's jc_.../ folder:
#      - junctions.fasta → junctions.fasta
#      - iso.sorted.bam → anc.sorted.bam (renamed: ancestor's "iso" is isolate's "anc")
#   3. Isolate creates its own iso.sorted.bam
#
# Coverage:
#   - Ancestor run (anc_path=None): RAW coverage (no normalization)
#   - Isolate run (anc_path=set): iso/anc RATIO

# Isolate without ancestor (RAW coverage - not normalized)
{output}/
└── {ref_name}/
    └── {iso_name}/                 # iso is its own "ancestor" for folder structure
        └── {iso_name}/             # uses RAW coverage (anc_path=None)
            ├── breseq/
            ├── tn_jc.csv
            ├── ...
            └── jc_{...}/           # per-candidate (no anc.sorted.bam needed)
                ├── junctions.fasta
                └── iso.sorted.bam
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
- Output organized by: `{output}/{ref_name}/{anc_name}/{iso_name}/`
- `ref_name` - reference genome (required), top-level grouping
- `iso_name` - sample name (required), the sample being analyzed
- `anc_name` - ancestor name (defaults to `iso_name` when `anc_path` is None)
- `anc_path` - determines coverage normalization (None → raw, set → normalized)

```bash
# Run isolate without ancestor (raw coverage, no normalization)
amplifinder -i my_isolate.fastq --ref U00096 --iso-name my_isolate
# Output: output/U00096/my_isolate/my_isolate/   (anc_name defaults to iso_name)
# Coverage: RAW (not normalized)

# Run ancestor standalone (for later use by isolates)
amplifinder -i my_ancestor.fastq --ref U00096 --iso-name my_ancestor
# Output: output/U00096/my_ancestor/my_ancestor/
# Coverage: RAW (same as isolate without ancestor)

# Run isolate with ancestor (normalized coverage)
amplifinder -i my_isolate.fastq --ref U00096 --iso-name my_isolate -a my_ancestor.fastq --anc-name my_ancestor
# Output: output/U00096/my_ancestor/my_isolate/
# Coverage: NORMALIZED (iso/anc ratio)
# Auto-runs ancestor to output/U00096/my_ancestor/my_ancestor/ if not already done

# Run another isolate with same ancestor (reuses cached ancestor run + anc.sorted.bam)
amplifinder -i my_isolate_2.fastq --ref U00096 --iso-name my_isolate_2 -a my_ancestor.fastq --anc-name my_ancestor
# Output: output/U00096/my_ancestor/my_isolate_2/
# Skips ancestor run (already exists at output/U00096/my_ancestor/my_ancestor/)
# For shared candidates: copies junctions.fasta + anc.sorted.bam from my_isolate_1
```

### Convergence Points (anc_path=`None` vs anc_path=`{path}`)

```
Step 7:  anc_path=None → raw coverage only
         anc_path=set  → load anc coverage, compute amplicon_coverage ratio

Steps 10-11: anc_path=None → create junctions, align THIS sample only
             anc_path=set  → create junctions, align BOTH iso reads AND anc reads

Step 12: anc_path=None → parse own BAM → jc_cov_* fields  
         anc_path=set  → parse BOTH BAMs → iso_jc_cov_* + anc_jc_cov_* fields

Step 13: anc_path=None → limited classification
         anc_path=set  → full iso vs anc pattern comparison
```

---

## Implementation Modules

> **IMPORTANT: Column Naming Consistency**
> 
> The original MATLAB code has inconsistent column naming where the same data appears with different 
> names in different output tables (e.g., `pos_Chr` vs `Positions_in_chromosome`, `amplicon_coverage` 
> vs `median_copy_number`). 
> 
> **Best practice for Python implementation:**
> 1. Use **record inheritance** (`TnJc2 → CoveredTnJc2 → ClassifiedTnJc2 → ...`) to ensure column 
>    names are consistent across all downstream tables
> 2. Define column names **once** in the base record and inherit them - never rename inherited fields
> 3. Only add **new** fields in derived records, never duplicate existing fields with different names
> 4. Apply renaming only at **export time** (Step 14) if user-friendly names are needed for Excel output

> **IMPORTANT: Include Reference and Ancestor Context**
> 
> In the MATLAB code, `Reference` and `Ancestor` columns are only added at export time, making 
> intermediate tables lack essential context. 
> 
> **Fix in Python implementation:**
> 1. Add `ref_name: str` to `TnJc2` (base record) - all junctions are relative to a reference
> 2. Add `iso_name: str` to `TnJc2` - identifies which sample the junctions came from
> 3. Add `anc_name: Optional[str]` to `CoveredTnJc2` - present when ancestor comparison is performed
> 
> This ensures every output table is self-documenting and traceable to its source data.
> 
> **Action required:** Update existing `TnJc2` in `data_types/record_types.py` to include:
> ```python
> class TnJc2(Record):
>     # Context fields (add these)
>     ref_name: str      # reference genome name
>     iso_name: str      # isolate/sample name
>     
>     # ... existing fields ...
> ```
> 
> And update `CoveredTnJc2` definition (Step 7) to include:
> ```python
> class CoveredTnJc2(TnJc2):
>     anc_name: Optional[str] = None  # ancestor name (None if standalone run)
>     # ... existing coverage fields ...
> ```

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
      # Context (inherited: ref_name, iso_name from TnJc2)
      anc_name: Optional[str] = None    # ancestor name (None if standalone run)
      
      # Coverage fields
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

> **NOTE: `raw_event` vs `isolate_architecture` (Step 13)**
> 
> Both fields use similar vocabulary ("flanked", "unflanked", etc.) but represent **different analyses**:
> - `raw_event` (Step 8): **Structural prediction** based on junction pair relationships (which junctions share single-locus IS elements)
> - `isolate_architecture` (Step 13): **Empirical verification** based on read alignment patterns to synthetic junctions
> 
> **Keep both fields** - comparing them can reveal discrepancies between predicted and observed architectures, which may indicate artifacts or complex rearrangements.

### 4. FilterCandidatesStep (Step 9)
- Module: `steps/filter_candidates.py`
- MATLAB: filtering logic in `curate_candidate_amplicons.m`
- **Record types:**
  ```python
  class CandidateTnJc2(ClassifiedTnJc2):
      analysis_directory: str  # "jc_{start}_{end}_{tn_id}_L{read_len}"
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
      event: Tuple[RawEvent, List[EventModifierAndSide]]
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


  class EventModifierAndSide(NamedTuple):
      event_modifier: EventModifier
      side: Side
  ```
- **input_files:** None
- **data_inputs:**
  - `analyzed_candidates`: RecordTypedDF[AnalyzedCandidateTnJc2]
- **Output:** RecordTypedDF[AnalyzedCandidateTnJc2] (with `event`, `isolate_architecture`, `ancestor_architecture` populated)
- File output: None (memory only)
- Logic: Match junction read patterns to expected architectures
- Tests: `test_steps/test_classify_candidates.py`
- Commit: "Implement ClassifyCandidatesStep"

> **NOTE: See Step 8 note on `raw_event` vs `isolate_architecture`**
> 
> `isolate_architecture` is derived from read alignment evidence, while `raw_event` (inherited from Step 8) 
> is derived from junction pairing logic. Both should be retained for comparison.

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
  - If `anc_path` set:
    - Check if ancestor run exists at `{output}/{anc_name}/{anc_name}/`
    - If missing → recursively run full pipeline with `iso_name=anc_name, anc_path=None`
    - Load ancestor config via `load_config_from_run()`
  - Run full pipeline for isolate
  - Pass ancestor results to convergence points (Steps 7, 10-12)
- Tests: `test_integration/test_pipeline.py`
- Commit: "Integrate remaining steps with ancestor support"

### 13. Full Integration Tests
- Tests: `test_integration/test_full_run.py`
- Commit: "Add full integration tests"


### 14. VisualizeCandidates (Step 15)
This is not a `Step`. It is a module that allows displaying the results of the pipeline after the fact. 
- Module: `visualization/visualize.py`
- Library: `matplotlib` (v3.7+) [need to pip install in correct conda env]
- **input_files:**
  - path to run folder
- **Output:** None (side effect: displays figures)
- File output (optional): `{analysis_dir}/coverage_plot.png`
- **User interaction modes:**
  - `amplifinder --visualize run_folder`: Show all candidates sequentially (with buttons on the figure to go NEXT and BACK)
  - `amplifinder --visualize run_folder`: Save plots to png files, instead of interactive display
- **Plot specification:**
      
      Layout:
      - Title: "{raw_event} | {amplicon_start}-{amplicon_end} | copy_number: {copy_number:.1f}x"
      - X-axis: genomic position
      - Y-axis: coverage depth
      - Blue line: isolate coverage
      - Gray line: ancestor coverage (if available)
      - Vertical dashed red lines: amplicon start/end positions
      - Shaded region: amplicon interior
      
      Plot range: [amplicon_start - flank, amplicon_end + flank]
      where flank = amplicon_length * flank_fraction

- Tests: `test_steps/test_visualize.py`
- Commit: "Implement VisualizeCandidates"

### "THE END"

- Complete your work witht he interactive GUI running - to show off!

---

## Development Workflow

For each implementation module, follow this scheme:
A. implement module
B. write unit tests - run/debug until passing
C. write integration test with Matlab comparison (test_integration, AmpliFinder_test) - run/debug
D. flake8 - run and fix
E. git commit
F. write a summary to the chat with: 
    (a) where we are in the overall process
    (b) what you did in this module
    (c) number of new code lines for this module/commit.
G. Continue to the next development module


## Test Patterns

Follow existing patterns in `tests/`:

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

Integration tests compare Python outputs against MATLAB reference outputs from `AmpliFinder_test/`:

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
    assert python_df["amplicon_end"].tolist() == matlab_expected["amplicon_end"]
```

Reference data location: `AmpliFinder_test/` (MATLAB outputs for known test cases)

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
- **matplotlib** (v3.7+) - for coverage visualization [need to `pip install`, when needed]
- **scipy** - for reading MATLAB files (`scipy.io.loadmat()`, `pip install` as needed)


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

# Ancestor run's candidate dir (CANONICAL SOURCE)
anc_jc_dir = anc_run_dir / jc_name
anc_jc_dir / "junctions.fasta"           # created during ancestor run
anc_jc_dir / "iso.sorted.bam"            # ancestor's alignment (becomes isolate's anc.sorted.bam)

# Isolate run's candidate dir
iso_jc_dir = iso_run_dir / jc_name
iso_jc_dir / "junctions.fasta"           # copied from anc_jc_dir
iso_jc_dir / "anc.sorted.bam"            # copied from anc_jc_dir/iso.sorted.bam (renamed)
iso_jc_dir / "iso.sorted.bam"            # isolate-specific alignment
iso_jc_dir / "coverage_plot.png"         # (if --save-plots)
```


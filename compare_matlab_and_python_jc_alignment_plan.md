# Plan: Compare MATLAB and Python Junction Alignment Results

## Overview
This plan details a systematic comparison between MATLAB and Python pipeline outputs for junction alignment analysis, focusing on identifying differences in BAM file alignments and read-type classification (left/right/green) for the junction `chr_U00096_873161F_890745R_IS_41R`.

## Target Junction
- **Junction ID**: `chr_U00096_873161F_890745R_IS_41R`
- **MATLAB BAM**: `/zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877/chr_U00096_873161F_890745R_IS_41R/alignment/alignment.sorted.bam`
- **Python BAM**: `/zdata/user-data/idan/small_projects/AmpliFinder/pyAmpliFinder/amplifinder/tests/test_outputs/integration/output/U00096/SRR25242906/SRR25242877/junctions/jc_U00096_873161-890745_U00096_3583428-3584195_S_302bp/sorted.bam`

---

## Part 1: BAM File Comparison

### 1.1 Convert BAM to SAM for Text Comparison

**Objective**: Convert both BAM files to SAM format for human-readable comparison, preserving all flags and tags.

**Steps**:
1. Use `samtools view` to convert MATLAB BAM to SAM (with all tags):
   ```bash
   samtools view -h <matlab_bam> > matlab_alignment.sam
   ```
   The `-h` flag includes the header. By default, `samtools view` outputs all tags.

2. Use `samtools view` to convert Python BAM to SAM (with all tags):
   ```bash
   samtools view -h <python_bam> > python_alignment.sam
   ```

3. **Verify completeness**: Check that SAM files contain all tags by:
   - Inspecting a few lines to ensure tags are present
   - Using `grep` to count lines with tags (should be all alignment lines)
   - Comparing tag presence between MATLAB and Python SAM files

4. Store SAM files in a comparison output directory.

**Note**: SAM format includes all tags in the last column (TAB-separated TAG:TYPE:VALUE format). All tags present in the BAM file will be included in the SAM output.

**Output**: 
- `matlab_alignment.sam` (with all flags and tags)
- `python_alignment.sam` (with all flags and tags)

### 1.2 Extract Alignment Records

**Objective**: Parse both BAM files and extract all alignment records with ALL fields and tags, including those not present for each read.

**All SAM Fields to Extract** (include all, use "N/A" or empty when not present):
- **Standard SAM fields**:
  - `QNAME`: Read ID
  - `FLAG`: SAM flag (integer)
  - `RNAME`: Reference name
  - `POS`: 1-based leftmost mapping position
  - `MAPQ`: Mapping quality
  - `CIGAR`: CIGAR string
  - `RNEXT`: Reference name of mate/next read
  - `PNEXT`: Position of mate/next read
  - `TLEN`: Template length
  - `SEQ`: Read sequence
  - `QUAL`: Quality scores (ASCII)

- **Flag-derived boolean fields** (extracted from FLAG):
  - `is_paired`: Read is paired
  - `is_proper_pair`: Read is in proper pair
  - `is_unmapped`: Read is unmapped
  - `mate_unmapped`: Mate is unmapped
  - `is_reverse`: Read is reverse strand
  - `mate_reverse`: Mate is reverse strand
  - `is_read1`: First read in pair
  - `is_read2`: Second read in pair
  - `is_secondary`: Secondary alignment
  - `is_qcfail`: QC failure
  - `is_duplicate`: PCR or optical duplicate
  - `is_supplementary`: Supplementary alignment

- **Coordinate fields** (for compatibility):
  - `ref_start`: 0-based start (for Python compatibility)
  - `ref_end`: 0-based end (exclusive)
  - `query_alignment_start`: Query alignment start (0-based)
  - `query_alignment_end`: Query alignment end (0-based, exclusive)
  - `query_length`: Full query length

- **All possible SAM tags** (extract all tags present in BAM, include columns for all common tags even if not present):
  - `AS`: Alignment score (i)
  - `NM`: Edit distance (i)
  - `nM`: Number of mismatches (i)
  - `NH`: Number of reported alignments (i)
  - `HI`: Query hit index (i)
  - `XS`: Suboptimal alignment score (i)
  - `MD`: String for mismatching positions (Z)
  - `MC`: CIGAR string for mate/next segment (Z)
  - `MQ`: Mapping quality of mate/next segment (i)
  - `SA`: Other canonical alignments (Z)
  - `RG`: Read group (Z)
  - `BC`: Barcode sequence (Z)
  - `OC`: Original CIGAR (Z)
  - `OP`: Original mapping position (i)
  - `OA`: Original alignment (Z)
  - Any other tags found in the BAM files

**Implementation**:
- Use `pysam` for Python BAM parsing
- For MATLAB BAM, also use `pysam` (BAM format is standard)
- For each alignment:
  - Extract all standard SAM fields
  - Parse FLAG to extract all boolean flag fields
  - Extract ALL tags present in the alignment using `read.get_tags()`
  - Create a comprehensive record with all fields
  - For tags not present, use "N/A" or empty string
- Store all alignments in structured format (e.g., pandas DataFrame) with columns for ALL possible fields and tags

**Note**: When extracting tags, iterate through all tags in each read using `read.get_tags()` to ensure no tags are missed. Create a comprehensive list of all unique tags found across all reads, and include columns for all of them in the output.

### 1.3 Compare Alignments Read-by-Read

**Objective**: Identify differences in alignments between MATLAB and Python BAM files.

**Comparison Strategy**:
1. **Group by Read ID**: Collect all alignments for each read ID from both BAM files
2. **Match Alignments**: For each read, attempt to match alignments between MATLAB and Python based on:
   - Reference name
   - Start position (with small tolerance, e.g., ±1 bp)
   - CIGAR string
   - Flag (primary/secondary/supplementary status)
3. **Identify Differences**:
   - Reads present in MATLAB but not in Python
   - Reads present in Python but not in MATLAB
   - Reads with different alignment positions
   - Reads with different CIGAR strings
   - Reads with different flags (secondary/supplementary status)
   - Reads with different alignment scores (AS, NM tags)
   - Reads with different number of hits (NH tag)

**Output Categories**:
- **Matched alignments**: Alignments that appear in both with identical key fields
- **Unmatched MATLAB alignments**: Alignments only in MATLAB BAM
- **Unmatched Python alignments**: Alignments only in Python BAM
- **Differing alignments**: Same read, but different alignment properties

### 1.4 Generate Comparison Outputs

#### 1.4.1 Summary Report

**Format**: Markdown or text file

**Contents**:
- Total number of reads in MATLAB BAM
- Total number of reads in Python BAM
- Number of unique read IDs in MATLAB
- Number of unique read IDs in Python
- Number of matched alignments
- Number of unmatched MATLAB alignments
- Number of unmatched Python alignments
- Number of reads with differing alignments
- Breakdown by alignment type (primary/secondary/supplementary)
- Breakdown by reference (junction type: 1-7)

#### 1.4.2 FASTQ File with Differing Reads

**Objective**: Extract read sequences for reads that have different mappings.

**Steps**:
1. Identify all read IDs with differing alignments
2. Extract original read sequences from source FASTQ files (if available) or from BAM SEQ field
3. Create FASTQ file with:
   - Read ID as header
   - Sequence
   - Quality scores

**Output**: `green_mismatch_reads.fastq` (or similar naming)

**Note**: If original FASTQ is not available, extract from BAM SEQ/QUAL fields, but note that BAM may contain reverse-complemented sequences for reverse-strand alignments.

#### 1.4.3 CSV File with All Differing Hits

**Objective**: Create detailed CSV listing all alignment records that differ between MATLAB and Python, including ALL flags and annotations.

**Format**: CSV with columns for ALL SAM fields and tags (use "N/A" or empty when not present):

**Standard SAM Fields**:
- `read_id` (QNAME): Read identifier
- `source`: "matlab" or "python"
- `flag`: SAM flag (integer)
- `ref_name` (RNAME): Reference name (junction type)
- `pos` (POS): 1-based position
- `mapq` (MAPQ): Mapping quality
- `cigar` (CIGAR): CIGAR string
- `rnext` (RNEXT): Reference name of mate/next read
- `pnext` (PNEXT): Position of mate/next read
- `tlen` (TLEN): Template length
- `sequence` (SEQ): Read sequence (optional, may be large)
- `quality` (QUAL): Quality scores (optional, may be large)

**Flag-Derived Boolean Fields** (all extracted from FLAG):
- `is_paired`
- `is_proper_pair`
- `is_unmapped`
- `mate_unmapped`
- `is_reverse`
- `mate_reverse`
- `is_read1`
- `is_read2`
- `is_secondary`
- `is_qcfail`
- `is_duplicate`
- `is_supplementary`

**Coordinate Fields**:
- `ref_start`: 0-based start (for Python compatibility)
- `ref_end`: 0-based end (exclusive)
- `query_alignment_start`: Query alignment start (0-based)
- `query_alignment_end`: Query alignment end (0-based, exclusive)
- `query_length`: Full query length

**All SAM Tags** (include columns for all tags found in either BAM file, use "N/A" when not present):
- `AS`: Alignment score (i)
- `NM`: Edit distance (i)
- `nM`: Number of mismatches (i)
- `NH`: Number of reported alignments (i)
- `HI`: Query hit index (i)
- `XS`: Suboptimal alignment score (i)
- `MD`: String for mismatching positions (Z)
- `MC`: CIGAR string for mate/next segment (Z)
- `MQ`: Mapping quality of mate/next segment (i)
- `SA`: Other canonical alignments (Z)
- `RG`: Read group (Z)
- `BC`: Barcode sequence (Z)
- `OC`: Original CIGAR (Z)
- `OP`: Original mapping position (i)
- `OA`: Original alignment (Z)
- Plus any other tags found in the BAM files (dynamically add columns as needed)

**Implementation Notes**:
1. **First pass**: Scan both BAM files to collect all unique tag names
2. **Create column list**: Include all standard SAM fields, all flag-derived fields, and all unique tags found
3. **For each alignment**: Populate all columns, using "N/A" for missing tags/fields
4. **Tag extraction**: Use `read.get_tags()` to get all tags, then populate corresponding columns

**Sorting**:
- Primary sort: `read_id` (alphabetical)
- Secondary sort: `source` (matlab first, then python)
- Tertiary sort: `ref_name`, then `pos`

**Grouping**: For each read ID, list all MATLAB hits first, then all Python hits, to facilitate side-by-side comparison.

**Output**: `bam_comparison_differences.csv`

---

## Part 2: Read-Type Counting Comparison

### 2.1 Extract MATLAB Read Counts

**Objective**: Extract read counts (left, right, green) from MATLAB output files.

**Data Sources**:
- **bamreads.mat**: Contains `nmbr_green_reads`, `nmbr_left_reads`, `nmbr_right_reads` (arrays of size 1x7, one per junction type)
- **biomap.mat**: Contains `RdStart`, `RdRef`, `RdLength`, `RdFlag`, `RdFull` (per-read alignment data)

**Location**: Same directory as MATLAB BAM file:
```
/zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877/chr_U00096_873161F_890745R_IS_41R/alignment/
```

**Steps**:
1. **Load MATLAB CSV files** (already exported in format `<matfile>__<varname>.csv`):
   - Load `bamreads__nmbr_green_reads.csv` to get summary counts per junction type
   - Load `bamreads__nmbr_left_reads.csv`
   - Load `bamreads__nmbr_right_reads.csv`
   - Load `biomap__RdStart.csv` - Start positions (1-based, indexed by read position in BAM)
   - Load `biomap__RdRef.csv` - Reference names (junction types, as strings '1'..'7')
   - Load `biomap__RdLength.csv` - Alignment lengths
   - Load `biomap__RdFlag.csv` - SAM flags
   - Load `biomap__RdFull.csv` - Boolean array (true if CIGAR-only-M)

**Note**: MATLAB variables are already exported to CSV files in the format `<matfile>__<varname>.csv`. These files should be located in the same directory as the MATLAB BAM file.

**Output Files** (to be loaded):
- `bamreads__nmbr_green_reads.csv`
- `bamreads__nmbr_left_reads.csv`
- `bamreads__nmbr_right_reads.csv`
- `biomap__RdStart.csv`
- `biomap__RdRef.csv`
- `biomap__RdLength.csv`
- `biomap__RdFlag.csv`
- `biomap__RdFull.csv`

### 2.2 Extract Python Read Counts

**Objective**: Extract or generate read counts (left, right, green) from Python pipeline.

**Data Sources**:
- Check if Python pipeline saves read counts to files (likely in `AnalyzedTnJc2` records or similar)
- If not saved, regenerate counts by:
  1. Parsing Python BAM file
  2. Applying Python classification logic (from `analyze_alignments.py` and `JunctionReadCounts`)
  3. Counting reads per junction type and category

**Python Classification Logic** (from codebase):
- Filter: `read.is_unmapped`, `read.is_secondary`, `read.is_supplementary` → skip
- Filter: Alignment length within `avg_read_length * (1 ± tolerance)` (default tolerance 0.1)
- Filter: CIGAR-only-M (`RdFull` equivalent)
- Filter: `NM <= 3` and `AS >= -25`
- Classification uses `JunctionReadCounts.add_read()` which calls `get_read_type()`:
  - **LEFT**: Read starts on left arm, ends before junction
  - **RIGHT**: Read starts after junction, ends on right arm
  - **MIDDLE (green)**: Read starts before junction and ends after junction (with min_overlap_len on both sides)

**Steps**:
1. Check for existing count files in Python output directory
2. If missing, create script to:
   - Load Python BAM file
   - For each junction type (1-7), apply classification
   - Count reads per category (left, right, green/spanning)
   - Save to CSV: `python_read_counts.csv` with columns: `junction_type`, `green_count`, `left_count`, `right_count`

**Output Files**:
- `python_read_counts.csv` (summary)
- `python_per_read_classification.csv` (optional, detailed per-read data)

### 2.3 Compare Summary Counts

**Objective**: Create side-by-side comparison of read counts.

**Format**: CSV or Markdown table

**Columns**:
- `junction_type`: Junction type (1-7)
- `matlab_green`: Count from MATLAB
- `python_green`: Count from Python
- `green_diff`: Difference (Python - MATLAB)
- `matlab_left`: Count from MATLAB
- `python_left`: Count from Python
- `left_diff`: Difference
- `matlab_right`: Count from MATLAB
- `python_right`: Count from Python
- `right_diff`: Difference

**Output**: `read_count_comparison_summary.csv`

### 2.4 Identify Differentially Labelled Reads

**Objective**: Find specific reads that are classified differently between MATLAB and Python.

**Challenge**: MATLAB stores alignment data by index (position in BAM file), not by read ID. We need to match MATLAB indices to actual read IDs from the BAM file, then compare classifications.

#### Step 2.4.1: Build MATLAB Index-to-ReadID Mapping

**Objective**: Create a mapping from MATLAB array indices to read IDs from the MATLAB BAM file.

**Steps**:
1. **Parse MATLAB BAM file**:
   - Open MATLAB BAM file using `pysam.AlignmentFile`
   - **First pass**: Collect all unique tag names across all reads
   - Iterate through all alignments in order (preserving BAM file order)
   - For each alignment at index `i` (0-based):
     - Extract all standard SAM fields: `read_id` (QNAME), `ref_name`, `start` (1-based), `end` (1-based), `cigar`, `flag`, `mapq`, `rnext`, `pnext`, `tlen`, `seq`, `qual`
     - Extract all flag-derived boolean fields (is_paired, is_secondary, etc.)
     - Extract ALL tags using `read.get_tags()` - store as dict with tag name as key
     - For tags not present, use "N/A" or None
     - Store in mapping: `matlab_index_to_readid[i] = read_id`
     - Store alignment details: `matlab_alignments[i] = {
         'read_id': read_id,
         'ref_name': ref_name,
         'start': start,
         'end': end,
         'cigar': cigar,
         'flag': flag,
         'mapq': mapq,
         'rnext': rnext,
         'pnext': pnext,
         'tlen': tlen,
         'seq': seq,
         'qual': qual,
         'is_paired': read.is_paired,
         'is_proper_pair': read.is_proper_pair,
         'is_unmapped': read.is_unmapped,
         'mate_unmapped': read.mate_unmapped,
         'is_reverse': read.is_reverse,
         'mate_reverse': read.mate_reverse,
         'is_read1': read.is_read1,
         'is_read2': read.is_read2,
         'is_secondary': read.is_secondary,
         'is_qcfail': read.is_qcfail,
         'is_duplicate': read.is_duplicate,
         'is_supplementary': read.is_supplementary,
         'tags': {tag: value for tag, value in read.get_tags()}  # All tags
     }`

2. **Create reverse mapping** (read_id + alignment properties → MATLAB index):
   - For each MATLAB index `i`:
     - Get read_id and alignment properties
     - Create key: `(read_id, ref_name, start, end, cigar, flag)`
     - Store: `matlab_alignment_key_to_index[key] = i`

**Output Data Structures**:
- `matlab_index_to_readid`: dict[int] → str (read_id)
- `matlab_alignments`: dict[int] → dict (full alignment details)
- `matlab_alignment_key_to_index`: dict[tuple] → int (for reverse lookup)

#### Step 2.4.2: Extract MATLAB Per-Read Classifications

**Objective**: Apply MATLAB classification logic to each read and store results.

**Steps**:
1. **Load MATLAB data from CSV files**:
   - Load `biomap__RdStart.csv` → array `RdStart` (1-based positions)
   - Load `biomap__RdRef.csv` → array `RdRef` (junction types as strings '1'..'7')
   - Load `biomap__RdLength.csv` → array `RdLength` (alignment lengths)
   - Load `biomap__RdFlag.csv` → array `RdFlag` (SAM flags)
   - Load `biomap__RdFull.csv` → array `RdFull` (boolean, CIGAR-only-M)

2. **Get classification parameters**:
   - Junction length (`seq_length`): From BAM header or FASTA file
   - Read length: From `fastq_read_length.mat` (exported as CSV) or config (typically 150)
   - `req_ovrlp` (min_overlap_len): From config or default (typically 12)
   - `min_bp_in_frame`: 10 (hardcoded in MATLAB)

3. **Apply MATLAB classification logic for each index `i`**:
   - Skip if `RdFull[i] == False` (not CIGAR-only-M)
   - Check read length: `okz_length = (RdLength[i] > read_length * 0.9) and (RdLength[i] < read_length * 1.1)`
   - If not `okz_length`, skip (classification = None)
   - Get junction type: `jct_type = RdRef[i]` (string '1'..'7')
   - Calculate junction midpoint: `jct_point = seq_length / 2`
   - Calculate read end (1-based inclusive): `read_end = RdStart[i] + RdLength[i] - 1`
   
   - **Green classification**:
     ```python
     okz_green = (RdStart[i] < jct_point - req_ovrlp) and (read_end > jct_point + req_ovrlp)
     ```
   
   - **Left classification**:
     ```python
     okz_left = (RdStart[i] + RdLength[i] > jct_point - read_length + min_bp_in_frame) and \
                (read_end < jct_point + req_ovrlp)
     ```
   
   - **Right classification**:
     ```python
     okz_right = (RdStart[i] > jct_point - req_ovrlp) and \
                 (RdStart[i] < jct_point + read_length - min_bp_in_frame)
     ```

4. **Store MATLAB classifications**:
   - For each index `i`:
     - Get read_id: `read_id = matlab_index_to_readid.get(i)`
     - If read_id is None, skip (alignment not in BAM or index mismatch)
     - Determine classification:
       - If `okz_green`: classification = "green"
       - Else if `okz_left`: classification = "left"
       - Else if `okz_right`: classification = "right"
       - Else: classification = None
     - Store: `matlab_classifications[(read_id, jct_type)] = classification`
     - Store details: `matlab_classification_details[(read_id, jct_type)] = {
         'classification': classification,
         'start': RdStart[i],
         'end': read_end,
         'length': RdLength[i],
         'flag': RdFlag[i],
         'matlab_index': i
     }`

**Output Data Structures**:
- `matlab_classifications`: dict[(read_id, jct_type)] → str (classification: "green"/"left"/"right"/None)
- `matlab_classification_details`: dict[(read_id, jct_type)] → dict (full details)

#### Step 2.4.3: Extract Python Per-Read Classifications

**Objective**: Parse Python BAM and apply Python classification logic.

**Steps**:
1. **Parse Python BAM file**:
   - Open Python BAM file using `pysam.AlignmentFile`
   - Get junction lengths from BAM header: `jct_lengths = dict(zip(bam.references, bam.lengths))`

2. **Get classification parameters**:
   - `avg_read_length`: From config or calculated (typically 150)
   - `min_overlap_len`: From config (typically 12)
   - `read_length_tolerance`: 0.1 (default)
   - `max_dist_from_junction`: 10 (default, equivalent to min_bp_in_frame)

3. **For each alignment in Python BAM**:
   - Skip if `read.is_unmapped`, `read.is_secondary`, or `read.is_supplementary`
   - Get `read_id = read.query_name`
   - Get `jct_name = read.reference_name`
   - Get `jct_type = JunctionType[jct_name]` (convert to string '1'..'7' for comparison)
   - Get `jct_len = jct_lengths[jct_name]`
   
   - **Apply Python filters**:
     - Alignment length: `alignment_length = read.query_alignment_length`
     - Check: `min_alignment_length <= alignment_length <= max_alignment_length`
     - Check CIGAR-only-M: `all(op == 0 for op, _ in read.cigartuples)`
     - Check alignment quality tags: `NM <= 3` and `AS >= -25`
     - If any filter fails, skip (classification = None)
   
   - **Apply Python classification**:
     - Convert to 1-based coordinates: `start_1 = read.reference_start + 1`, `end_1 = read.reference_end`
     - Calculate `arm_len = jct_len // 2`
     - Use `JunctionReadCounts.get_read_type()` logic:
       - Calculate indices: `idx_L_1 = arm_len - end_1`, `idx_L_2 = arm_len - start_1`
       - Calculate indices: `idx_R_1 = start_1 - (arm_len + 1)`, `idx_R_2 = end_1 - (arm_len + 1)`
       - Check distance filters
       - Determine classification:
         - If `idx_L_1 >= 0`: classification = "left"
         - Else if `idx_R_1 >= 0`: classification = "right"
         - Else if spanning both sides: classification = "green" (MIDDLE)
         - Else: classification = None
   
   - **Store Python classifications**:
     - Extract ALL tags using `read.get_tags()` - store as dict
     - Store: `python_classifications[(read_id, jct_type)] = classification`
     - Store details: `python_classification_details[(read_id, jct_type)] = {
         'classification': classification,
         'start': start_1,
         'end': end_1,
         'length': alignment_length,
         'cigar': read.cigarstring,
         'flag': read.flag,
         'mapq': read.mapping_quality,
         'rnext': read.next_reference_name,
         'pnext': read.next_reference_start,
         'tlen': read.template_length,
         'seq': read.query_sequence,
         'qual': read.query_qualities,
         'is_paired': read.is_paired,
         'is_proper_pair': read.is_proper_pair,
         'is_unmapped': read.is_unmapped,
         'mate_unmapped': read.mate_unmapped,
         'is_reverse': read.is_reverse,
         'mate_reverse': read.mate_reverse,
         'is_read1': read.is_read1,
         'is_read2': read.is_read2,
         'is_secondary': read.is_secondary,
         'is_qcfail': read.is_qcfail,
         'is_duplicate': read.is_duplicate,
         'is_supplementary': read.is_supplementary,
         'tags': {tag: value for tag, value in read.get_tags()}  # All tags
     }`

**Output Data Structures**:
- `python_classifications`: dict[(read_id, jct_type)] → str
- `python_classification_details`: dict[(read_id, jct_type)] → dict

#### Step 2.4.4: Match and Compare Classifications

**Objective**: Match reads between MATLAB and Python and identify classification differences.

**Steps**:
1. **Find common read-junction pairs**:
   - Get all keys from both `matlab_classifications` and `python_classifications`
   - Find intersection: `common_keys = set(matlab_classifications.keys()) & set(python_classifications.keys())`
   - Find MATLAB-only: `matlab_only = set(matlab_classifications.keys()) - set(python_classifications.keys())`
   - Find Python-only: `python_only = set(python_classifications.keys()) - set(matlab_classifications.keys())`

2. **Compare classifications for common reads**:
   - For each `(read_id, jct_type)` in `common_keys`:
     - Get MATLAB classification: `matlab_class = matlab_classifications[(read_id, jct_type)]`
     - Get Python classification: `python_class = python_classifications[(read_id, jct_type)]`
     - If `matlab_class != python_class`, mark as mismatch
     - Get details from both sources

3. **Handle MATLAB-only reads**:
   - For each `(read_id, jct_type)` in `matlab_only`:
     - MATLAB has classification, Python does not
     - Mark as mismatch

4. **Handle Python-only reads**:
   - For each `(read_id, jct_type)` in `python_only`:
     - Python has classification, MATLAB does not
     - Mark as mismatch

5. **Generate output CSV**:
   - Collect all unique tag names from both MATLAB and Python alignment details
   - For each mismatch, create row with ALL fields:
     - **Basic fields**:
       - `read_id`
       - `junction_type`
       - `matlab_classification` (or "N/A" if not in MATLAB)
       - `python_classification` (or "N/A" if not in Python)
     - **MATLAB alignment fields** (or "N/A" if not in MATLAB):
       - `matlab_start`, `matlab_end`, `matlab_length`
       - `matlab_cigar`, `matlab_flag`, `matlab_mapq`
       - `matlab_rnext`, `matlab_pnext`, `matlab_tlen`
       - All MATLAB flag-derived booleans (is_paired, is_secondary, etc.)
       - All MATLAB tags (AS, NM, nM, NH, HI, XS, MD, MC, MQ, SA, RG, BC, OC, OP, OA, etc.)
     - **Python alignment fields** (or "N/A" if not in Python):
       - `python_start`, `python_end`, `python_length`
       - `python_cigar`, `python_flag`, `python_mapq`
       - `python_rnext`, `python_pnext`, `python_tlen`
       - All Python flag-derived booleans
       - All Python tags (AS, NM, nM, NH, HI, XS, MD, MC, MQ, SA, RG, BC, OC, OP, OA, etc.)
     - `difference_reason`: Brief explanation (e.g., "MATLAB: green, Python: left", "Only in MATLAB", etc.)
   
   **Note**: Include columns for ALL tags found in either BAM file, using "N/A" when a tag is not present for a particular read.

**Output**: 
- `differentially_labelled_reads.csv` with columns as specified above
- Sorted by `read_id`, then `junction_type`

**Validation**:
- Verify that MATLAB indices align correctly with BAM file order
- Check for edge cases where multiple alignments of same read might cause index mismatches
- Handle cases where MATLAB has more/fewer alignments than Python BAM

### 2.5 Create Per-Junction Read Lists

**Objective**: For each junction type, list all reads counted in each category (left, right, green) for both MATLAB and Python.

**Format**: Separate CSV files per junction type, or single CSV with junction_type column.

**For MATLAB**:
- Use `bamreads.mat` to get read indices for green/left/right reads per junction
- Cross-reference with `biomap.mat` to get read details
- Extract read IDs from BAM file (may require parsing BAM Signature field or matching by position/CIGAR)

**For Python**:
- Parse BAM file and apply classification
- Store read IDs per category

**Output Structure**:
- `matlab_reads_junction_<N>.csv` (N=1..7) with columns:
  - `read_id`
  - `classification` (green/left/right)
  - `start`
  - `end`
  - `length`
  - `cigar`
  - `flag`
- `python_reads_junction_<N>.csv` with same structure

**Alternative**: Single file `per_junction_read_lists.csv` with columns:
- `junction_type`
- `source` (matlab/python)
- `read_id`
- `classification`
- `start`
- `end`
- `length`
- `cigar`
- `flag`

---

## Implementation Details

### Tools and Libraries Required

1. **Python**:
   - `pysam`: For BAM/SAM file parsing
   - `pandas`: For data manipulation and CSV export
   - `scipy.io` or `h5py`: For loading MATLAB .mat files
   - `samtools`: Command-line tool (via subprocess or direct binary)

2. **MATLAB** (if needed for .mat export):
   - MATLAB or Octave for exporting .mat to CSV (if CSVs don't exist)

### Script Structure

**Suggested script organization**:
```
compare_jc_alignment/
├── compare_bam_files.py          # Part 1: BAM comparison
├── compare_read_counts.py          # Part 2: Read count comparison
├── extract_matlab_data.py         # Helper: Extract MATLAB .mat to CSV
├── extract_python_data.py         # Helper: Extract Python classifications
├── utils/
│   ├── bam_parser.py              # BAM parsing utilities
│   ├── read_classifier.py         # Classification logic (MATLAB and Python)
│   └── matlab_loader.py           # MATLAB .mat file loader
└── outputs/
    ├── bam_comparison/
    │   ├── matlab_alignment.sam
    │   ├── python_alignment.sam
    │   ├── bam_comparison_differences.csv
    │   ├── green_mismatch_reads.fastq
    │   └── comparison_summary.txt
    └── read_count_comparison/
        ├── read_count_comparison_summary.csv
        ├── differentially_labelled_reads.csv
        └── per_junction_read_lists.csv
```

### Key Parameters

**From MATLAB code**:
- `read_length`: From `fastq_read_length.mat` (typically 150)
- `req_ovrlp` (min_overlap_len): From `prms.REQ_OVERLAP` (typically 12)
- `min_bp_in_frame`: 10
- Read length tolerance: 0.9 to 1.1 (hardcoded in `count_reads.m`)

**From Python code**:
- `avg_read_length`: From config or calculated (typically 150)
- `min_overlap_len`: From config (typically 12)
- `read_length_tolerance`: 0.1 (default)
- `min_bp_in_frame`: Not explicitly used in Python, but similar logic in `max_dist_from_junction` (default 10)

**Junction-specific**:
- Junction length: From FASTA file or BAM header
- Junction types: 1-7 (corresponding to different junction orientations)

---

## Expected Outputs Summary

### Part 1 Outputs:
1. `matlab_alignment.sam` - SAM conversion of MATLAB BAM
2. `python_alignment.sam` - SAM conversion of Python BAM
3. `bam_comparison_differences.csv` - Detailed CSV of all differing alignments
4. `green_mismatch_reads.fastq` - FASTQ of reads with different mappings
5. `comparison_summary.txt` - Text summary of differences

### Part 2 Outputs:
1. `read_count_comparison_summary.csv` - Side-by-side count comparison
2. `differentially_labelled_reads.csv` - Reads with different classifications
3. `per_junction_read_lists.csv` - Complete lists of reads per category per junction
4. MATLAB data CSVs (if not already present):
   - `bamreads__*.csv`
   - `biomap__*.csv`

---

## Validation and Quality Checks

1. **BAM File Validation**:
   - Verify both BAM files are valid (use `samtools quickcheck`)
   - Check BAM headers match (same references, same lengths)

2. **Read ID Matching**:
   - Ensure read IDs are extracted correctly from both BAM files
   - Handle potential differences in read ID formatting

3. **Classification Logic Verification**:
   - Manually verify a few example reads to ensure classification logic matches MATLAB/Python code
   - Check edge cases (reads at exact junction boundaries)

4. **Data Completeness**:
   - Verify all junction types (1-7) are present in both pipelines
   - Check that all expected reads are accounted for

5. **Tag Completeness**:
   - Verify that all tags present in BAM files are captured in CSV outputs
   - Check that tags not present for a read are marked as "N/A" (not missing columns)
   - Compare tag presence between MATLAB and Python BAM files to identify tag differences
   - Ensure SAM files contain all tags (verify with `samtools view` output)

---

## Notes and Considerations

1. **Coordinate Systems**:
   - MATLAB BioMap uses 1-based coordinates
   - Python pysam uses 0-based coordinates
   - Ensure consistent conversion when comparing

2. **CIGAR String Interpretation**:
   - MATLAB checks for CIGAR-only-M via string parsing
   - Python checks via `cigartuples` (operation code 0 = M)
   - Verify these are equivalent

3. **Read Filtering**:
   - MATLAB filters by `RdFull` (CIGAR-only-M)
   - Python filters by checking `cigartuples` for non-M operations
   - Both filter by alignment quality (AS, NM tags)
   - Verify filtering criteria match exactly

4. **Secondary/Supplementary Alignments**:
   - MATLAB may handle these differently than Python
   - Check if both pipelines skip secondary/supplementary or count them

5. **Junction Type Mapping**:
   - Verify junction type names match between MATLAB and Python
   - MATLAB uses '1'..'7' as strings, Python uses `JunctionType` enum

6. **Bowtie2 Parameters**:
   - Verify both pipelines use identical bowtie2 parameters
   - Check: `--local` vs `--end-to-end`, `-k`, `--score-min`, `--mp`, `--reorder`

---

## Next Steps After Comparison

Once differences are identified:
1. Analyze root causes (bowtie2 differences, BAM parsing differences, classification logic differences)
2. Determine which pipeline behavior is correct (if applicable)
3. Propose fixes to align Python behavior with MATLAB (or vice versa)
4. Re-run comparison after fixes to verify alignment

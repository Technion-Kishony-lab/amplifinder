# TnJc2 Naming Convention Renaming Plan

## Overview
Standardize naming to use TN terminology consistently, follow verb-based naming pattern, and support singular/plural variable naming.

## Naming Convention

### Pattern
- **Class names**: `[Verb]TnJc2` (PascalCase, e.g., `CoveredTnJc2`)
- **Singular variables**: `[verb]_tnjc2` (snake_case, e.g., `covered_tnjc2`)
- **Plural variables**: `[verb]_tnjc2s` (snake_case with 's', e.g., `covered_tnjc2s`)
- **File names**: `tnjc2_[verb].csv` (snake_case, e.g., `tnjc2_covered.csv`)

## All Step Classes

| Step Name | Renamed? | New Name |
|-----------|----------|----------|
| `InitializingStep` | No | - |
| `GetRefGenomeStep` | No | - |
| `LocateTNsUsingGenbankStep` | No | - |
| `LocateTNsUsingISfinderStep` | No | - |
| `BreseqStep` | No | - |
| `CreateRefTnJcStep` | No | - |
| `CreateRefTnEndSeqsStep` | No | - |
| `CreateTnJcStep` | No | - |
| `CreateTnJc2Step` | **Yes** | `PairTnJc2Step` |
| `CalcAmpliconCoverageStep` | **Yes** | `CalcTnJc2AmpliconCoverageStep` |
| `ClassifyStructureStep` | **Yes** | `ClassifyTnJc2StructureStep` |
| `FilterCandidatesStep` | **Yes** | `FilterTnJc2CandidatesStep` |
| `CreateSyntheticJunctionsStep` | No | - |
| `AlignReadsToJunctionsStep` | No | - |
| `AnalyzeAlignmentsStep` | **Yes** | `AnalyzeTnJc2AlignmentsStep` |
| `ClassifyCandidatesStep` | **Yes** | `ClassifyTnJc2CandidatesStep` |
| `ExportStep` | **Yes** | `ExportTnJc2Step` |

**Note**: `ISfinderStep` is a deprecated alias for `LocateTNsUsingISfinderStep` (not renamed).

## Changes Summary

| Step (Old → New) | Record (Old → New) | Variable (singular) | Variable (plural) | File Name (Old → New) |
|------------------|-------------------|---------------------|-------------------|----------------------|
| `CreateTnJc2Step` → `PairTnJc2Step` | `TnJc2` → `RawTnJc2` | `tnjc2` → `raw_tnjc2` | `tnjc2` → `raw_tnjc2s` | `tn_jc2.csv` → `tnjc2_raw.csv` |
| `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep` | `CoveredTnJc2` → `CoveredTnJc2` | `covered_tnjc2` → `covered_tnjc2` | `covered_tnjc2` → `covered_tnjc2s` | `tn_jc2_covered.csv` → `tnjc2_covered.csv` |
| `ClassifyStructureStep` → `ClassifyTnJc2StructureStep` | `ClassifiedTnJc2` → `ClassifiedTnJc2` | `classified_tnjc2` → `classified_tnjc2` | `classified_tnjc2` → `classified_tnjc2s` | `tn_jc2_classified.csv` → `tnjc2_classified.csv` |
| `FilterCandidatesStep` → `FilterTnJc2CandidatesStep` | `CandidateTnJc2` → `FilteredTnJc2` | `candidate` → `filtered_tnjc2` | `candidates` → `filtered_tnjc2s` | `tn_jc2_candidates.csv` → `tnjc2_filtered.csv` |
| `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep` | `AnalyzedTnJc2` → `AnalyzedTnJc2` | `analyzed_tnjc2` → `analyzed_tnjc2` | `analyzed_tnjc2` → `analyzed_tnjc2s` | `tn_jc2_analyzed.csv` → `tnjc2_analyzed.csv` |
| `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep` | `AnalyzedTnJc2` → `AnalyzedTnJc2` | `analyzed_tnjc2` → `analyzed_tnjc2` | `analyzed_tnjc2` → `analyzed_tnjc2s` | (no file change) |
| `ExportStep` → `ExportTnJc2Step` | `ISJC2Export` → `ExportedTnJc2` | `isjc2` → `exported_tnjc2` | `isjc2` → `exported_tnjc2s` | `ISJC2.csv` → `tnjc2_exported.csv` |

## Detailed Changes

### 1. Step Classes

#### `CreateTnJc2Step` → `PairTnJc2Step`
- **File**: `amplifinder/steps/create_tnjc2.py`
- **Class name**: `CreateTnJc2Step` → `PairTnJc2Step`
- **Reason**: More accurately describes the operation (pairing junctions), includes TnJc2 in name

#### `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
- **File**: `amplifinder/steps/calc_amplicon_coverage.py`
- **Class name**: `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
- **Reason**: Clarifies that it calculates coverage for TnJc2 data

#### `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`
- **File**: `amplifinder/steps/classify_structure.py`
- **Class name**: `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`
- **Reason**: Clarifies that it classifies TnJc2 structure

#### `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
- **File**: `amplifinder/steps/filter_candidates.py`
- **Class name**: `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
- **Reason**: Clarifies that it filters TnJc2 candidates

#### `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`
- **File**: `amplifinder/steps/analyze_alignments.py`
- **Class name**: `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`
- **Reason**: Clarifies that it analyzes alignments for TnJc2 candidates

#### `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`
- **File**: `amplifinder/steps/classify_candidates.py`
- **Class name**: `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`
- **Reason**: Clarifies that it classifies TnJc2 candidates

#### `ExportStep` → `ExportTnJc2Step`
- **File**: `amplifinder/steps/export.py`
- **Class name**: `ExportStep` → `ExportTnJc2Step`
- **Reason**: Clarifies that it exports TnJc2 data

### 2. Record Types

#### `TnJc2` → `RawTnJc2`
- **File**: `amplifinder/data_types/record_types.py`
- **Class name**: `TnJc2` → `RawTnJc2`
- **Reason**: Follows verb pattern, indicates raw/unprocessed pairs
- **Inheritance chain**: All other types inherit from this

#### `CandidateTnJc2` → `FilteredTnJc2`
- **File**: `amplifinder/data_types/record_types.py`
- **Class name**: `CandidateTnJc2` → `FilteredTnJc2`
- **Reason**: Follows verb pattern (from `FilterCandidatesStep`)
- **Inheritance**: `ClassifiedTnJc2` → `FilteredTnJc2` → `AnalyzedTnJc2`

#### `ISJC2Export` → `ExportedTnJc2`
- **File**: `amplifinder/data_types/record_types.py`
- **Class name**: `ISJC2Export` → `ExportedTnJc2`
- **Reason**: Follows verb pattern, uses TN terminology

### 3. File Names

#### Base file pattern
- `tn_jc2.csv` → `tnjc2_raw.csv`
- `tn_jc2_covered.csv` → `tnjc2_covered.csv`
- `tn_jc2_classified.csv` → `tnjc2_classified.csv`
- `tn_jc2_candidates.csv` → `tnjc2_filtered.csv`
- `tn_jc2_analyzed.csv` → `tnjc2_analyzed.csv`

#### Export file
- `ISJC2.csv` → `tnjc2_exported.csv`
- **Note**: Legacy compatibility may require keeping `ISJC2.csv` as an alias

### 4. Variable Names

#### Singular variables
- `tnjc2` → `raw_tnjc2`
- `covered_tnjc2` → `covered_tnjc2` (no change)
- `classified_tnjc2` → `classified_tnjc2` (no change)
- `candidate` → `filtered_tnjc2`
- `candidates` → `filtered_tnjc2` (singular) or `filtered_tnjc2s` (plural)
- `analyzed_tnjc2` → `analyzed_tnjc2` (no change)
- `isjc2` (export) → `exported_tnjc2`

#### Plural variables
- Use `s` suffix: `[verb]_tnjc2s`
- Examples: `raw_tnjc2s`, `covered_tnjc2s`, `filtered_tnjc2s`, etc.

## Files to Update

### Core Data Types
- `amplifinder/data_types/record_types.py`
  - Rename `TnJc2` → `RawTnJc2`
  - Rename `CandidateTnJc2` → `FilteredTnJc2`
  - Rename `ISJC2Export` → `ExportedTnJc2`
  - Update all inheritance chains

- `amplifinder/data_types/__init__.py`
  - Update exports for renamed classes

### Steps
- `amplifinder/steps/create_tnjc2.py`
  - Rename `CreateTnJc2Step` → `PairTnJc2Step`
  - Update type hints: `TnJc2` → `RawTnJc2`
  - Update variable names

- `amplifinder/steps/calc_amplicon_coverage.py`
  - Rename `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
  - Update type hints: `TnJc2` → `RawTnJc2`
  - Update variable names

- `amplifinder/steps/classify_structure.py`
  - Rename `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`
  - Update variable names to plural where appropriate

- `amplifinder/steps/filter_candidates.py`
  - Rename `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
  - Update type hints: `CandidateTnJc2` → `FilteredTnJc2`
  - Update variable names: `candidate` → `filtered_tnjc2`, `candidates` → `filtered_tnjc2s`

- `amplifinder/steps/analyze_alignments.py`
  - Rename `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`
  - Update variable names

- `amplifinder/steps/classify_candidates.py`
  - Rename `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`
  - Update variable names

- `amplifinder/steps/export.py`
  - Rename `ExportStep` → `ExportTnJc2Step`
  - Update type hints: `ISJC2Export` → `ExportedTnJc2`
  - Update file names: `ISJC2.csv` → `tnjc2_exported.csv`
  - Update variable names

- `amplifinder/steps/io_naming.py`
  - Update `DEFAULT_FILENAMES` mapping for all renamed types
  - Update file name patterns

### Pipeline
- `amplifinder/pipeline.py`
  - Update all type hints
  - Update all step class names:
    - `CreateTnJc2Step` → `PairTnJc2Step`
    - `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
    - `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`
    - `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
    - `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`
    - `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`
    - `ExportStep` → `ExportTnJc2Step`
  - Update variable names throughout

### Tests
- `tests/test_steps/test_create_tnjc2.py`
  - Update step class name: `CreateTnJc2Step` → `PairTnJc2Step`
  - Update type hints: `TnJc2` → `RawTnJc2`
  - Update variable names

- `tests/test_steps/test_calc_coverage.py`
  - Update step class name: `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
  - Update type hints: `TnJc2` → `RawTnJc2`

- `tests/test_steps/test_classify_structure.py`
  - Update step class name: `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`

- `tests/test_steps/test_filter_candidates.py`
  - Update step class name: `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
  - Update type hints: `CandidateTnJc2` → `FilteredTnJc2`

- `tests/test_steps/test_analyze_alignments.py`
  - Update step class name: `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`

- `tests/test_steps/test_classify_candidates.py`
  - Update step class name: `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`

- `tests/test_steps/test_export.py`
  - Update step class name: `ExportStep` → `ExportTnJc2Step`
  - Update type hints: `ISJC2Export` → `ExportedTnJc2`
  - Update file name assertions

- `tests/test_integration/test_pipeline.py`
  - Update all step class names
  - Update all type hints and variable names
  - Update file name checks

### Other Files
- `amplifinder/steps/__init__.py`
  - Update all step class exports:
    - `CreateTnJc2Step` → `PairTnJc2Step`
    - `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
    - `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`
    - `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
    - `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`
    - `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`
    - `ExportStep` → `ExportTnJc2Step`

- `amplifinder/visualization/visualize.py`
  - Update type hints

## Migration Strategy

### Phase 1: Type Definitions
1. Rename classes in `record_types.py`
2. Update `__init__.py` exports
3. Update `io_naming.py` file mappings

### Phase 2: Step Classes
1. Rename all TnJc2-related step classes:
   - `CreateTnJc2Step` → `PairTnJc2Step`
   - `CalcAmpliconCoverageStep` → `CalcTnJc2AmpliconCoverageStep`
   - `ClassifyStructureStep` → `ClassifyTnJc2StructureStep`
   - `FilterCandidatesStep` → `FilterTnJc2CandidatesStep`
   - `AnalyzeAlignmentsStep` → `AnalyzeTnJc2AlignmentsStep`
   - `ClassifyCandidatesStep` → `ClassifyTnJc2CandidatesStep`
   - `ExportStep` → `ExportTnJc2Step`
2. Update all step type hints
3. Update step variable names

### Phase 3: Pipeline Integration
1. Update `pipeline.py` with new names
2. Update all variable names to follow convention

### Phase 4: Tests
1. Update all test files
2. Run test suite to verify changes

### Phase 5: File Names
1. Update file name generation in `io_naming.py`
2. Consider backward compatibility for export files

## Notes
- All changes follow the `[Verb]TnJc2` pattern for consistency
- Variable names use snake_case with `s` suffix for plurals
- File names use snake_case pattern: `tnjc2_[verb].csv`
- Step names reflect their operations (verb-based)


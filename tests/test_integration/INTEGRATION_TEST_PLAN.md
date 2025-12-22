# Full Integration Test Plan

## Overview
Run complete Python pipeline and compare outputs with MATLAB reference outputs.

## Test Data
- **Test isolate**: SRR25242877
- **Ancestor**: SRR25242906  
- **Reference**: U00096 (E. coli K-12)
- **MATLAB outputs**: `/zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877/`
- **Python outputs**: `/zdata/user-data/rkishony/AmpliFinder_test/python_outputs/output/U00096/SRR25242877/SRR25242877/`
- **Breseq outputs**: `/zdata/user-data/idan/small_projects/breseq_to_amplification_from_dropbox/Breseq1/SRR25242877/`
- **FASTQ files**: `/zdata/user-data/idan/small_projects/breseq_to_amplification_from_dropbox/morbidostat_sra/SRR25242877/`

## MATLAB Output Files to Compare
1. **ISJC2.xlsx** - All analyzed candidates (155 rows)
   - Columns: `iso`, `Reference`, `Positions_in_chromosome_1/2`, `Direction_in_chromosome_1/2`, `amplicon_length`, `IS_element`, `IS_direction`, `median_copy_number`, `mode_copy_number`, `Ancestor`
   
2. **classified_amplifications.xlsx** - Filtered candidates
   
3. **candidate_amplifications.xlsx** - Additional filtered candidates

## Implementation Steps

### 1. Create Full Pipeline Test Function
**File**: `tests/test_integration/test_pipeline.py`

#### Output Directory Strategy

**Outputs placed next to MATLAB outputs for side-by-side comparison:**
- MATLAB outputs: `/zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877/`
- Python outputs: `/zdata/user-data/rkishony/AmpliFinder_test/python_outputs/output/U00096/SRR25242877/SRR25242877/`

This allows easy comparison of outputs and inspection of results.

#### Full Test Implementation
```python
@pytest.mark.slow
def test_full_pipeline_matches_matlab(isolate_srr25242877):
    """Run full pipeline and compare with MATLAB outputs."""
    from amplifinder.config import Config, get_run_dir
    from amplifinder.pipeline import Pipeline
    
    # Setup output directory next to MATLAB outputs
    matlab_output_dir = isolate_srr25242877["matlab_output"]
    # matlab_output_dir = /zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877
    # Go up to AmpliFinder_test, then create python_outputs directory
    test_output_root = matlab_output_dir.parent.parent.parent / "python_outputs"
    # test_output_root = /zdata/user-data/rkishony/AmpliFinder_test/python_outputs
    
    # Setup config using pre-computed breseq output
    config = Config(
        iso_path=isolate_srr25242877["fastq_path"],
        ref_name="U00096",
        anc_path=None,  # Can add ancestor later
        iso_name="SRR25242877",
        output_dir=test_output_root / "output",
        ref_path=test_output_root / "genomesDB",
        iso_breseq_path=isolate_srr25242877["breseq_path"],
        ncbi=True,
        use_isfinder=False,  # Use GenBank for consistency
    )
    
    # Run full pipeline
    pipeline = Pipeline(config)
    result = pipeline.run()
    
    # Verify outputs exist
    run_dir = get_run_dir(config)
    # run_dir = test_output_root / "output" / "U00096" / "SRR25242877" / "SRR25242877"
    assert (run_dir / "ISJC2.csv").exists()
    assert (run_dir / "candidate_amplifications.csv").exists()
    
    # Load Python outputs
    python_isjc2 = pd.read_csv(run_dir / "ISJC2.csv")
    
    # Load MATLAB outputs
    matlab_isjc2 = pd.read_excel(
        isolate_srr25242877["matlab_output"] / "ISJC2.xlsx"
    )
    
    # Compare outputs (1-to-1 matching required)
    compare_isjc2_outputs(python_isjc2, matlab_isjc2)
```

### 2. Implement Comparison Functions
**File**: `tests/test_integration/matlab_compare.py`

#### Column Mapping
- MATLAB `Positions_in_chromosome_1/2` → Python `Positions_in_chromosome` (split by `-`)
- MATLAB `Direction_in_chromosome_1/2` → Python `Direction_in_chromosome` (split by `/`)
- MATLAB `IS_direction` → Not in Python export (may need to derive from `tn_orientations`)

#### Comparison Strategy
1. **1-to-1 matching**: Every MATLAB junction must have exactly one Python match (and vice versa)
2. **Position matching**: Match junctions by chromosome positions (±tolerance)
3. **Amplicon length**: Compare lengths (±tolerance) for matched pairs
4. **IS elements**: Compare IS_element lists (order-independent) for matched pairs
5. **Copy numbers**: Compare median/mode copy numbers (±tolerance for NaN handling) for matched pairs

#### Matching Algorithm
```python
def match_junctions(python_df, matlab_df, pos_tolerance=5):
    """Match Python and MATLAB junctions by position (1-to-1 matching).
    
    Returns:
        matches: List of (python_idx, matlab_idx) tuples
        python_matched: Set of matched Python indices
        matlab_matched: Set of matched MATLAB indices
    """
    matches = []
    python_matched = set()
    matlab_matched = set()
    
    # Build candidate matches (may have multiple candidates per junction)
    candidates = []
    for i, matlab_row in matlab_df.iterrows():
        matlab_pos1 = matlab_row['Positions_in_chromosome_1']
        matlab_pos2 = matlab_row['Positions_in_chromosome_2']
        
        for j, python_row in python_df.iterrows():
            # Parse Python positions
            py_positions = python_row['Positions_in_chromosome'].split('-')
            py_pos1, py_pos2 = int(py_positions[0]), int(py_positions[1])
            
            # Check if positions match within tolerance
            if (abs(py_pos1 - matlab_pos1) <= pos_tolerance and
                abs(py_pos2 - matlab_pos2) <= pos_tolerance):
                candidates.append((j, i, abs(py_pos1 - matlab_pos1) + abs(py_pos2 - matlab_pos2)))
    
    # Sort by distance (best matches first)
    candidates.sort(key=lambda x: x[2])
    
    # Greedy matching: assign best matches first, ensuring 1-to-1
    for py_idx, matlab_idx, distance in candidates:
        if py_idx not in python_matched and matlab_idx not in matlab_matched:
            matches.append((py_idx, matlab_idx))
            python_matched.add(py_idx)
            matlab_matched.add(matlab_idx)
    
    return matches, python_matched, matlab_matched
```

### 3. Comparison Metrics
- **1-to-1 match requirement**: Every junction must have exactly one match
  - All MATLAB junctions must match exactly one Python junction
  - All Python junctions must match exactly one MATLAB junction
  - No duplicates, no missing junctions
- **Match count**: Should be `len(matches) == len(matlab_df) == len(python_df)`
- **Position accuracy**: Position differences for matched junctions (should be ≤tolerance)
- **Amplicon length accuracy**: Length differences for matched junctions
- **IS element agreement**: IS_element sets should match for matched junctions
- **Copy number correlation**: Correlation coefficient for matched junctions (handling NaN)

### 4. Test Assertions
```python
def compare_isjc2_outputs(python_df, matlab_df, 
                         pos_tolerance=5,
                         length_tolerance=10):
    """Compare ISJC2 outputs with 1-to-1 matching requirement."""
    matches, py_matched, matlab_matched = match_junctions(
        python_df, matlab_df, pos_tolerance
    )
    
    # Require 1-to-1 match: every junction must have exactly one match
    assert len(matches) == len(matlab_df), (
        f"Not all MATLAB junctions matched: {len(matches)}/{len(matlab_df)} matched. "
        f"Missing MATLAB junctions: {len(matlab_df) - len(matlab_matched)}"
    )
    
    assert len(matches) == len(python_df), (
        f"Not all Python junctions matched: {len(matches)}/{len(python_df)} matched. "
        f"Extra Python junctions: {len(python_df) - len(py_matched)}"
    )
    
    # Verify no duplicates
    assert len(py_matched) == len(matches), "Duplicate Python junction matches"
    assert len(matlab_matched) == len(matches), "Duplicate MATLAB junction matches"
    
    # Compare matched junctions
    for py_idx, matlab_idx in matches:
        py_row = python_df.iloc[py_idx]
        matlab_row = matlab_df.iloc[matlab_idx]
        
        # Compare amplicon length
        py_len = py_row['amplicon_length']
        matlab_len = matlab_row['amplicon_length']
        assert abs(py_len - matlab_len) <= length_tolerance, (
            f"Length mismatch: Python={py_len}, MATLAB={matlab_len}"
        )
        
        # Compare IS elements (order-independent)
        py_is = set(str(py_row['IS_element']).split(','))
        matlab_is = set(str(matlab_row['IS_element']).split(','))
        assert py_is == matlab_is, (
            f"IS elements mismatch: Python={py_is}, MATLAB={matlab_is}"
        )
```

### 5. Handle Differences
- **NaN copy numbers**: MATLAB has NaN for some candidates (coverage calculation differences)
- **Column name differences**: Map between formats
- **IS_direction**: May need to derive from `tn_orientations` field
- **Position format**: Parse `Positions_in_chromosome` string vs separate columns

### 6. Test Configuration
- Use `@pytest.mark.slow` marker
- Skip if MATLAB outputs not available
- Use pre-computed breseq output (don't run breseq in test)
- Output directory: `/zdata/user-data/rkishony/AmpliFinder_test/python_outputs/`
  - Allows side-by-side comparison with MATLAB outputs
  - Requires write access to test data directory

### 7. Ancestor Support (Optional)
If ancestor data available:
```python
config = Config(
    iso_path=isolate_srr25242877["fastq_path"],
    anc_path=isolate_srr25242906["fastq_path"],
    anc_name="SRR25242906",
    anc_breseq_path=isolate_srr25242906["breseq_path"],
    # ... rest of config
)
```

## Expected Outcomes
- Full pipeline runs successfully
- **1-to-1 junction matching**: Every junction must have exactly one match
  - All MATLAB junctions matched to Python junctions
  - All Python junctions matched to MATLAB junctions
  - No duplicates, no missing junctions
- **Matched junction properties**:
  - Position differences ≤5bp
  - Amplicon length differences ≤10bp
  - IS element sets match (order-independent)
  - Copy numbers correlate (handling NaN differences)

## Files to Modify
1. `tests/test_integration/test_pipeline.py` - Add full pipeline test
2. `tests/test_integration/matlab_compare.py` - Implement comparison functions
3. Update `test_pipeline_matches_matlab_tn_count` to actually run pipeline

## Testing Order
1. Implement comparison functions first (test with existing partial outputs)
2. Add full pipeline test
3. Run and debug differences
4. Adjust tolerances/assertions based on results


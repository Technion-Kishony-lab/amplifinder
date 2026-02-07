#!/usr/bin/env python3
"""Create a test coverage plot with schematic junctions."""
import sys
sys.path.insert(0, 'tests')
import pytest
from pathlib import Path
from conftest import *  # Import all fixtures

def test_create_coverage_plot_with_schematic(analyzed_tnjc2_record, tiny_genome):
    """Create a coverage plot we can actually view."""
    import numpy as np
    from amplifinder.data_types.jc_types import JunctionType, JcCall
    
    # Import directly to avoid circular import
    import sys
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "plot_amp_coverage",
        "/zdata/user-data/rkishony/amplifinder/amplifinder/visualization/plot_amp_coverage.py"
    )
    plot_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(plot_module)
    plot_amplicon_coverage = plot_module.plot_amplicon_coverage
    
    # Create dummy coverage
    iso_coverage = {}
    anc_coverage = {}
    for scaffold_name, scaffold in tiny_genome.scaffolds.items():
        scaffold_length = len(scaffold)
        iso_coverage[scaffold_name] = np.random.poisson(10, scaffold_length).astype(float)
        anc_coverage[scaffold_name] = np.random.poisson(10, scaffold_length).astype(float)
    
    # Make amplicon region higher coverage (so we can see it)
    amp_start = analyzed_tnjc2_record.left
    amp_end = analyzed_tnjc2_record.right
    iso_coverage[analyzed_tnjc2_record.scaf][amp_start:amp_end] *= 3
    
    # Create junction calls - ALL POSITIVE so we see schematic elements
    jc_calls = {jt: JcCall.POS for jt in JunctionType}
    jc_calls_anc = {jt: JcCall.POS for jt in JunctionType}
    
    output_path = Path("/tmp/test_coverage_with_schematic.png")
    
    plot_amplicon_coverage(
        tnjc2=analyzed_tnjc2_record,
        jc_calls=jc_calls,
        jc_calls_anc=jc_calls_anc,
        iso_scafs_to_covs=iso_coverage,
        anc_scafs_to_covs=anc_coverage,
        output_path=output_path,
    )
    
    assert output_path.exists(), f"Plot not created at {output_path}"
    size = output_path.stat().st_size
    print(f"\n✓ Plot created successfully!")
    print(f"  Path: {output_path}")
    print(f"  Size: {size:,} bytes")
    print(f"\nView with: display {output_path}")
    
if __name__ == "__main__":
    # Run via pytest to use fixtures
    sys.exit(pytest.main([__file__, "-v", "-s"]))

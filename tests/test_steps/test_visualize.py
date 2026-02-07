"""Tests for visualization module."""

from amplifinder.visualization.plot_amp_coverage import plot_amplicon_coverage
from amplifinder.data_types import JunctionType, JcCall
import numpy as np


def test_plot_amplicon_coverage(analyzed_tnjc2_record, tmp_path, tiny_genome):
    """Should create coverage plot."""
    # Create dummy coverage dict (scaffold name -> coverage array)
    iso_coverage = {}
    anc_coverage = {}
    for scaffold_name, scaffold in tiny_genome.scaffolds.items():
        scaffold_length = len(scaffold)
        iso_coverage[scaffold_name] = np.random.poisson(10, scaffold_length).astype(float)
        anc_coverage[scaffold_name] = np.random.poisson(10, scaffold_length).astype(float)

    # Create junction calls (all positive for visualization)
    jc_calls = {jt: JcCall.POS for jt in JunctionType}
    jc_calls_anc = {jt: JcCall.POS for jt in JunctionType}

    output_path = tmp_path / "coverage_plot.png"

    plot_amplicon_coverage(
        tnjc2=analyzed_tnjc2_record,
        jc_calls=jc_calls,
        jc_calls_anc=jc_calls_anc,
        iso_scafs_to_covs=iso_coverage,
        anc_scafs_to_covs=anc_coverage,
        output_path=output_path,
    )

    assert output_path.exists()

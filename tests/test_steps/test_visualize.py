"""Tests for visualization module."""

from amplifinder.visualization import plot_amplicon_coverage
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

    output_path = tmp_path / "coverage_plot.png"

    plot_amplicon_coverage(
        iso_scafs_to_covs=iso_coverage,
        anc_scafs_to_covs=anc_coverage,
        tnjc2=analyzed_tnjc2_record,
        output_path=output_path,
    )

    assert output_path.exists()

"""Tests for visualization module."""

from amplifinder.visualization import plot_candidate_coverage
import numpy as np


def test_plot_candidate_coverage(analyzed_tnjc2_record, tmp_path, tiny_genome):
    """Should create coverage plot."""
    # Create dummy coverage arrays
    genome_length = len(tiny_genome)
    iso_coverage = np.random.poisson(10, genome_length).astype(float)
    anc_coverage = np.random.poisson(10, genome_length).astype(float)

    output_path = tmp_path / "coverage_plot.png"

    plot_candidate_coverage(
        candidate=analyzed_tnjc2_record,
        iso_coverage=iso_coverage,
        anc_coverage=anc_coverage,
        genome=tiny_genome,
        output_path=output_path,
        show=False,
    )

    assert output_path.exists()

"""Tests for visualization module."""

import pytest
from pathlib import Path
from amplifinder.visualization import plot_candidate_coverage
from amplifinder.data_types import AnalyzedTnJc2, RawEvent, Orientation, Coverage
import numpy as np


@pytest.fixture
def sample_candidate(tmp_path):
    """Create a sample AnalyzedTnJc2."""
    return AnalyzedTnJc2(
        jc_num_L=1, jc_num_R=2,
        scaf="chr1",
        pos_scaf_L=100, pos_scaf_R=200,
        pos_tn_L=10, pos_tn_R=20,
        dir_scaf_L=Orientation.FORWARD, dir_scaf_R=Orientation.REVERSE,
        dir_tn_L=Orientation.FORWARD, dir_tn_R=Orientation.REVERSE,
        tn_ids=[1], tn_orientations=[Orientation.FORWARD],
        span_origin=False,
        amplicon_length=100, complementary_length=900,
        ref_name="U00096", iso_name="sample1",
        amplicon_coverage=2.0, scaf_coverage=Coverage(mean=1.0, median=1.0, mode=1.0),
        copy_number=2.0, amplicon_coverage_mode=2.0,
        raw_event=RawEvent.FLANKED,
        shared_tn_ids=[1], chosen_tn_id=1,
        analysis_dir="jc_100_200_001_L150",
        jc_cov_left=[0]*7, jc_cov_right=[0]*7, jc_cov_spanning=[0]*7,
        isolate_architecture=RawEvent.FLANKED,
        ancestor_architecture=None,
        event="flanked",
        event_modifiers=[],
    )


def test_plot_candidate_coverage(sample_candidate, tmp_path):
    """Should create coverage plot."""
    # Create dummy coverage arrays
    genome_length = 1000
    iso_coverage = np.random.poisson(10, genome_length).astype(float)
    anc_coverage = np.random.poisson(10, genome_length).astype(float)
    
    output_path = tmp_path / "coverage_plot.png"
    
    plot_candidate_coverage(
        candidate=sample_candidate,
        iso_coverage=iso_coverage,
        anc_coverage=anc_coverage,
        genome_length=genome_length,
        output_path=output_path,
        show=False,
    )
    
    assert output_path.exists()

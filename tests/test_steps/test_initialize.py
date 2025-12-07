"""Tests for InitializingStep."""

import pytest
from amplifinder.steps import InitializingStep


@pytest.fixture
def init_step(tmp_path):
    return InitializingStep(output_dir=tmp_path / "output", iso_name="sample1")


def test_creates_directories(init_step, tmp_path):
    """Should create output and isolate directories."""
    result = init_step.run_and_read_outputs()
    
    assert (tmp_path / "output").exists()
    assert (tmp_path / "output" / "sample1").exists()
    assert result == tmp_path / "output" / "sample1"


def test_skips_if_exists(init_step):
    """Should skip if output already exists."""
    assert init_step.run() is True
    assert init_step.run() is False


def test_force_reruns(tmp_path):
    """Force=True should re-run even if output exists."""
    step = InitializingStep(output_dir=tmp_path / "output", iso_name="sample1")
    step.run()
    
    step_force = InitializingStep(output_dir=tmp_path / "output", iso_name="sample1", force=True)
    assert step_force.run() is True


"""Tests for InitializingStep."""

import pytest
from amplifinder.steps import InitializingStep


@pytest.fixture
def init_step(tmp_path):
    return InitializingStep(output_dir=tmp_path / "output", anc_name="ancestor1", iso_name="sample1")


def test_creates_directories(init_step, tmp_path):
    """Should create output, ancestor, and isolate directories."""
    result = init_step.run()

    assert (tmp_path / "output").exists()
    assert (tmp_path / "output" / "ancestor1").exists()
    assert (tmp_path / "output" / "ancestor1" / "sample1").exists()
    assert result == tmp_path / "output" / "ancestor1" / "sample1"


def test_self_ancestor(tmp_path):
    """When iso_name == anc_name, creates {name}/{name}/ structure."""
    step = InitializingStep(output_dir=tmp_path / "output", anc_name="sample1", iso_name="sample1")
    result = step.run()

    assert (tmp_path / "output" / "sample1" / "sample1").exists()
    assert result == tmp_path / "output" / "sample1" / "sample1"


def test_skips_if_exists(init_step):
    """Should skip if output already exists."""
    init_step.run()
    assert init_step.run_count == 1
    init_step.run()
    assert init_step.run_count == 1  # didn't run again


def test_force_reruns(tmp_path):
    """Force=True should re-run even if output exists."""
    step = InitializingStep(output_dir=tmp_path / "output", anc_name="ancestor1", iso_name="sample1")
    step.run()

    step_force = InitializingStep(output_dir=tmp_path / "output", anc_name="ancestor1", iso_name="sample1", force=True)
    step_force.run()
    assert step_force.run_count == 1

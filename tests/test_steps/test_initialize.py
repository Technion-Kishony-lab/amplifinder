"""Tests for InitializingStep."""

import pytest
from amplifinder.steps import InitializingStep
from amplifinder.config import Config


@pytest.fixture
def sample_config(tmp_path):
    """Create a sample config for testing."""
    return Config(
        iso_fastq_path=tmp_path / "isolate.fastq",
        anc_fastq_path=tmp_path / "ancestor.fastq",  # Need fastq_path for has_ancestor
        ref_name="U00096",
        iso_name="sample1",
        anc_name="ancestor1",
        output_dir=tmp_path / "output",
    )


@pytest.fixture
def init_step_factory(sample_config):
    """Factory to create InitializingStep instances with the same config."""
    def make():
        return InitializingStep(config=sample_config)
    return make


def test_creates_directories(init_step_factory, tmp_path):
    """Should create output directories with ref_name/anc_name/iso_name structure."""
    init_step = init_step_factory()
    iso_run_dir, anc_run_dir = init_step.run()

    assert (tmp_path / "output").exists()
    assert (tmp_path / "output" / "U00096").exists()
    assert (tmp_path / "output" / "U00096" / "ancestor1").exists()
    assert (tmp_path / "output" / "U00096" / "ancestor1" / "sample1").exists()
    assert iso_run_dir == tmp_path / "output" / "U00096" / "ancestor1" / "sample1"


def test_creates_config_file(init_step_factory, tmp_path):
    """Should create run_config.yaml file."""
    init_step = init_step_factory()
    result = init_step.run()
    iso_run_dir, anc_run_dir = result

    config_path = iso_run_dir / "run_config.yaml"
    assert config_path.exists()


def test_self_ancestor(tmp_path):
    """When iso_name == anc_name, creates {ref_name}/{name}/{name}/ structure."""
    config = Config(
        iso_fastq_path=tmp_path / "isolate.fastq",
        ref_name="U00096",
        iso_name="sample1",
        anc_name="sample1",  # same as iso_name
        output_dir=tmp_path / "output",
    )
    step = InitializingStep(config=config)
    iso_run_dir, anc_run_dir = step.run()

    assert (tmp_path / "output" / "U00096" / "sample1" / "sample1").exists()
    assert iso_run_dir == tmp_path / "output" / "U00096" / "sample1" / "sample1"


def test_skips_if_exists(init_step_factory):
    """Should skip if output already exists."""
    step1 = init_step_factory()
    step1.run()
    assert step1._artifacts_generated is True

    # New instance should skip due to cached artifacts
    step2 = init_step_factory()
    step2.run()
    assert step2._artifacts_generated is False


def test_force_reruns(tmp_path):
    """Force=True should re-run even if output exists."""
    config = Config(
        iso_fastq_path=tmp_path / "isolate.fastq",
        ref_name="U00096",
        iso_name="sample1",
        anc_name="ancestor1",
        output_dir=tmp_path / "output",
    )
    step = InitializingStep(config=config)
    step.run()

    step_force = InitializingStep(config=config, force=True)
    step_force.run()
    assert step_force._artifacts_generated is True

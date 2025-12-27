"""Tests for config save/load utilities."""

import pytest
from pathlib import Path

from amplifinder.config import Config
from amplifinder.utils.file_utils import ensure_dir


@pytest.fixture
def sample_config(tmp_path):
    """Create a sample config for testing."""
    return Config(
        iso_path=tmp_path / "isolate.fastq",
        ref_name="U00096",
        iso_name="sample1",
        anc_name="ancestor1",
        output_dir=tmp_path / "output",
    )


def test_get_run_dir(sample_config):
    """Should return correct run directory path."""
    run_dir = sample_config.iso_run_dir
    expected = sample_config.output_dir / "U00096" / "ancestor1" / "sample1"
    assert run_dir == expected


def test_save_config(sample_config, tmp_path):
    """Should save config to run_config.yaml."""
    run_dir = tmp_path / "output" / "U00096" / "ancestor1" / "sample1"
    sample_config.save(run_dir)
    
    config_path = run_dir / "run_config.yaml"
    assert config_path.exists()
    
    # Check that it's valid YAML
    import yaml
    with open(config_path) as f:
        config_dict = yaml.safe_load(f)
    
    assert config_dict["ref_name"] == "U00096"
    assert config_dict["iso_name"] == "sample1"
    assert config_dict["anc_name"] == "ancestor1"


def test_load_config_from_run(sample_config, tmp_path):
    """Should load config from run_config.yaml."""
    run_dir = tmp_path / "output" / "U00096" / "ancestor1" / "sample1"
    sample_config.save(run_dir)
    
    loaded_config = Config.load_from_run(run_dir)
    
    assert loaded_config.ref_name == sample_config.ref_name
    assert loaded_config.iso_name == sample_config.iso_name
    assert loaded_config.anc_name == sample_config.anc_name
    assert loaded_config.output_dir == sample_config.output_dir


def test_load_config_missing_file(tmp_path):
    """Should raise FileNotFoundError if config file doesn't exist."""
    run_dir = ensure_dir(tmp_path / "output" / "U00096" / "ancestor1" / "sample1")
    
    with pytest.raises(FileNotFoundError):
        Config.load_from_run(run_dir)

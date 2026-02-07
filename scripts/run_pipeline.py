#!/usr/bin/env python3
"""Simple pipeline runner for debugging - edit CONFIG_FILE below and run."""

from pathlib import Path

from amplifinder.config import Config, load_config
from amplifinder.pipeline import Pipeline

# Edit this to point to your config file
CONFIG_FILE = Path(
    "/zdata/user-data/rkishony/amplifinder/examples/run_all/output/"
    "U00096/SRR5182940/SRR5182911/run_config.yaml"
)

# Load config
config_dict = load_config(CONFIG_FILE)
config = Config(**config_dict)

pipeline = Pipeline(config, verbose=True)
pipeline.run()

#!/bin/bash
# Example: Creating and using config files

# Step 1: Create a complete config file from CLI args
# This generates a YAML with ALL parameters (including defaults)
amplifinder \
    -i data/isolate.fastq \
    -r U00096 \
    -a data/ancestor.fastq \
    --create-config my_experiment.yaml

echo "Config created: my_experiment.yaml"
echo ""
echo "You can now:"
echo "  1. Edit my_experiment.yaml to customize parameters"
echo "  2. Run with: amplifinder --config my_experiment.yaml"
echo ""

# Step 2: Run with the config
# amplifinder --config my_experiment.yaml

# Step 3: Rerun using the saved run_config.yaml
# amplifinder --config output/U00096/ancestor/isolate/run_config.yaml

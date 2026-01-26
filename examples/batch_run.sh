#!/bin/bash
# Example batch processing script
# Processes multiple isolates against a shared ancestor

# Shared parameters
PARAMS="params_example.yaml"
ANCESTOR="data/ancestor.fastq"
REF_NAME="U00096"

# Process all isolates
for ISOLATE in data/isolates/*.fastq; do
    ISO_NAME=$(basename "$ISOLATE" .fastq)
    echo "Processing $ISO_NAME..."
    
    amplifinder \
        -i "$ISOLATE" \
        -a "$ANCESTOR" \
        -r "$REF_NAME" \
        --iso-name "$ISO_NAME" \
        --config "$PARAMS"
done

echo "Batch processing complete!"

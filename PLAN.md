# AmpliFinder Implementation Plan

## Current Package Structure
```
amplifinder/
├── steps/       # Pipeline Step classes (base.py, init.py, get_reference.py, isfinder.py, run_breseq.py)
├── tools/       # External tool runners + parsers (blast.py, breseq.py)
├── utils/       # General utilities (fasta.py)
├── data_types/  # Data structures (genome.py, junction.py, schema.py)
├── data/        # Bundled data (ISfinderDB, field definitions)
├── cli.py       # CLI entry point
├── config.py    # Configuration
└── logger.py    # Logging
```

## Implementation Phases

1. **Scaffold + Config** ✅: Project structure, `pyproject.toml`, config loading, logger, CLI skeleton
2. **Reference + breseq** ✅: Reference genome handling, breseq runner/parser → `JC` table output
3. **IS detection**: Find IS elements in reference, assign IS to junctions → `ISJC` table
4. **Junction pairs**: Combine ISJC pairs, coverage calculation → `ISJC2` table
5. **Classification + Export**: Event classification, FASTA/alignment/BAM analysis, Excel export

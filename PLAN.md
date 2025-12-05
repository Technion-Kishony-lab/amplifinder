# AmpliFinder Implementation Plan

1. **Scaffold + Config**: Project structure, `pyproject.toml`, config loading, logger, CLI skeleton
2. **Reference + breseq**: Reference genome handling, breseq runner/parser → `JC` table output
3. **IS detection**: Find IS elements in reference, assign IS to junctions → `ISJC` table
4. **Junction pairs**: Combine ISJC pairs, coverage calculation → `ISJC2` table
5. **Classification + Export**: Event classification, FASTA/alignment/BAM analysis, Excel export

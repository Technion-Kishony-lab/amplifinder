# AmpliFinder

A MATLAB tool for detecting gene amplifications and deletions via Insertion Sequence (IS) junction analysis from whole-genome sequencing data.

## Overview

AmpliFinder analyzes paired-end sequencing data to identify IS-mediated genomic rearrangements by:
- Aligning reads using breseq
- Detecting junction reads spanning IS elements
- Classifying candidate amplicons/deletions based on coverage

## Usage

```matlab
AmpliFinder('-iso_path', 'path/to/isolate.fastq', ...
            '-anc_path', 'path/to/ancestor.fastq', ...
            '-ref_name', 'CP000000')
```

### Key Parameters

| Parameter | Description |
|-----------|-------------|
| `-iso_path` | Path to isolate FASTQ files |
| `-anc_path` | Path to ancestor FASTQ files |
| `-ref_name` | Reference genome name (NCBI accession) |
| `-iso_name` | Isolate name (default: "isolate") |
| `-NCBI` | Use NCBI reference (default: true) |
| `-ISfinder` | Identify IS elements via ISfinder DB (default: false) |

## Requirements

- MATLAB
- breseq
- bowtie2
- samtools
- BLAST+

## Output

- `ISJC2.xlsx` - Candidate amplicons/deletions table
- `filteredISJC2.xlsx` - Filtered high-confidence calls


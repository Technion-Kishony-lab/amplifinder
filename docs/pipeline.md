# Pipeline

AmpliFinder processes whole-genome FASTQ reads (typically paired-end Illumina) through 16 steps to detect IS-mediated gene amplifications and deletions. **Inputs are directories:** each of `iso_fastq_path` and optional `anc_fastq_path` points to a folder whose `*.fastq*` files are all used for that sample (multiple lanes or both PE mates together). Synthetic junction alignment uses Bowtie2 in unpaired (`-U`) mode; junction tables may include a “paired” category when mate pairs map with opposite orientations on the synthetic sequence—single-end data will usually show zero there without indicating an error.

## Flow Diagram

```
{xxx}  input
xxx[]  array
- - -  optional
(?)    optional

                                      INPUTS
   ┌────────────────────────────────────┬────────────────────────────────────┐
{FASTQ}                             {ref_name}                          {anc_path}
   │                                    │                               (optional)
   │                                    ▼                                    │
   │                          ┌───────────────────┐                          │
   │                          │ 1. GetRefGenome   │                          │
   │                          └─────────┬─────────┘                          │
   │                                    ▼                                    │
   │           ┌──────────────────    Genome    ──────────────┐              │
   │           │                  (FASTA + GBK)               │              │
   │           │                        │                     │              │
   │           │                        ▼                     │              │
   │           │          ┌────────────────────────────┐      │              │
   │           │          │ 2a. LocateTNsUsingGenbank  │      │              │
   │           │          │ 2b. LocateTNsUsingISfinder │      │              │
   │           │          │ 2c. LocateTNsUsingISEScan  │      │              │
   │           │          └─────────────┬──────────────┘      │              │
   │           │                        ▼                     │              │
   │           │                     RefTn[]                  │              │
   │           │                        │                     │              │
   │           │                        ▼                     │              │
   │           │           ┌────────────────────────┐         │              │
   │           │           │  3. CreateRefTnJc      │         │              │
   │           │           └────────────┬───────────┘         │              │
   │           │                        ▼                     │              │
   │           │             ┌── RefTnJunction[]              │              │
   │           │             │                                │              │
   │           ▼             │                                │              │
   │     ┌───────────┐       │                                │              │
   │────►│ 4. Breseq │       │                                │              │
   │     └─────┬─────┘       │                                │              │
   │           ▼             │                                │              │
   │     BreseqJunction[]     │                                │              │
   │      + coverage         ▼                                │              │
   │            └─────┬──────┘                                │              │
   │                  ▼            ┌───────────────────┐      │              │
   │               all Jc ────────►│  5. CreateTnJc    │◄─────┘              │
   │                               └─────────┬─────────┘                     │
   │                                         ▼                               │
   │                                    TnJunction[]                         │
   │                                         │                               │
   │                                         ▼                               │
   │                            ┌─────────────────────────┐                  │
   │                            │  6. PairTnJcToRawTnJc2  │                  │
   │                            └────────────┬────────────┘                  │
   │                                         ▼                               │
   │                                     RawTnJc2[]                          │
   │                                         │                               │
   │                                         ▼                               │
   │                     ┌──────────────────────────────────────┐            │
   │                     │ 7. LinkTnJc2ToSingleLocusPairs       │            │
   │                     └──────────────────┬───────────────────┘            │
   │                                        ▼                                │
   │                                 SingleLocusLinkedTnJc2[]                │
   │                                        │                                │
   │                                        ▼                                │
   │                        ┌───────────────────────────────┐                │
   │                        │ 8. CalcTnJc2AmpliconCoverage  │◄ - - - - - - - ┤
   │                        └───────────────┬───────────────┘     anc        │
   │                                        ▼                  coverage      │
   │                                   CoveredTnJc2[]                        │
   │                                        │                                │
   │                                        ▼                                │
   │                          ┌───────────────────────────┐                  │
   │                          │ 9. FilterTnJc2Candidates  │                  │
   │                          └─────────────┬─────────────┘                  │
   │                                        ▼                                │
   │                                 CoveredTnJc2[] (filtered)               │
   │                                        │                                │
   │     ┌──────────────────┐               │                                │
   │────►│  10. ReadLen     │               │                                │
   │     └────────┬─────────┘               │                                │
   │              ▼                         │                                │
   │         ReadLengths                    │                                │
   │              │                         ▼                                │
   │              │       ┌──────────────────────────────┐                   │
   │              └──────►│ 11. CreateSyntheticJunctions │                   │
   │                         └──────────────┬───────────────┘                │
   │                                        ▼                                │
   │                                 SynJctsTnJc2[]                          │
   │                          (7 JunctionTypes per candidate)                │
   │                                        │                                │
   │                                        ▼                                │
   │                     ┌──────────────────────────────────────┐            │
   └────────────────────►│ 12. AlignReadsToJunctions            │◄- - - - - -┘
            {FASTQ}      └──────────────────┬───────────────────┘ {anc_FASTQ}
                                            ▼
                                     iso.bam, anc.bam(?)
                                            │
                                            ▼
                             ┌─────────────────────────────┐
                             │ 13. AnalyzeTnJc2Alignments  │
                             └───────────────┬─────────────┘
                                             ▼
                                  AnalyzedTnJc2[] (junction coverage)
                                             │
                                             ▼
                             ┌──────────────────────────────┐
                             │ 14. ClassifyTnJc2Candidates  │
                             └───────────────┬──────────────┘
                                             ▼
                                   ClassifiedTnJc2[] (event + modifiers)
                                             │
                                             ▼
                             ┌──────────────────────────────┐
                             │ 15. PlotTnJc2Coverage        │
                             └───────────────┬──────────────┘
                                             ▼
                                      coverage plots (optional)
                                             │
                                             ▼
                                   ┌───────────────────┐
                                   │ 16. ExportTnJc2   │
                                   └─────────┬─────────┘
                                             ▼
                                        summary.yaml
                                        summary.json
```

## Steps

| # | Step | Description |
|---|------|-------------|
| 1 | GetRefGenome | Download reference FASTA and GenBank from NCBI |
| 2a | LocateTNsUsingGenbank | Find IS elements from GenBank annotations |
| 2b | LocateTNsUsingISfinder | Find IS elements via BLAST against ISfinder DB |
| 2c | LocateTNsUsingISEScan | Find IS elements via ISEScan |
| 3 | CreateRefTnJc | Build reference junctions at IS element boundaries |
| 4 | Breseq | Align reads and detect novel junctions |
| 5 | CreateTnJc | Match breseq junctions to IS elements |
| 6 | PairTnJcToRawTnJc2 | Pair left/right junctions into amplicon candidates |
| 7 | LinkTnJc2ToSingleLocusPairs | Classify junction pair structures (flanked, unflanked, etc.) |
| 8 | CalcTnJc2AmpliconCoverage | Calculate amplicon coverage (raw or normalized by ancestor) |
| 9 | FilterTnJc2Candidates | Filter by copy number and amplicon length |
| 10 | ReadLen | Determine read lengths from FASTQ and calculate junction arm lengths |
| 11 | CreateSyntheticJunctions | Build 7 synthetic junction sequences per candidate |
| 12 | AlignReadsToJunctions | Align isolate (and ancestor) reads to synthetic junctions |
| 13 | AnalyzeTnJc2Alignments | Count spanning/left/right reads at each junction |
| 14 | ClassifyTnJc2Candidates | Classify events from read patterns (de novo, ancestral, etc.) |
| 15 | PlotTnJc2Coverage | Generate junction and amplicon coverage plots |
| 16 | ExportTnJc2 | Write summary.yaml and summary.json |

## Coverage Modes

- **Without ancestor** (`-a` not provided): raw coverage depth from breseq.
- **With ancestor** (`-a` provided): isolate/ancestor coverage ratio. Ancestor breseq output and junction alignments are cached and shared across isolates.

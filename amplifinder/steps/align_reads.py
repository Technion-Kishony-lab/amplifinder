"""Step 11: Align reads to synthetic junctions."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2
from amplifinder.steps.base import Step
from amplifinder.tools.bowtie2 import align_reads_to_fasta


class AlignReadsToJunctionsStep(Step):
    """Align reads to synthetic junction sequences."""
    is_ancestor: bool = False

    def __init__(
        self,
        synjcs_tnjc2s: RecordTypedDf[SynJctsTnJc2],
        output_dir: Path,
        fastq_path: Path,
        threads: int = 1,
        score_min: Optional[str] = None,  # None = use default (G,0,-0.25 for local)
        num_alignments: int = 100,  # TODO: should this be '1' ?
        force: Optional[bool] = None,
    ):
        self.synjcs_tnjc2s = synjcs_tnjc2s
        self.output_dir = Path(output_dir)
        self.fastq_path = Path(fastq_path)
        self.threads = threads
        self.score_min = score_min
        self.num_alignments = num_alignments

        # Build list of expected output BAM files
        output_files = []
        input_files = [fastq_path]
        
        for filtered_tnjc2 in synjcs_tnjc2s:
            # Junction FASTA file is input
            input_files.append(filtered_tnjc2.fasta_path(self.output_dir, is_ancestor=self.is_ancestor))
            # BAM file is output
            output_files.append(filtered_tnjc2.bam_path(self.output_dir, is_ancestor=self.is_ancestor))

        super().__init__(
            input_files=input_files,
            artifact_files=output_files,
            force=force,
        )

    def _artifact_labels(self) -> list[str]:
        """Summarize outputs as count."""
        n = len(self.synjcs_tnjc2s)
        return [f"{n} BAM files of synthetic junctions"]

    def _generate_artifacts(self) -> None:
        """Align reads to synthetic junctions; skip BAMs that already exist."""
        analysis_dir_names = [tnjc2.analysis_dir_name(is_ancestor=self.is_ancestor) for tnjc2 in self.synjcs_tnjc2s]
        max_name_length = max(len(name) for name in analysis_dir_names)
        for filtered_tnjc2 in self.synjcs_tnjc2s:
            junctions_fasta = filtered_tnjc2.fasta_path(self.output_dir, is_ancestor=self.is_ancestor)
            assert junctions_fasta.exists()

            bam_path = filtered_tnjc2.bam_path(self.output_dir, is_ancestor=self.is_ancestor)
            print(f"{filtered_tnjc2.analysis_dir_name(is_ancestor=self.is_ancestor):<{max_name_length}}: ", end="", flush=True)
            if bam_path.exists():
                print("file exists, skipping")
                continue
            align_reads_to_fasta(
                ref_fasta=junctions_fasta,
                fastq_path=self.fastq_path,
                output_bam=bam_path,
                threads=self.threads,
                score_min=self.score_min,
                num_alignments=self.num_alignments,
            )
            assert bam_path.exists()


class AncAlignReadsToJunctionsStep(AlignReadsToJunctionsStep):
    """Align ancestor reads to synthetic junctions."""
    is_ancestor: bool = True

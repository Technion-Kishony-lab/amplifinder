"""Step 11: Align reads to synthetic junctions."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, SynJctsTnJc2
from amplifinder.config import BowtieParams
from amplifinder.logger import logger
from amplifinder.steps.base import Step
from amplifinder.tools.bowtie2 import align_reads_to_fasta
from amplifinder.utils.file_lock import locked_resource


class AlignReadsToJunctionsStep(Step):
    """Align reads to synthetic junction sequences."""
    is_ancestor: bool = False

    def __init__(
        self,
        synjcs_tnjc2s: RecordTypedDf[SynJctsTnJc2],
        output_dir: Path,
        fastq_path: Path,
        threads: int = 4,
        bowtie_params: BowtieParams = None,
        force: Optional[bool] = None,
    ):
        self.synjcs_tnjc2s = synjcs_tnjc2s
        self.output_dir = Path(output_dir)
        self.fastq_path = Path(fastq_path)
        self.threads = threads
        self.bowtie_params = bowtie_params or BowtieParams()

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
            name = filtered_tnjc2.analysis_dir_name(is_ancestor=self.is_ancestor)
            logger.print_progress(f"{name:<{max_name_length}}: ", end="")

            # Lock per-junction for ancestor (None for isolate = no lock)
            lock_path = bam_path.parent if self.is_ancestor else None

            with locked_resource(lock_path, "junction_bam", timeout=1800):
                if bam_path.exists():
                    logger.print_progress("file exists, skipping")
                    continue
                align_reads_to_fasta(
                    ref_fasta=junctions_fasta,
                    fastq_path=self.fastq_path,
                    output_bam=bam_path,
                    threads=self.threads,
                    score_min=self.bowtie_params.score_min,
                    mismatch_penalty=self.bowtie_params.mismatch_penalty,
                    num_alignments=self.bowtie_params.num_alignments,
                    local=self.bowtie_params.local,
                    min_qlen=self.bowtie_params.min_qlen,
                )
                assert bam_path.exists()


class AncAlignReadsToJunctionsStep(AlignReadsToJunctionsStep):
    """Align ancestor reads to synthetic junctions."""
    is_ancestor: bool = True

"""Step 11: Align reads to synthetic junctions."""

from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDF, CandidateTnJc2
from amplifinder.steps.base import Step
from amplifinder.tools.bowtie2 import align_reads_to_fasta
from amplifinder.logger import info


class AlignReadsToJunctionsStep(Step[RecordTypedDF[CandidateTnJc2]]):
    """Align reads to synthetic junction sequences.
    
    Alignment depends on run type:
    - anc_fastq_path=None: align isolate reads only → iso.sorted.bam
    - anc_fastq_path=set: align both isolate and ancestor reads → iso.sorted.bam + anc.sorted.bam
    """

    def __init__(
        self,
        candidates: RecordTypedDF[CandidateTnJc2],
        output_dir: Path,
        iso_fastq_path: Path,
        anc_fastq_path: Optional[Path] = None,
        threads: int = 1,
        score_min: str = "L,0,-0.6",
        num_alignments: int = 10,
        force: Optional[bool] = None,
    ):
        self.candidates = candidates
        self.output_dir = Path(output_dir)
        self.iso_fastq_path = Path(iso_fastq_path)
        self.anc_fastq_path = Path(anc_fastq_path) if anc_fastq_path else None
        self.threads = threads
        self.score_min = score_min
        self.num_alignments = num_alignments
        
        # Build list of expected output BAM files
        output_files = []
        for cand in candidates:
            analysis_dir = output_dir / cand.analysis_dir
            output_files.append(analysis_dir / "iso.sorted.bam")
            if anc_fastq_path:
                output_files.append(analysis_dir / "anc.sorted.bam")
        
        input_files = [iso_fastq_path]
        if anc_fastq_path:
            input_files.append(anc_fastq_path)
        # Junction FASTA files are inputs
        for cand in candidates:
            input_files.append(output_dir / cand.analysis_dir / "junctions.fasta")
        
        super().__init__(
            input_files=input_files,
            output_files=output_files,
            force=force,
        )

    @property
    def has_ancestor(self) -> bool:
        """True if ancestor reads should be aligned."""
        return self.anc_fastq_path is not None

    def _calculate_output(self) -> RecordTypedDF[CandidateTnJc2]:
        """Align reads to synthetic junctions for each candidate."""
        for cand in self.candidates:
            analysis_dir = self.output_dir / cand.analysis_dir
            junctions_fasta = analysis_dir / "junctions.fasta"
            
            if not junctions_fasta.exists():
                info(f"Skipping {cand.analysis_dir}: no junctions.fasta")
                continue
            
            # Align isolate reads
            iso_bam = analysis_dir / "iso.sorted.bam"
            align_reads_to_fasta(
                ref_fasta=junctions_fasta,
                fastq_path=self.iso_fastq_path,
                output_bam=iso_bam,
                threads=self.threads,
                score_min=self.score_min,
                num_alignments=self.num_alignments,
            )
            
            # Align ancestor reads if provided
            if self.has_ancestor:
                anc_bam = analysis_dir / "anc.sorted.bam"
                align_reads_to_fasta(
                    ref_fasta=junctions_fasta,
                    fastq_path=self.anc_fastq_path,
                    output_bam=anc_bam,
                    threads=self.threads,
                    score_min=self.score_min,
                    num_alignments=self.num_alignments,
                )
        
        return self.candidates

    def _save_output(self, output: RecordTypedDF[CandidateTnJc2]) -> None:
        """BAM files already created in _calculate_output."""
        pass

    def load_outputs(self) -> RecordTypedDF[CandidateTnJc2]:
        """Return candidates (BAM files are side effects)."""
        return self.candidates


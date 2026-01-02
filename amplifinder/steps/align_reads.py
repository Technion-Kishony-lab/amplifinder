"""Step 11: Align reads to synthetic junctions."""

import shutil
from pathlib import Path
from typing import Optional

from amplifinder.data_types import RecordTypedDf, FilteredTnJc2
from amplifinder.steps.base import Step
from amplifinder.tools.bowtie2 import align_reads_to_fasta
from amplifinder.utils import ensure_dir


class AlignReadsToJunctionsStep(Step[RecordTypedDf[FilteredTnJc2]]):
    """Align reads to synthetic junction sequences.

    Alignment depends on run type:
    - anc_fastq_path=None: align isolate reads only → iso.sorted.bam
    - anc_fastq_path=set: align both isolate and ancestor reads → iso.sorted.bam + anc.sorted.bam
    """

    def __init__(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        output_dir: Path,
        iso_fastq_path: Path,
        anc_fastq_path: Optional[Path] = None,
        threads: int = 1,
        score_min: Optional[str] = None,  # None = use default (G,0,-0.25 for local)
        num_alignments: int = 10,
        force: Optional[bool] = None,
    ):
        self.filtered_tnjc2s = filtered_tnjc2s
        self.output_dir = Path(output_dir)
        self.iso_fastq_path = Path(iso_fastq_path)
        self.anc_fastq_path = Path(anc_fastq_path) if anc_fastq_path else None
        self.threads = threads
        self.score_min = score_min
        self.num_alignments = num_alignments

        # Build list of expected output BAM files
        # Only process candidates with valid chosen_tn_id (those with junction files)
        output_files = []
        input_files = [iso_fastq_path]
        if anc_fastq_path:
            input_files.append(anc_fastq_path)
        
        for filtered_tnjc2 in filtered_tnjc2s:
            if filtered_tnjc2.chosen_tn_id is None:
                continue
            analysis_dir = output_dir / "junctions" / filtered_tnjc2.analysis_dir
            # Junction FASTA file is input
            input_files.append(analysis_dir / "junctions.fasta")
            # BAM files are outputs
            output_files.append(analysis_dir / "iso.sorted.bam")
            if anc_fastq_path:
                output_files.append(analysis_dir / "anc.sorted.bam")

        super().__init__(
            input_files=input_files,
            output_files=output_files,
            force=force,
        )

    def _output_labels(self) -> list[str]:
        """Summarize outputs as count."""
        if self.output_files:
            n = len(self.filtered_tnjc2s)
            return [f"{n} junctions"]
        return []

    @property
    def has_ancestor(self) -> bool:
        """True if ancestor reads should be aligned."""
        return self.anc_fastq_path is not None

    def _calculate_output(self) -> RecordTypedDf[FilteredTnJc2]:
        """Align reads to synthetic junctions for each candidate."""
        for filtered_tnjc2 in self.filtered_tnjc2s:
            analysis_dir = self.output_dir / "junctions" / filtered_tnjc2.analysis_dir
            junctions_fasta = analysis_dir / "junctions.fasta"

            if not junctions_fasta.exists():
                self.log(f"Skipping {filtered_tnjc2.analysis_dir}: no junctions.fasta")
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

        return self.filtered_tnjc2s

    def _save_output(self, output: RecordTypedDf[FilteredTnJc2]) -> None:
        """BAM files already created in _calculate_output."""
        pass

    def load_outputs(self) -> RecordTypedDf[FilteredTnJc2]:
        """Return candidates (BAM files are side effects)."""
        return self.filtered_tnjc2s


class AncAlignReadsToJunctionsStep(AlignReadsToJunctionsStep):
    """Align ancestor reads to synthetic junctions.
    
    Copies junction files from isolate folder to ancestor folder before alignment.
    """
    
    def __init__(self, iso_output_dir: Path, **kwargs):
        self.iso_output_dir = Path(iso_output_dir)
        super().__init__(**kwargs)
    
    def _calculate_output(self) -> RecordTypedDf[FilteredTnJc2]:
        """Copy junctions then align reads."""
        # Copy junctions first (only if not already there)
        for filtered_tnjc2 in self.filtered_tnjc2s:
            iso_jc_dir = self.iso_output_dir / "junctions" / filtered_tnjc2.analysis_dir
            anc_jc_dir = self.output_dir / "junctions" / filtered_tnjc2.analysis_dir
            anc_junctions = anc_jc_dir / "junctions.fasta"
            
            # Only copy if not already in ancestor folder
            if anc_junctions.exists():
                continue
            
            # Copy junctions.fasta from isolate to ancestor folder
            iso_junctions = iso_jc_dir / "junctions.fasta"
            if iso_junctions.exists():
                ensure_dir(anc_jc_dir)
                shutil.copy2(iso_junctions, anc_junctions)
        
        # Then perform alignment
        return super()._calculate_output()

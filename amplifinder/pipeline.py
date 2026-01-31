"""Pipeline orchestration for AmpliFinder."""
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Tuple, Optional

from amplifinder.config import Config
from amplifinder.env import BRESEQ_DOCKER, ISDB_PATH
from amplifinder.exceptions import PrematureTerminationError
from amplifinder.data_types import (
    Genome, RecordTypedDf, RefTn, RefTnJunction, BreseqJunction, TnJunction, RawTnJc2,
    CoveredTnJc2, SingleLocusLinkedTnJc2, SynJctsTnJc2, AnalyzedTnJc2, ClassifiedTnJc2,
)
from amplifinder.logger import logger, c, setup_logger
from amplifinder.utils.file_utils import ensure_dir

from amplifinder.steps import (
    InitializingStep,
    GetRefGenomeStep,
    LocateTNsUsingISfinderStep,
    LocateTNsUsingGenbankStep,
    BreseqStep, AncBreseqStep,
    CreateRefTnJcStep,
    CreateTnJcStep,
    PairTnJcToRawTnJc2Step,
    CalcTnJc2AmpliconCoverageStep,
    LinkTnJc2ToSingleLocusPairsStep,
    FilterTnJc2CandidatesStep,
    CreateSyntheticJunctionsStep, AncCreateSyntheticJunctionsStep,
    AlignReadsToJunctionsStep, AncAlignReadsToJunctionsStep,
    AnalyzeTnJc2AlignmentsStep, AncAnalyzeTnJc2AlignmentsStep, PlotTnJc2CoverageStep,
    ClassifyTnJc2CandidatesStep,
    ExportTnJc2Step,
    ReadLenStep, ReadLengths,
)
from amplifinder.data import get_builtin_isfinder_db_path
from amplifinder.steps.locate_tns import compare_tn_locations


@dataclass
class Pipeline:
    """AmpliFinder pipeline with phased execution."""

    config: Config
    verbose: bool = False

    def _setup_logger(self) -> None:
        """Set up logger with file in run directory."""
        ensure_dir(self.config.iso_run_dir)
        log_file = self.config.iso_run_dir / "run_log.txt"
        setup_logger(log_path=log_file, use_colors=True, verbose=self.verbose)

    def _calc_read_lengths(self) -> ReadLengths:
        """Calculate read lengths and junction lengths."""
        return ReadLenStep(
            iso_fastq_path=self.config.iso_fastq_path,
            anc_fastq_path=self.config.anc_fastq_path if self.config.has_ancestor else None,
            iso_read_length=self.config.iso_read_length,
            anc_read_length=self.config.anc_read_length,
            min_num_bases=self.config.min_num_bases,
        ).run()

    def _log_run_info(self, include_no_ancestor_warning: bool = True) -> None:
        """Log run information: reference, isolate, and ancestor."""
        if self.config.anc_fastq_path:
            anc_msg = f"Ancestor: {self.config.anc_fastq_path}"
        elif include_no_ancestor_warning:
            anc_msg = "No ancestor assigned; using raw (non-normalized) coverage analysis"
        else:
            anc_msg = "No ancestor"
        logger.info(
            f"\nReference: {c(self.config.ref_name, 'cyan')}\n"
            f"Isolate: {c(str(self.config.iso_fastq_path), 'magenta')}\n{anc_msg}")

    def run(self) -> Optional[RecordTypedDf[ClassifiedTnJc2]]:
        """Run full pipeline with exception handling and status tracking."""
        self._setup_logger()
        self._log_run_info()

        # Initialize output directories
        iso_output, anc_output = self._initialize()

        try:
            # Clear old status files and mark start
            self._clear_status_files(iso_output)
            self._write_status_file(iso_output, 'started')

            # Run the pipeline
            classified_tnjc2s = self._run(iso_output, anc_output)

            # Mark successful completion
            self._write_status_file(iso_output, 'completed')
            return classified_tnjc2s

        except PrematureTerminationError as e:
            # Graceful early termination - not an error
            self._write_status_file(iso_output, 'terminated', e.get_detailed_message())
            logger.info(f"Pipeline terminated early:\n{e}")
            return None  # Graceful exit

        except Exception as e:
            # Actual errors
            self._write_status_file(iso_output, 'failed', str(e))
            raise

    def _run(self, iso_output: Path, anc_output: Optional[Path]) -> RecordTypedDf[ClassifiedTnJc2]:
        """Execute pipeline steps, return classified candidates."""

        # Load reference genome (needed for ancestor breseq and isolate pipeline)
        genome = self._load_reference()

        # Run pipeline
        ref_tns = self._locate_tns_in_reference(genome)
        ref_tnjcs = self._create_ref_tn_junctions(ref_tns, genome, iso_output)
        breseq_jcs = self._run_breseq(genome)
        tnjcs = self._create_tnjcs(breseq_jcs, ref_tnjcs, genome, iso_output)
        raw_tnjc2s = self._create_tnjc2s(tnjcs, genome, iso_output)
        covered_tnjc2s = self._calc_amplicon_coverage(raw_tnjc2s, genome, iso_output)
        linked_tnjc2s = self._link_side_of_tnjc2s_to_single_locus_pairs(covered_tnjc2s, genome, ref_tns, iso_output)
        filtered_tnjc2s = self._filter_candidates(linked_tnjc2s, iso_output)
        read_lengths = self._calc_read_lengths()
        synjct_tnjc2s = self._create_synthetic_junctions(
            filtered_tnjc2s, genome, ref_tns, iso_output, anc_output, read_lengths)
        self._align_reads(synjct_tnjc2s, iso_output, anc_output)
        analyzed_tnjc2s, iso_alignment_cache, anc_alignment_cache = self._analyze_alignments(
            synjct_tnjc2s, iso_output, anc_output, read_lengths)
        classified_tnjc2s = self._classify_candidates(analyzed_tnjc2s, iso_output)
        self._plot_coverage(classified_tnjc2s, iso_output, anc_output, read_lengths,
                            iso_alignment_cache, anc_alignment_cache)
        self._export(classified_tnjc2s, ref_tns, iso_output, read_lengths)

        return classified_tnjc2s

    def run_breseq_only(self) -> None:
        """Run only breseq steps (ancestor and isolate), then exit."""
        self._setup_logger()
        self._log_run_info(include_no_ancestor_warning=False)

        # Initialize output directories
        iso_output, anc_output = self._initialize()

        genome = self._load_reference()
        self._run_breseq(genome)
        logger.info("breseq-only mode completed")

    @property
    def ref_tn_source(self) -> str:
        """Get the source of the reference TN data."""
        return self.config.use_isfinder and "isfinder" or "genbank"

    def _write_status_file(self, iso_output: Path, status: str, reason: str = "") -> None:
        """Write status marker file to isolate run directory.

        Args:
            iso_output: Isolate run directory
            status: Status name (started, completed, terminated, failed)
            reason: Optional reason/message to include in the status file
        """
        status_file = iso_output / f"run.{status}"
        with open(status_file, 'w') as f:
            f.write(f"timestamp: {datetime.now().isoformat()}\n")
            f.write(f"reference: {self.config.ref_name}\n")
            f.write(f"isolate: {self.config.iso_name}\n")
            f.write(f"ancestor: {self.config.anc_name}\n")
            if reason:
                f.write(f"reason: {reason}\n")

    def _clear_status_files(self, iso_output: Path) -> None:
        """Remove all status marker files from isolate run directory."""
        for status in ['started', 'completed', 'terminated', 'failed']:
            status_file = iso_output / f"run.{status}"
            if status_file.exists():
                status_file.unlink()

    def _initialize(self) -> Tuple[Path, Optional[Path]]:
        """Step 0: Initialize output directories."""
        return InitializingStep(config=self.config).run()

    def _load_reference(self) -> Genome:
        """Step 1: Get reference genome."""
        return GetRefGenomeStep(
            ref_name=self.config.ref_name,
            ref_path=self.config.ref_path,
            ncbi=self.config.ncbi,
        ).run()

    def _locate_tns_in_reference(self, genome: Genome) -> RecordTypedDf[RefTn]:
        """Step 2: Locate TN elements using GenBank and/or ISfinder."""
        cfg = self.config

        tn_loc_dir = cfg.ref_path / "tn_loc" / genome.name

        # 2a: GenBank annotations
        ref_tns_genbank = LocateTNsUsingGenbankStep(
            genome=genome,
            output_dir=tn_loc_dir,
        ).run()

        # 2b: ISfinder database
        ref_tns_isfinder = LocateTNsUsingISfinderStep(
            genome=genome,
            output_dir=tn_loc_dir,
            isdb_path=ISDB_PATH or get_builtin_isfinder_db_path(),
            evalue=cfg.isfinder_evalue,
            critical_coverage=cfg.isfinder_critical_coverage,
        ).run()

        # Compare sources
        if ref_tns_genbank is not None:
            compare_tn_locations(ref_tns_genbank, ref_tns_isfinder, output_file=tn_loc_dir / "tn_location_diffs.txt")

        # Select source
        if ref_tns_genbank is None and not cfg.use_isfinder:
            raise ValueError("No TN locations - provide GenBank or set --use-isfinder")
        tn_loc = ref_tns_isfinder if cfg.use_isfinder else ref_tns_genbank
        assert tn_loc is not None
        return tn_loc

    def _create_ref_tn_junctions(
        self, ref_tn_locs: RecordTypedDf[RefTn], genome: Genome, iso_output: Path
    ) -> RecordTypedDf[RefTnJunction]:
        """Step 3: Create reference junctions."""
        cfg = self.config

        return CreateRefTnJcStep(
            ref_tn_locs=ref_tn_locs,
            genome=genome,
            output_dir=iso_output,
            source=self.ref_tn_source,
            reference_IS_out_span=cfg.reference_IS_out_span,
            reference_IS_in_span=cfg.reference_IS_in_span,
        ).run()

    def _run_breseq(self, genome: Genome) -> RecordTypedDf[BreseqJunction]:
        """Step 4: Run breseq on isolate to get junctions."""
        cfg = self.config

        if cfg.has_ancestor:
            AncBreseqStep(
                fastq_path=cfg.anc_fastq_path,
                ref_file=genome.genbank_path or genome.fasta_path,
                breseq_path=cfg.get_anc_breseq_path(),
                docker=BRESEQ_DOCKER,
                threads=cfg.threads,
                breseq_output_size_threshold=cfg.breseq_output_size_threshold,
                sample_name=f"ancestor ({cfg.anc_name})",
                remove_jc_breseq_reject=cfg.remove_jc_breseq_reject,
            ).run()

        return BreseqStep(
            fastq_path=cfg.iso_fastq_path,
            ref_file=genome.genbank_path or genome.fasta_path,
            breseq_path=cfg.get_iso_breseq_path(),
            docker=BRESEQ_DOCKER,
            threads=cfg.threads,
            breseq_output_size_threshold=cfg.breseq_output_size_threshold,
            sample_name=f"isolate ({cfg.iso_name})",
            remove_jc_breseq_reject=cfg.remove_jc_breseq_reject,
        ).run()

    def _create_tnjcs(
        self,
        breseq_jcs: RecordTypedDf[BreseqJunction],
        ref_tnjcs: RecordTypedDf[RefTnJunction],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDf[TnJunction]:
        """Step 5: Match junctions to TN elements."""
        cfg = self.config

        # Combine breseq junctions with reference TN junctions as a list to preserve types
        junctions = list(ref_tnjcs) + list(breseq_jcs)

        return CreateTnJcStep(
            junctions=junctions,
            ref_tnjcs=ref_tnjcs,
            genome=genome,
            output_dir=iso_output,
            max_dist_to_tn=cfg.max_dist_to_IS,
            trim_jc_flanking=cfg.trim_jc_flanking,
        ).run()

    def _create_tnjc2s(
        self,
        tnjcs: RecordTypedDf[TnJunction],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDf[RawTnJc2]:
        """Step 6: Combine junction pairs (RawTnJc2)."""
        return PairTnJcToRawTnJc2Step(
            tnjcs=tnjcs,
            genome=genome,
            output_dir=iso_output,
        ).run()

    def _calc_amplicon_coverage(
        self,
        raw_tnjc2s: RecordTypedDf[RawTnJc2],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDf[CoveredTnJc2]:
        """Step 7: Calculate amplicon coverage."""
        cfg = self.config
        iso_breseq_path, anc_breseq_path = cfg.get_breseq_paths()

        return CalcTnJc2AmpliconCoverageStep(
            raw_tnjc2s=raw_tnjc2s,
            output_dir=iso_output,
            iso_breseq_path=iso_breseq_path,
            anc_breseq_path=anc_breseq_path,
            average_method=cfg.average_method,
            min_amplicon_length=cfg.min_amplicon_length,
            max_amplicon_length=cfg.max_amplicon_length,
        ).run()

    def _link_side_of_tnjc2s_to_single_locus_pairs(
        self,
        covered_tnjc2s: RecordTypedDf[CoveredTnJc2],
        genome: Genome,
        tn_loc: RecordTypedDf[RefTn],
        iso_output: Path,
    ) -> RecordTypedDf[SingleLocusLinkedTnJc2]:
        """Step 8: Classify junction pair structures."""
        return LinkTnJc2ToSingleLocusPairsStep(
            covered_tnjc2s=covered_tnjc2s,
            genome=genome,
            tn_locs=tn_loc,
            output_dir=iso_output,
            transposition_threshold=self.config.min_amplicon_length,
        ).run()

    def _filter_candidates(
        self,
        linked_tnjc2s: RecordTypedDf[SingleLocusLinkedTnJc2],
        iso_output: Path,
    ) -> RecordTypedDf[SingleLocusLinkedTnJc2]:
        """Step 9: Filter candidates by amplicon length and copy number."""
        return FilterTnJc2CandidatesStep(
            linked_tnjc2s=linked_tnjc2s,
            output_dir=iso_output,
            min_amplicon_length=self.config.min_amplicon_length,
            max_amplicon_length=self.config.max_amplicon_length,
            replication_copy_number_threshold=self.config.replication_copy_number_threshold,
            deletion_copy_number_threshold=self.config.deletion_copy_number_threshold,
        ).run()

    def _create_synthetic_junctions(
        self,
        filtered_tnjc2s: RecordTypedDf[SingleLocusLinkedTnJc2],
        genome: Genome,
        ref_tns: RecordTypedDf[RefTn],
        iso_output: Path,
        anc_output: Optional[Path],
        read_lengths: ReadLengths,
    ) -> RecordTypedDf[SynJctsTnJc2]:
        """Step 10: Create synthetic junction sequences."""

        # Create junctions for isolate
        syn_tnjc2s = CreateSyntheticJunctionsStep(
            filtered_tnjc2s=filtered_tnjc2s,
            genome=genome,
            ref_tns=ref_tns,
            output_dir=iso_output,
            jc_arm_len=read_lengths.jc_arm_len_iso,
        ).run()

        # Create junctions for ancestor if needed
        if anc_output:
            syn_tnjc2s = AncCreateSyntheticJunctionsStep(
                filtered_tnjc2s=syn_tnjc2s,
                genome=genome,
                ref_tns=ref_tns,
                output_dir=anc_output,
                jc_arm_len=read_lengths.jc_arm_len_anc,
            ).run()

        return syn_tnjc2s

    def _align_reads(
        self,
        syn_tnjc2s: RecordTypedDf[SynJctsTnJc2],
        iso_output: Path,
        anc_output: Optional[Path],
    ) -> None:
        """Step 11: Align reads to synthetic junctions."""
        cfg = self.config

        # Align isolate reads in isolate folder
        AlignReadsToJunctionsStep(
            synjcs_tnjc2s=syn_tnjc2s,
            output_dir=iso_output,
            fastq_path=cfg.iso_fastq_path,
            threads=cfg.threads,
            bowtie_params=cfg.bowtie_params,
        ).run()

        # Align ancestor reads in ancestor folder
        if anc_output:
            AncAlignReadsToJunctionsStep(
                synjcs_tnjc2s=syn_tnjc2s,
                output_dir=anc_output,
                fastq_path=cfg.anc_fastq_path,
                threads=cfg.threads,
                bowtie_params=cfg.bowtie_params,
            ).run()

    def _analyze_alignments(
        self,
        synjct_tnjc2s: RecordTypedDf[SynJctsTnJc2],
        iso_output: Path,
        anc_output: Optional[Path],
        read_lengths: ReadLengths,
    ) -> tuple[RecordTypedDf[AnalyzedTnJc2], dict, dict]:
        """Step 12: Analyze read alignments.

        Returns:
            Tuple of (analyzed_tnjc2s, iso_alignment_cache, anc_alignment_cache)
        """
        cfg = self.config

        # Analyze isolate alignments
        analyzed_tnjc2s, iso_alignment_cache = AnalyzeTnJc2AlignmentsStep(
            tnjc2s=synjct_tnjc2s,
            output_dir=iso_output,
            arm_len=read_lengths.jc_arm_len_iso,
            read_length=read_lengths.read_len_iso,
            alignment_classify_params=cfg.alignment_analysis_params,
            alignment_filter_params=cfg.alignment_filter_params,
            jc_call_params=cfg.jc_call_params,
        ).run_with_cache()

        # Analyze ancestor alignments if present
        anc_alignment_cache = {}
        if anc_output:
            analyzed_tnjc2s, anc_alignment_cache = AncAnalyzeTnJc2AlignmentsStep(
                tnjc2s=analyzed_tnjc2s,
                output_dir=anc_output,
                arm_len=read_lengths.jc_arm_len_anc,
                read_length=read_lengths.read_len_anc,
                alignment_classify_params=cfg.alignment_analysis_params,
                alignment_filter_params=cfg.alignment_filter_params,
                jc_call_params=cfg.jc_call_params,
            ).run_with_cache()

        return analyzed_tnjc2s, iso_alignment_cache, anc_alignment_cache

    def _classify_candidates(
        self,
        analyzed_tnjc2s: RecordTypedDf[AnalyzedTnJc2],
        iso_output: Path,
    ) -> RecordTypedDf[ClassifiedTnJc2]:
        """Step 13: Final classification of candidates."""
        return ClassifyTnJc2CandidatesStep(
            analyzed_tnjc2s=analyzed_tnjc2s,
            output_dir=iso_output,
        ).run()

    def _plot_coverage(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        iso_output: Path,
        anc_output: Optional[Path],
        read_lengths: ReadLengths,
        iso_alignment_cache: dict,
        anc_alignment_cache: dict,
    ) -> None:
        """Plot junction/amplicon coverage (post-classification)."""
        cfg = self.config
        if not cfg.create_plots:
            return
        iso_breseq_path, anc_breseq_path = cfg.get_breseq_paths()

        PlotTnJc2CoverageStep(
            classified_tnjc2s=classified_tnjc2s,
            output_dir=iso_output,
            iso_breseq_path=iso_breseq_path,
            read_lengths=read_lengths,
            anc_output_dir=anc_output,
            anc_breseq_path=anc_breseq_path,
            alignment_classify_params=cfg.alignment_analysis_params,
            alignment_filter_params=cfg.alignment_filter_params,
            iso_alignments=iso_alignment_cache,
            anc_alignments=anc_alignment_cache,
        ).run()

    def _export(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        ref_tns: RecordTypedDf[RefTn],
        iso_output: Path,
        read_lengths: ReadLengths,
    ) -> None:
        """Step 14: Export results to YAML."""
        ExportTnJc2Step(
            classified_tnjc2s=classified_tnjc2s,
            output_dir=iso_output,
            ref_name=self.config.ref_name,
            iso_name=self.config.iso_name,
            read_lengths=read_lengths,
            ref_tns=ref_tns,
            anc_name=self.config.anc_name,
        ).run()


def run_pipeline(config: Config, verbose: bool = False) -> Optional[RecordTypedDf[ClassifiedTnJc2]]:
    """Run the AmpliFinder pipeline."""
    return Pipeline(config, verbose=verbose).run()

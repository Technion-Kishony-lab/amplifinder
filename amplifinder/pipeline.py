"""Pipeline orchestration for AmpliFinder."""
import shutil
import pandas as pd
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional

from amplifinder.config import Config
from amplifinder.data_types import (
    Genome, RecordTypedDf, RefTnLoc, RefTnJunction, Junction, BreseqJunction, TnJunction, RawTnJc2,
    CoveredTnJc2, ClassifiedTnJc2, FilteredTnJc2, AnalyzedTnJc2,
)
from amplifinder.logger import info
from amplifinder.utils.file_utils import ensure_dir
from amplifinder.utils.fasta import get_read_length
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
    ClassifyTnJc2StructureStep,
    FilterTnJc2CandidatesStep,
    CreateSyntheticJunctionsStep,
    AlignReadsToJunctionsStep, AncAlignReadsToJunctionsStep,
    AnalyzeTnJc2AlignmentsStep,
    ClassifyTnJc2CandidatesStep,
    ExportTnJc2Step,
)
from amplifinder.data import get_builtin_isfinder_db_path
from amplifinder.steps.locate_tns import compare_tn_locations


@dataclass
class Pipeline:
    """AmpliFinder pipeline with phased execution."""

    config: Config
    _iso_read_length: Optional[int] = None
    _anc_read_length: Optional[int] = None

    def _get_iso_read_length(self) -> int:
        """Get isolate read length from config or auto-detect from FASTQ."""
        if self._iso_read_length is None:
            self._iso_read_length = get_read_length(
                fastq_dir=self.config.iso_path,
                provided_length=self.config.iso_read_length,
            )
        return self._iso_read_length

    def _get_anc_read_length(self) -> int:
        """Get ancestor read length from config or auto-detect from FASTQ."""
        if self._anc_read_length is None:
            self._anc_read_length = get_read_length(
                fastq_dir=self.config.anc_path,
                provided_length=self.config.anc_read_length,
            )
        return self._anc_read_length

    def run(self) -> RecordTypedDf[AnalyzedTnJc2]:
        """Run full pipeline, return analyzed candidates."""

        info(
            f"Running AmpliFinder pipeline, reference: {self.config.ref_name}, "
            f"isolate: {self.config.iso_name}, ancestor: {self.config.anc_name}\n"
        )
        iso_output, anc_output = self._initialize()

        # Load reference genome (needed for ancestor breseq and isolate pipeline)
        genome = self._load_reference()

        # Handle ancestor breseq if needed
        self._ancestor_breseq(genome, anc_output)

        # Run isolate pipeline
        tn_loc = self._locate_tns_in_reference(genome)
        ref_tnjc = self._create_ref_tn_junctions(tn_loc, genome, iso_output)
        breseq_jc = self._run_breseq(genome, iso_output)
        tnjcs = self._create_tnjc(breseq_jc, ref_tnjc, genome, iso_output)
        raw_tnjc2s = self._create_tnjc2(tnjcs, genome, iso_output)
        covered_tnjc2s = self._calc_amplicon_coverage(raw_tnjc2s, genome, iso_output)
        classified_tnjc2s = self._classify_structure(covered_tnjc2s, genome, tn_loc, iso_output)
        filtered_tnjc2s = self._filter_candidates(classified_tnjc2s, iso_output)
        self._create_synthetic_junctions(filtered_tnjc2s, genome, tn_loc, iso_output)
        self._align_reads(filtered_tnjc2s, iso_output)
        analyzed_tnjc2s = self._analyze_alignments(filtered_tnjc2s, iso_output, anc_output)
        analyzed_tnjc2s = self._classify_candidates(analyzed_tnjc2s, iso_output)
        self._export(analyzed_tnjc2s, genome, iso_output)

        return analyzed_tnjc2s

    def run_breseq_only(self) -> None:
        """Run only breseq steps (ancestor and isolate), then exit."""
        iso_output, anc_output = self._initialize()
        genome = self._load_reference()
        self._ancestor_breseq(genome, anc_output)
        self._run_breseq(genome, iso_output)
        info("breseq-only mode completed")

    def _ancestor_breseq(self, genome: Genome, anc_output: Path) -> None:
        """Ensure ancestor breseq output exists, run breseq if needed."""
        if not self.config.has_ancestor:
            return

        # Determine breseq path (provided or default)
        if self.config.anc_breseq_path:
            anc_breseq_path = Path(self.config.anc_breseq_path)
        else:
            anc_breseq_path = anc_output / "breseq"

        # Run breseq (BreseqStep handles caching - skips if output exists)
        AncBreseqStep(
            output_path=anc_breseq_path,
            fastq_path=self.config.anc_path,
            ref_file=genome.genbank_path or genome.fasta_path,
            docker=self.config.breseq_docker,
            threads=self.config.breseq_threads,
        ).run()

    @property
    def ref_tn_source(self) -> str:
        """Get the source of the reference TN data."""
        return self.config.use_isfinder and "isfinder" or "genbank"

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

    def _locate_tns_in_reference(self, genome: Genome) -> RecordTypedDf[RefTnLoc]:
        """Step 2: Locate TN elements using GenBank and/or ISfinder."""
        cfg = self.config

        tn_loc_dir = cfg.ref_path / "tn_loc" / genome.name

        # 2a: GenBank annotations
        tn_loc_genbank = LocateTNsUsingGenbankStep(
            genome=genome,
            output_dir=tn_loc_dir,
        ).run()

        # 2b: ISfinder database
        tn_loc_isfinder = LocateTNsUsingISfinderStep(
            genome=genome,
            output_dir=tn_loc_dir,
            isdb_path=cfg.isdb_path or get_builtin_isfinder_db_path(),
            evalue=cfg.isfinder_evalue,
            critical_coverage=cfg.isfinder_critical_coverage,
        ).run()

        # Compare sources
        if tn_loc_genbank is not None:
            compare_tn_locations(tn_loc_genbank, tn_loc_isfinder, output_file=tn_loc_dir / "tn_location_diffs.txt")

        # Select source
        if tn_loc_genbank is None and not cfg.use_isfinder:
            raise ValueError("No TN locations - provide GenBank or set --use-isfinder")
        tn_loc = tn_loc_isfinder if cfg.use_isfinder else tn_loc_genbank
        assert tn_loc is not None
        return tn_loc

    def _create_ref_tn_junctions(
        self, ref_tn_locs: RecordTypedDf[RefTnLoc], genome: Genome, iso_output: Path
    ) -> RecordTypedDf[RefTnJunction]:
        """Step 3: Create reference junctions."""
        cfg = self.config

        return CreateRefTnJcStep(
            ref_tn_locs=ref_tn_locs,
            genome=genome,
            output_dir=iso_output,
            source=self.ref_tn_source,
            reference_tn_out_span=cfg.reference_IS_out_span,
        ).run()

    def _run_breseq(self, genome: Genome, iso_output: Path) -> RecordTypedDf[BreseqJunction]:
        """Step 4: Run breseq on isolate to get junctions."""
        cfg = self.config

        breseq_output = BreseqStep(
            fastq_path=cfg.iso_path,
            ref_file=genome.genbank_path or genome.fasta_path,
            output_path=cfg.iso_breseq_path or iso_output / "breseq",
            docker=cfg.breseq_docker,
            threads=cfg.breseq_threads,
        ).run()

        breseq_jc_df = breseq_output["JC"]
        breseq_jcs = RecordTypedDf(breseq_jc_df, BreseqJunction)
        return breseq_jcs

    def _create_tnjc(
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

    def _create_tnjc2(
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

        # Get ancestor breseq path if ancestor exists
        if cfg.has_ancestor:
            anc_breseq_path = cfg.get_anc_breseq_path()
            anc_name = cfg.anc_name
        else:
            anc_breseq_path = None
            anc_name = None

        iso_breseq_path = cfg.get_iso_breseq_path()

        return CalcTnJc2AmpliconCoverageStep(
            raw_tnjc2s=raw_tnjc2s,
            genome=genome,
            output_dir=iso_output,
            ref_name=cfg.ref_name,
            iso_breseq_path=iso_breseq_path,
            iso_name=cfg.iso_name,
            anc_breseq_path=anc_breseq_path,
            anc_name=anc_name,
            ncp_min=cfg.ncp_min,
            ncp_max=cfg.ncp_max,
            ncp_n=cfg.ncp_n,
            average_method=cfg.average_method,
        ).run()

    def _classify_structure(
        self,
        covered_tnjc2s: RecordTypedDf[CoveredTnJc2],
        genome: Genome,
        tn_loc: RecordTypedDf[RefTnLoc],
        iso_output: Path,
    ) -> RecordTypedDf[ClassifiedTnJc2]:
        """Step 8: Classify junction pair structures."""
        return ClassifyTnJc2StructureStep(
            covered_tnjc2s=covered_tnjc2s,
            genome=genome,
            tn_locs=tn_loc,
            output_dir=iso_output,
            transposition_threshold=self.config.min_amplicon_length,
        ).run()

    def _filter_candidates(
        self,
        classified_tnjc2s: RecordTypedDf[ClassifiedTnJc2],
        iso_output: Path,
    ) -> RecordTypedDf[FilteredTnJc2]:
        """Step 9: Filter candidates by amplicon length."""
        return FilterTnJc2CandidatesStep(
            classified_tnjc2s=classified_tnjc2s,
            output_dir=iso_output,
            min_amplicon_length=self.config.min_amplicon_length,
            max_amplicon_length=self.config.max_amplicon_length,
        ).run()

    def _create_synthetic_junctions(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        genome: Genome,
        tn_loc: RecordTypedDf[RefTnLoc],
        iso_output: Path,
    ) -> None:
        """Step 10: Create synthetic junction sequences."""
        return CreateSyntheticJunctionsStep(
            filtered_tnjc2s=filtered_tnjc2s,
            genome=genome,
            tn_locs=tn_loc,
            output_dir=iso_output,
            read_length=self._get_iso_read_length(),
        ).run()

    def _align_reads(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        iso_output: Path,
    ) -> None:
        """Step 11: Align reads to synthetic junctions."""
        cfg = self.config

        # Align isolate reads in isolate folder
        AlignReadsToJunctionsStep(
            filtered_tnjc2s=filtered_tnjc2s,
            output_dir=iso_output,
            iso_fastq_path=cfg.iso_path,
            anc_fastq_path=None,  # Only align isolate reads here
            threads=cfg.breseq_threads,
        ).run()

        if cfg.has_ancestor:
            AncAlignReadsToJunctionsStep(
                filtered_tnjc2s=filtered_tnjc2s,
                output_dir=cfg.anc_run_dir,
                iso_output_dir=iso_output,    # Source folder with junctions
                iso_fastq_path=cfg.anc_path,  # Ancestor reads aligned as "iso" in ancestor folder
                anc_fastq_path=None,
                threads=cfg.breseq_threads,
            ).run()

    def _analyze_alignments(
        self,
        filtered_tnjc2s: RecordTypedDf[FilteredTnJc2],
        iso_output: Path,
        anc_output: Optional[Path],
    ) -> RecordTypedDf[AnalyzedTnJc2]:
        """Step 12: Analyze read alignments."""
        return AnalyzeTnJc2AlignmentsStep(
            filtered_tnjc2s=filtered_tnjc2s,
            output_dir=iso_output,
            anc_output_dir=anc_output,  # Ancestor BAM files are read from here
            read_length=self._get_iso_read_length(),
            anc_read_length=self._get_anc_read_length() if self.config.has_ancestor else None,
            min_overlap=self.config.min_overlap,
            min_jct_cov=self.config.min_jct_cov,
            has_ancestor=self.config.has_ancestor,
        ).run()

    def _classify_candidates(
        self,
        analyzed_tnjc2s: RecordTypedDf[AnalyzedTnJc2],
        iso_output: Path,
    ) -> RecordTypedDf[AnalyzedTnJc2]:
        """Step 13: Final classification of candidates."""
        return ClassifyTnJc2CandidatesStep(
            analyzed_tnjc2s=analyzed_tnjc2s,
            output_dir=iso_output,
            has_ancestor=self.config.has_ancestor,
            min_jct_cov=self.config.min_jct_cov,
        ).run()

    def _export(
        self,
        analyzed_tnjc2s: RecordTypedDf[AnalyzedTnJc2],
        genome: Genome,
        iso_output: Path,
    ) -> None:
        """Step 14: Export results to CSV."""
        ExportTnJc2Step(
            analyzed_tnjc2s=analyzed_tnjc2s,
            genome=genome,
            output_dir=iso_output,
            ref_name=self.config.ref_name,
            iso_name=self.config.iso_name,
            anc_name=self.config.anc_name,
            copy_number_threshold=self.config.copy_number_threshold,
            del_copy_number_threshold=self.config.del_copy_number_threshold,
            filter_amplicon_length=self.config.filter_amplicon_length,
        ).run()


def run_pipeline(config: Config) -> RecordTypedDf[AnalyzedTnJc2]:
    """Run the AmpliFinder pipeline."""
    return Pipeline(config).run()

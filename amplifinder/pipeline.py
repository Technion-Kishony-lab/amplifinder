"""Pipeline orchestration for AmpliFinder."""
import shutil
import pandas as pd
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Tuple, Optional

from amplifinder.config import Config, get_iso_run_dir, get_anc_run_dir, load_config_from_run, save_config
from amplifinder.data_types import (
    Genome, RecordTypedDF, RefTnLoc, RefTnJunction, SeqRefTnSide, Junction, TnJunction, TnJc2,
    CoveredTnJc2, ClassifiedTnJc2, CandidateTnJc2, AnalyzedTnJc2,
)
from amplifinder.logger import info
from amplifinder.steps import (
    InitializingStep,
    GetRefGenomeStep,
    LocateTNsUsingISfinderStep,
    LocateTNsUsingGenbankStep,
    BreseqStep,
    CreateRefTnJcStep,
    CreateRefTnEndSeqsStep,
    CreateTnJcStep,
    CreateTnJc2Step,
    CalcAmpliconCoverageStep,
    ClassifyStructureStep,
    FilterCandidatesStep,
    CreateSyntheticJunctionsStep,
    AlignReadsToJunctionsStep,
    AnalyzeAlignmentsStep,
    ClassifyCandidatesStep,
    ExportStep,
)
from amplifinder.data import get_builtin_isfinder_db_path
from amplifinder.utils.tn_loc import compare_tn_locations


@dataclass
class Pipeline:
    """AmpliFinder pipeline with phased execution."""

    config: Config

    def run(self) -> RecordTypedDF[AnalyzedTnJc2]:
        """Run full pipeline, return analyzed candidates."""
        iso_output, anc_output = self._initialize()

        # Load reference genome (needed for ancestor breseq and isolate pipeline)
        genome = self._load_reference()
        
        # Handle ancestor breseq if needed
        self._ancestor_breseq(genome, anc_output)
        
        # Run isolate pipeline
        tn_loc = self._locate_tns_in_reference(genome)
        ref_tn_jc, ref_tn_end_seqs = self._create_reference_tn_junctions(tn_loc, genome, iso_output)
        breseq_jc = self._run_breseq(genome, iso_output)
        tnjc = self._create_tnjc(breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output)
        tnjc2 = self._create_tnjc2(tnjc, genome, iso_output)
        covered = self._calc_amplicon_coverage(tnjc2, genome, iso_output)
        classified = self._classify_structure(covered, tn_loc, iso_output)
        candidates = self._filter_candidates(classified, iso_output)
        self._create_synthetic_junctions(candidates, genome, tn_loc, iso_output)
        self._align_reads(candidates, iso_output)
        analyzed = self._analyze_alignments(candidates, iso_output, anc_output)
        classified_final = self._classify_candidates(analyzed, iso_output)
        self._export(classified_final, iso_output)
        
        return classified_final
    
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
        BreseqStep(
            output_path=anc_breseq_path,
            fastq_path=self.config.anc_path,
            ref_file=genome.genbank_path or genome.fasta_path,
            docker=self.config.breseq_docker,
            threads=self.config.breseq_threads,
        ).run()
        
        info("Ancestor breseq completed")
    
    def _copy_junctions_to_ancestor(
        self,
        candidates: RecordTypedDF[CandidateTnJc2],
        iso_output: Path,
    ) -> None:
        """Copy junction files from isolate to ancestor folder (only if not already there).
        
        This allows ancestor alignments to be stored in the ancestor folder
        and shared across multiple isolate runs.
        """
        if not self.config.has_ancestor:
            return
        
        anc_run_dir = get_anc_run_dir(self.config)
        
        for candidate in candidates:
            iso_jc_dir = iso_output / candidate.analysis_dir
            anc_jc_dir = anc_run_dir / candidate.analysis_dir
            anc_junctions = anc_jc_dir / "junctions.fasta"
            
            # Only copy if not already in ancestor folder
            if anc_junctions.exists():
                continue
            
            # Copy junctions.fasta from isolate to ancestor folder
            iso_junctions = iso_jc_dir / "junctions.fasta"
            if iso_junctions.exists():
                anc_jc_dir.mkdir(parents=True, exist_ok=True)
                shutil.copy2(iso_junctions, anc_junctions)
    

    @property
    def ref_tn_source(self) -> str:
        """Get the source of the reference TN data."""
        return self.config.use_isfinder and "isfinder" or "genbank"

    def _initialize(self) -> Tuple[Path, Optional[Path]]:
        """Step 0: Initialize output directories."""
        iso_output, anc_output = InitializingStep(config=self.config).run()
        return iso_output, anc_output

    def _load_reference(self) -> Genome:
        """Step 1: Get reference genome."""
        genome = GetRefGenomeStep(
            ref_name=self.config.ref_name,
            ref_path=self.config.ref_path,
            ncbi=self.config.ncbi,
        ).run()
        info(f"Reference: {genome.name} ({genome.length:,} bp)")
        return genome

    def _locate_tns_in_reference(self, genome: Genome) -> RecordTypedDF[RefTnLoc]:
        """Step 2: Locate TN elements using GenBank and/or ISfinder."""
        cfg = self.config

        tn_loc_dir = cfg.ref_path / "tn_loc" / genome.name

        # 2a: GenBank annotations
        tn_loc_genbank = LocateTNsUsingGenbankStep(
            genome=genome,
            output_dir=tn_loc_dir,
        ).run()

        if tn_loc_genbank is not None:
            info(f"GenBank: found {len(tn_loc_genbank)} TN elements")
        else:
            info("No GenBank file provided - skipping GenBank TN annotation")

        # 2b: ISfinder database
        tn_loc_isfinder = LocateTNsUsingISfinderStep(
            genome=genome,
            output_dir=tn_loc_dir,
            isdb_path=cfg.isdb_path or get_builtin_isfinder_db_path(),
            evalue=cfg.isfinder_evalue,
            critical_coverage=cfg.isfinder_critical_coverage,
        ).run()

        info(f"ISfinder: found {len(tn_loc_isfinder)} TN elements")

        # Compare sources
        if tn_loc_genbank is not None:
            compare_tn_locations(tn_loc_genbank, tn_loc_isfinder)

        # Select source
        if tn_loc_genbank is None and not cfg.use_isfinder:
            raise ValueError("No TN locations - provide GenBank or set --use-isfinder")
        tn_loc = tn_loc_isfinder if cfg.use_isfinder else tn_loc_genbank
        assert tn_loc is not None
        return tn_loc

    def _create_reference_tn_junctions(
        self, tn_loc: RecordTypedDF[RefTnLoc], genome: Genome, iso_output: Path
    ) -> Tuple[RecordTypedDF[RefTnJunction], RecordTypedDF[SeqRefTnSide]]:
        """Step 3: Create reference junctions and TN end sequences."""
        cfg = self.config

        ref_tn_jc = CreateRefTnJcStep(
            tn_loc=tn_loc,
            genome=genome,
            output_dir=iso_output,
            source=self.ref_tn_source,
            reference_tn_out_span=cfg.reference_IS_out_span,
        ).run()

        ref_tn_end_seqs = CreateRefTnEndSeqsStep(
            ref_tn_jc=ref_tn_jc,
            tn_loc=tn_loc,
            genome=genome,
            output_dir=iso_output,
            source=self.ref_tn_source,
            max_dist_to_tn=cfg.max_dist_to_IS,
        ).run()

        return ref_tn_jc, ref_tn_end_seqs

    def _run_breseq(self, genome: Genome, iso_output: Path) -> pd.DataFrame:
        """Step 4: Run breseq on isolate to get junctions."""
        cfg = self.config

        breseq_output = BreseqStep(
            fastq_path=cfg.iso_path,
            ref_file=genome.genbank_path or genome.fasta_path,
            output_path=cfg.iso_breseq_path or iso_output / "breseq",
            docker=cfg.breseq_docker,
            threads=cfg.breseq_threads,
        ).run()

        breseq_jc = breseq_output["JC"]
        info(f"breseq: {len(breseq_jc)} junctions")
        return breseq_jc

    def _create_tnjc(
        self,
        breseq_jc: pd.DataFrame,
        ref_tn_jc: RecordTypedDF[RefTnJunction],
        ref_tn_end_seqs: RecordTypedDF[SeqRefTnSide],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDF[TnJunction]:
        """Step 5: Match junctions to TN elements."""
        cfg = self.config

        # Combine breseq junctions with reference TN junctions
        all_jc_df = pd.concat([breseq_jc, ref_tn_jc.df], ignore_index=True)
        all_jc = RecordTypedDF(all_jc_df, Junction)

        tnjc = CreateTnJcStep(
            jc_df=all_jc,
            ref_tn_end_seqs=ref_tn_end_seqs,
            genome=genome,
            output_dir=iso_output,
            max_dist_to_tn=cfg.max_dist_to_IS,
            trim_jc_flanking=cfg.trim_jc_flanking,
        ).run()
        info(f"TnJc: {len(tnjc)} TN-associated junctions")
        return tnjc

    def _create_tnjc2(
        self,
        tnjc: RecordTypedDF[TnJunction],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDF[TnJc2]:
        """Step 6: Combine junction pairs (TnJc2)."""
        tnjc2 = CreateTnJc2Step(
            tnjc=tnjc,
            genome=genome,
            output_dir=iso_output,
        ).run()
        info(f"TnJc2: {len(tnjc2)} junction pairs")
        return tnjc2

    def _calc_amplicon_coverage(
        self,
        tnjc2: RecordTypedDF[TnJc2],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDF[CoveredTnJc2]:
        """Step 7: Calculate amplicon coverage."""
        cfg = self.config
        
        # Get ancestor breseq path if ancestor exists
        anc_breseq_path = None
        anc_name = None
        if cfg.has_ancestor:
            if cfg.anc_breseq_path:
                # User provided breseq path
                anc_breseq_path = cfg.anc_breseq_path
            else:
                # Use breseq from ancestor run directory
                anc_run_dir = get_anc_run_dir(cfg)
                anc_breseq_path = anc_run_dir / "breseq"
            anc_name = cfg.anc_name
        
        iso_breseq_path = cfg.iso_breseq_path or (iso_output / "breseq")
        
        covered = CalcAmpliconCoverageStep(
            tnjc2=tnjc2,
            genome=genome,
            iso_breseq_path=iso_breseq_path,
            output_dir=iso_output,
            ref_name=cfg.ref_name,
            iso_name=cfg.iso_name,
            anc_breseq_path=anc_breseq_path,
            anc_name=anc_name,
            min_amplicon_length=cfg.min_amplicon_length,
            max_amplicon_length=cfg.max_amplicon_length,
            ncp_limit1=cfg.ncp_limit1,
            ncp_limit2=cfg.ncp_limit2,
            ncp_n=cfg.ncp_n,
        ).run()
        info(f"Coverage: {len(covered)} candidates")
        return covered
    
    def _classify_structure(
        self,
        covered: RecordTypedDF[CoveredTnJc2],
        tn_loc: RecordTypedDF[RefTnLoc],
        iso_output: Path,
    ) -> RecordTypedDF[ClassifiedTnJc2]:
        """Step 8: Classify junction pair structures."""
        classified = ClassifyStructureStep(
            covered_tnjc2=covered,
            tn_locs=tn_loc,
            output_dir=iso_output,
            min_amplicon_length=self.config.min_amplicon_length,
        ).run()
        info(f"Classification: {len(classified)} candidates")
        return classified
    
    def _filter_candidates(
        self,
        classified: RecordTypedDF[ClassifiedTnJc2],
        iso_output: Path,
    ) -> RecordTypedDF[CandidateTnJc2]:
        """Step 9: Filter candidates by amplicon length."""
        candidates = FilterCandidatesStep(
            classified_tnjc2=classified,
            output_dir=iso_output,
            min_amplicon_length=self.config.min_amplicon_length,
            max_amplicon_length=self.config.max_amplicon_length,
        ).run()
        info(f"Filtered: {len(candidates)} candidates")
        return candidates
    
    def _create_synthetic_junctions(
        self,
        candidates: RecordTypedDF[CandidateTnJc2],
        genome: Genome,
        tn_loc: RecordTypedDF[RefTnLoc],
        iso_output: Path,
    ) -> None:
        """Step 10: Create synthetic junction sequences."""
        if len(candidates) == 0:
            return
        
        # Always create junction files for isolate candidates
        read_length = self.config.iso_read_length or 150
        
        CreateSyntheticJunctionsStep(
            candidates=candidates,
            genome=genome,
            tn_locs=tn_loc,
            output_dir=iso_output,
            read_length=read_length,
        ).run()
        info(f"Created synthetic junctions for {len(candidates)} candidates")
    
    def _align_reads(
        self,
        candidates: RecordTypedDF[CandidateTnJc2],
        iso_output: Path,
    ) -> None:
        """Step 11: Align reads to synthetic junctions."""
        if len(candidates) == 0:
            return
        
        cfg = self.config
        
        # Align isolate reads in isolate folder
        AlignReadsToJunctionsStep(
            candidates=candidates,
            output_dir=iso_output,
            iso_fastq_path=cfg.iso_path,
            anc_fastq_path=None,  # Only align isolate reads here
            threads=cfg.breseq_threads,
        ).run()
        info(f"Aligned isolate reads for {len(candidates)} candidates")
        
        # If ancestor exists, copy junctions to ancestor folder (if not already there),
        # then align ancestor reads in ancestor folder
        # This allows ancestor alignments to be shared across multiple isolate runs
        if cfg.has_ancestor:
            anc_run_dir = get_anc_run_dir(cfg)
            # Copy junctions first (only if not already there)
            self._copy_junctions_to_ancestor(candidates, iso_output)
            # Then align ancestor reads in ancestor folder
            AlignReadsToJunctionsStep(
                candidates=candidates,
                output_dir=anc_run_dir,
                iso_fastq_path=cfg.anc_path,  # Ancestor reads aligned as "iso" in ancestor folder
                anc_fastq_path=None,
                threads=cfg.breseq_threads,
            ).run()
            info(f"Aligned ancestor reads for {len(candidates)} candidates")
    
    def _analyze_alignments(
        self,
        candidates: RecordTypedDF[CandidateTnJc2],
        iso_output: Path,
        anc_output: Optional[Path],
    ) -> RecordTypedDF[AnalyzedTnJc2]:
        """Step 12: Analyze read alignments."""
        if len(candidates) == 0:
            return RecordTypedDF.empty(AnalyzedTnJc2)
        
        analyzed = AnalyzeAlignmentsStep(
            candidates=candidates,
            output_dir=iso_output,
            anc_output_dir=anc_output,  # Ancestor BAM files are read from here
            read_length=self.config.iso_read_length or 150,
            req_overlap=self.config.req_overlap,
            min_jct_cov=self.config.min_jct_cov,
            has_ancestor=self.config.has_ancestor,
        ).run()
        info(f"Analyzed: {len(analyzed)} candidates")
        return analyzed
    
    def _classify_candidates(
        self,
        analyzed: RecordTypedDF[AnalyzedTnJc2],
        iso_output: Path,
    ) -> RecordTypedDF[AnalyzedTnJc2]:
        """Step 13: Final classification of candidates."""
        if len(analyzed) == 0:
            return analyzed
        
        classified = ClassifyCandidatesStep(
            analyzed=analyzed,
            output_dir=iso_output,
            has_ancestor=self.config.has_ancestor,
            min_jct_cov=self.config.min_jct_cov,
        ).run()
        info(f"Final classification: {len(classified)} candidates")
        return classified
    
    def _export(
        self,
        analyzed: RecordTypedDF[AnalyzedTnJc2],
        iso_output: Path,
    ) -> None:
        """Step 14: Export results to CSV."""
        ExportStep(
            analyzed_candidates=analyzed,
            output_dir=iso_output,
            copy_number_threshold=self.config.copy_number_threshold,
            del_copy_number_threshold=self.config.del_copy_number_threshold,
            filter_amplicon_length=self.config.filter_amplicon_length,
        ).run()
        info("Exported results to CSV")


def run_pipeline(config: Config) -> RecordTypedDF[AnalyzedTnJc2]:
    """Run the AmpliFinder pipeline."""
    return Pipeline(config).run()

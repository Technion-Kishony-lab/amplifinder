"""Pipeline orchestration for AmpliFinder."""

import pandas as pd
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

from amplifinder.config import Config
from amplifinder.data_types import (
    Genome, RecordTypedDF, TnLoc, RefTnJunction, TnEndSeq, Junction, TnJunction, TnJunctionPair,
)
from amplifinder.logger import info
from amplifinder.steps import (
    InitializingStep,
    GetReferenceStep,
    LocateTNsUsingISfinderStep,
    LocateTNsUsingGenbankStep,
    BreseqStep,
    CreateRefTnJcsStep,
    CreateRefTnEndSeqsStep,
    CreateTNJCStep,
    CreateTNJC2Step,
)
from amplifinder.data import get_builtin_isfinder_db_path
from amplifinder.utils.tn_loc import compare_tn_locations


@dataclass
class Pipeline:
    """AmpliFinder pipeline with phased execution."""

    config: Config

    def run(self) -> RecordTypedDF[TnJunctionPair]:
        """Run full pipeline, return TNJC2 results."""
        iso_output = self._initialize()
        genome = self._load_reference()
        tn_loc = self._locate_tns_in_reference(genome)
        ref_tn_jc, ref_tn_end_seqs = self._create_reference_tn_junctions(tn_loc, genome, iso_output)
        breseq_jc = self._run_breseq(genome, iso_output)
        tnjc = self._create_tnjc(breseq_jc, ref_tn_jc, ref_tn_end_seqs, genome, iso_output)
        tnjc2 = self._create_tnjc2(tnjc, genome, iso_output)
        return tnjc2

    def _initialize(self) -> Path:
        """Step 0: Initialize output directories."""
        return InitializingStep(
            output_dir=self.config.output_dir,
            iso_name=self.config.iso_name,
        ).run()

    def _load_reference(self) -> Genome:
        """Step 1: Get reference genome."""
        genome = GetReferenceStep(
            ref_name=self.config.ref_name,
            ref_path=self.config.ref_path,
            ncbi=self.config.ncbi,
        ).run()
        info(f"Reference: {genome.name} ({genome.length:,} bp)")
        return genome

    def _locate_tns_in_reference(self, genome: Genome) -> RecordTypedDF[TnLoc]:
        """Step 2: Locate TN elements using GenBank and/or ISfinder."""
        cfg = self.config

        # 2a: GenBank annotations
        tn_loc_genbank = LocateTNsUsingGenbankStep(
            genbank_path=genome.genbank_path,
            ref_name=genome.name,
            ref_path=cfg.ref_path,
        ).run()

        if tn_loc_genbank is not None:
            info(f"GenBank: found {len(tn_loc_genbank)} TN elements")
        else:
            info("No GenBank file provided - skipping GenBank TN annotation")

        # 2b: ISfinder database
        tn_loc_isfinder = LocateTNsUsingISfinderStep(
            ref_fasta=genome.fasta_path,
            ref_name=genome.name,
            ref_path=cfg.ref_path,
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
        self, tn_loc: RecordTypedDF[TnLoc], genome: Genome, iso_output: Path
    ) -> Tuple[RecordTypedDF[RefTnJunction], RecordTypedDF[TnEndSeq]]:
        """Step 3: Create reference junctions and TN end sequences."""
        cfg = self.config

        ref_tn_jc = CreateRefTnJcsStep(
            tn_loc=tn_loc,
            ref_name=genome.name,
            output_dir=iso_output,
            reference_tn_out_span=cfg.reference_IS_out_span,
        ).run()

        ref_tn_end_seqs = CreateRefTnEndSeqsStep(
            ref_tn_jc=ref_tn_jc,
            genome=genome,
            ref_path=cfg.ref_path,
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
        ref_tn_end_seqs: RecordTypedDF[TnEndSeq],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDF[TnJunction]:
        """Step 5: Match junctions to TN elements."""
        cfg = self.config

        # Combine breseq junctions with reference TN junctions
        all_jc_df = pd.concat([breseq_jc, ref_tn_jc.df], ignore_index=True)
        all_jc = RecordTypedDF(all_jc_df, Junction)

        tnjc = CreateTNJCStep(
            jc_df=all_jc,
            ref_tn_end_seqs=ref_tn_end_seqs,
            genome=genome,
            output_dir=iso_output,
            max_dist_to_tn=cfg.max_dist_to_IS,
            trim_jc_flanking=cfg.trim_jc_flanking,
        ).run()
        info(f"TNJC: {len(tnjc)} TN-associated junctions")
        return tnjc

    def _create_tnjc2(
        self,
        tnjc: RecordTypedDF[TnJunction],
        genome: Genome,
        iso_output: Path,
    ) -> RecordTypedDF[TnJunctionPair]:
        """Step 6: Combine junction pairs (TNJC2)."""
        tnjc2 = CreateTNJC2Step(
            tnjc=tnjc,
            genome=genome,
            output_dir=iso_output,
        ).run()
        info(f"TNJC2: {len(tnjc2)} junction pairs")
        return tnjc2

    # TODO: Step 7 - Classification + Export


def run_pipeline(config: Config) -> RecordTypedDF[TnJunctionPair]:
    """Run the AmpliFinder pipeline."""
    return Pipeline(config).run()

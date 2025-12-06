"""Pipeline orchestration for AmpliFinder."""

import pandas as pd

from amplifinder.config import Config
from amplifinder.logger import info
from amplifinder.steps import (
    InitializingStep,
    GetReferenceStep,
    LocateTNsUsingISfinder,
    LocateTNsUsingGenbank,
    BreseqStep,
    CreateReferenceJunctionsStep,
    CreateTNJCStep,
)
from amplifinder.data import get_builtin_isfinder_db_path
from amplifinder.utils.tn_loc import compare_tn_locations


def run_pipeline(config: Config) -> None:
    """Run the AmpliFinder pipeline."""
    # Step 0: Initialize output directories
    iso_output = InitializingStep(
        output_dir=config.output_dir,
        iso_name=config.iso_name,
    ).run_and_read_outputs()

    # Step 1: Get reference genome
    genome = GetReferenceStep(
        ref_name=config.ref_name,
        ref_path=config.ref_path,
        ncbi=config.ncbi,
    ).run_and_read_outputs()
    info(f"Reference: {genome.name} ({genome.length:,} bp)")

    """Locate TN elements using GenBank and/or ISfinder databases."""

    # Step 2a: Locate TN elements using GenBank annotations
    tn_loc_genbank = LocateTNsUsingGenbank(
        genbank_path=genome.genbank_path,
        ref_name=genome.name,
        ref_path=config.ref_path,
    ).run_and_read_outputs()
    if tn_loc_genbank is not None:
        info(f"GenBank: found {len(tn_loc_genbank)} TN elements")
    else:
        info("No GenBank file provided - skipping GenBank TN annotation")

    # Step 2b: Locate TN elements using ISfinder database
    tn_loc_isfinder = LocateTNsUsingISfinder(
        ref_fasta=genome.fasta_path,
        ref_name=genome.name,
        ref_path=config.ref_path,
        isdb_path=config.isdb_path or get_builtin_isfinder_db_path(),
    ).run_and_read_outputs()
    info(f"ISfinder: found {len(tn_loc_isfinder)} TN elements")

    # Compare TN locations from both sources
    if tn_loc_genbank is not None:
        compare_tn_locations(tn_loc_genbank, tn_loc_isfinder)

    # Select which tn_loc to use based on config (fallback to isfinder if no genbank)
    if tn_loc_genbank is None and not config.use_isfinder:
        raise ValueError("No TN locations found - please provide a GenBank file or set --use-isfinder")
    tn_loc = tn_loc_isfinder if config.use_isfinder else tn_loc_genbank
    assert tn_loc is not None

    # Step 3: Create reference TN junctions
    ref_tn_jc = CreateReferenceJunctionsStep(
        tn_loc=tn_loc,
        ref_name=genome.name,
        output_dir=iso_output,
    ).run_and_read_outputs()

    """Run breseq on isolate to get junctions."""

    # Step 4: Run breseq on isolate to get junctions
    breseq_output = BreseqStep(
        fastq_path=config.iso_path,
        ref_file=genome.genbank_path or genome.fasta_path,
        output_path=config.iso_breseq_path or iso_output / "breseq",
        docker=config.breseq_docker,
    ).run_and_read_outputs()
    breseq_jc = breseq_output["JC"]
    info(f"breseq: {len(breseq_jc)} junctions")

    # Combine all junctions
    all_jc = pd.concat([breseq_jc, ref_tn_jc], ignore_index=True)

    # Step 5: Match junctions to TN elements → TNJC
    tnjc = CreateTNJCStep(
        jc_df=all_jc,
        tn_loc=tn_loc,
        ref_fasta=genome.fasta_path,
        output_dir=iso_output,
    ).run_and_read_outputs()
    info(f"TNJC: {len(tnjc)} TN-associated junctions")

    # TODO: Step 6 - Combine junction pairs (TNJC2)
    # TODO: Step 7 - Classification + Export

"""Pipeline orchestration for AmpliFinder."""

from amplifinder.config import Config
from amplifinder.logger import info
from amplifinder.steps import (
    InitializingStep,
    GetReferenceStep,
    LocateISsUsingISfinder,
    LocateISsUsingGenbank,
    BreseqStep,
)
from amplifinder.data import get_builtin_isfinder_db_path


def run_pipeline(config: Config) -> None:
    """Run the AmpliFinder pipeline."""
    # Step 0: Initialize output directories
    initializing_step = InitializingStep(output_dir=config.output_dir, iso_name=config.iso_name)
    iso_output = initializing_step.run_and_read_outputs()
    
    # Step 1: Get reference genome
    get_reference_step = GetReferenceStep(
        ref_name=config.ref_name,
        ref_path=config.ref_path,
        ncbi=config.ncbi,
    )
    genome = get_reference_step.run_and_read_outputs()
    info(f"Reference: {genome.name} ({genome.length:,} bp)")
    
    # Step 2a: Locate IS elements using GenBank annotations
    locIS_genbank_step = LocateISsUsingGenbank(
        genbank_path=genome.genbank_path,
        ref_name=genome.name,
        ref_path=config.ref_path,
    )
    IS_loc_genbank = locIS_genbank_step.run_and_read_outputs()
    info(f"GenBank: found {len(IS_loc_genbank)} IS elements")
    
    # Step 2b: Locate IS elements using ISfinder database
    isdb_path = config.isdb_path or get_builtin_isfinder_db_path()
    locIS_isfinder_step = LocateISsUsingISfinder(
        ref_fasta=genome.fasta_path,
        ref_name=genome.name,
        ref_path=config.ref_path,
        isdb_path=isdb_path,
    )
    IS_loc_isfinder = locIS_isfinder_step.run_and_read_outputs()
    info(f"ISfinder: found {len(IS_loc_isfinder)} IS elements")
    
    # Select which IS_loc to use based on config
    IS_loc = IS_loc_isfinder if config.use_ISfinder else IS_loc_genbank
    
    # Step 3: Run breseq on isolate
    iso_breseq = config.iso_breseq_path or iso_output / "breseq"
    breseq_step = BreseqStep(
        fastq_path=config.iso_path,
        ref_file=genome.genbank_path or genome.fasta_path,
        output_path=iso_breseq,
        docker=config.breseq_docker,
    )
    breseq_step.run()
    
    # TODO: Step 3 - IS detection (uses IS_loc)
    # TODO: Step 4 - Junction pairs
    # TODO: Step 5 - Classification + Export

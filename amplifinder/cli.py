"""Command-line interface for AmpliFinder."""

import logging
from pathlib import Path
from typing import Optional

import click

from amplifinder import __version__
from amplifinder.config import Config, load_config, merge_config
from amplifinder.logger import setup_logger, info, warning, error
from amplifinder.pipeline import Pipeline


@click.command()
@click.option(
    "-i", "--iso-path",
    type=click.Path(exists=True, path_type=Path),
    required=False,
    default=None,
    help="Path to isolate FASTQ file(s) or directory.",
)
@click.option(
    "-r", "--ref-name",
    type=str,
    required=True,
    help="Reference genome name (e.g., U00096 for E. coli K-12).",
)
@click.option(
    "-a", "--anc-path",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to ancestor FASTQ file(s) or directory.",
)
@click.option(
    "--iso-name",
    type=str,
    default=None,
    help="Isolate name (default: derived from path).",
)
@click.option(
    "--anc-name",
    type=str,
    default=None,
    help="Ancestor name (default: derived from path).",
)
@click.option(
    "--ref-path",
    type=click.Path(path_type=Path),
    default=Path("genomesDB"),
    help="Path to reference genome files (default: genomesDB).",
)
@click.option(
    "-o", "--output-dir",
    type=click.Path(path_type=Path),
    default=Path("output"),
    help="Output directory (default: output).",
)
@click.option(
    "--iso-breseq-path",
    type=click.Path(path_type=Path),
    default=None,
    help="Path to existing breseq output for isolate.",
)
@click.option(
    "--anc-breseq-path",
    type=click.Path(path_type=Path),
    default=None,
    help="Path to existing breseq output for ancestor.",
)
@click.option(
    "--ncbi/--no-ncbi",
    default=True,
    help="Fetch reference from NCBI (default: True).",
)
@click.option(
    "--use-isfinder/--no-use-isfinder",
    "use_isfinder",
    default=False,
    help="Use ISfinder database for IS detection (default: False).",
)
@click.option(
    "--config",
    "config_file",
    type=click.Path(exists=True, path_type=Path),
    default=None,
    help="Path to YAML/JSON config file.",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
    help="Logging level (default: INFO).",
)
@click.option(
    "--fetch-only",
    is_flag=True,
    default=False,
    help="Only fetch reference genome and create BLAST DB, then exit.",
)
@click.option(
    "--breseq-only",
    is_flag=True,
    default=False,
    help="Only run through breseq step, then exit.",
)
@click.option(
    "--visualize",
    type=click.Path(path_type=Path),
    default=None,
    help="Visualize results from a completed run directory.",
)
@click.option(
    "--save-plots",
    is_flag=True,
    default=False,
    help="Save plots to PNG files instead of displaying (use with --visualize).",
)
@click.version_option(version=__version__)
def main(
    iso_path: Path,
    ref_name: str,
    anc_path: Optional[Path],
    iso_name: Optional[str],
    anc_name: Optional[str],
    ref_path: Path,
    output_dir: Path,
    iso_breseq_path: Optional[Path],
    anc_breseq_path: Optional[Path],
    ncbi: bool,
    use_isfinder: bool,
    config_file: Optional[Path],
    log_level: str,
    fetch_only: bool,
    breseq_only: bool,
    visualize: Optional[Path],
    save_plots: bool,
) -> None:
    """AmpliFinder: Detect IS-mediated gene amplifications from WGS data."""
    # Setup logger early (use output_dir directly for fetch-only)
    log_dir = output_dir / ref_name
    log_dir.mkdir(parents=True, exist_ok=True)
    setup_logger(log_path=log_dir / "amplifinder.log", level=getattr(logging, log_level.upper()))

    info(f"AmpliFinder v{__version__}")

    try:
        # Visualization mode
        if visualize is not None:
            info(f"Visualizing results from: {visualize}")
            from amplifinder.visualization import visualize_candidates
            visualize_candidates(
                run_dir=visualize,
                save_plots=save_plots,
                interactive=not save_plots,
            )
            info("Done")
            return
        
        info(f"Reference: {ref_name}")

        # Fetch-only mode: only need ref_name and ref_path
        if fetch_only:
            info("Running fetch-only mode (reference + TN location)")
            from amplifinder.pipeline import run_fetch_only
            run_fetch_only(
                ref_name=ref_name,
                ref_path=ref_path,
                ncbi=ncbi,
                use_isfinder=use_isfinder,
            )
            info("Done")
            return

        # Full pipeline modes require iso_path
        if iso_path is None:
            raise click.ClickException("--iso-path is required (unless using --fetch-only)")

        # Load config file if provided
        file_config = None
        if config_file is not None:
            file_config = load_config(config_file)

        # Collect CLI arguments
        cli_args = {
            "iso_path": iso_path,
            "ref_name": ref_name,
            "anc_path": anc_path,
            "iso_name": iso_name,
            "anc_name": anc_name,
            "ref_path": ref_path,
            "output_dir": output_dir,
            "iso_breseq_path": iso_breseq_path,
            "anc_breseq_path": anc_breseq_path,
            "ncbi": ncbi,
            "use_isfinder": use_isfinder,
        }

        # Merge configurations
        merged = merge_config(cli_args, file_config)

        try:
            config = Config(**merged)
        except ValueError as e:
            raise click.ClickException(str(e))

        info(f"Isolate: {config.iso_path}")
        if config.anc_path:
            info(f"Ancestor: {config.anc_path}")
        else:
            warning("No ancestor assigned; using raw (non-normalized) coverage analysis")

        pipeline = Pipeline(config)
        if breseq_only:
            info("Running breseq-only mode")
            pipeline.run_breseq_only()
        else:
            pipeline.run()

    except Exception as e:
        error(f"Pipeline failed: {e}")
        raise click.ClickException(str(e))

    info("Done")


if __name__ == "__main__":
    main()

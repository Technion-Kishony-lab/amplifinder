"""Command-line interface for AmpliFinder."""

from pathlib import Path
from typing import Optional

import click

from amplifinder import __version__
from amplifinder.config import Config, load_config, merge_config
from amplifinder.pipeline import Pipeline


@click.command()
@click.option(
    "-i", "--iso-path",
    type=click.Path(path_type=Path),
    required=False,
    default=None,
    help="Path to isolate FASTQ file(s) or directory.",
)
@click.option(
    "-r", "--ref-name",
    type=str,
    required=False,
    default=None,
    help="Reference genome name (e.g., U00096 for E. coli K-12).",
)
@click.option(
    "-a", "--anc-path",
    type=click.Path(path_type=Path),
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
    type=click.Path(path_type=Path),
    default=None,
    help="Path to YAML/JSON config file.",
)
@click.option(
    "--create-config",
    "create_config_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Create a config file with current settings and exit (does not run pipeline).",
)
@click.option(
    "--log-level",
    type=click.Choice(["DEBUG", "INFO", "WARNING", "ERROR"], case_sensitive=False),
    default="INFO",
    help="Logging level (default: INFO).",
)
@click.option(
    "--breseq-only",
    is_flag=True,
    default=False,
    help="Only run through breseq step, then exit.",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Report which step is running and its output files.",
)
@click.option(
    "--create-plots/--no-create-plots",
    default=True,
    help="Create junction and amplicon coverage plots (default: True).",
)
@click.version_option(version=__version__)
def main(
    iso_path: Optional[Path],
    ref_name: Optional[str],
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
    create_config_path: Optional[Path],
    log_level: str,
    breseq_only: bool,
    verbose: bool,
    create_plots: bool,
) -> None:
    """AmpliFinder: Detect IS-mediated gene amplifications from WGS data."""
    try:
        # Load config file if provided
        file_config = None
        if config_file is not None:
            if not config_file.exists():
                raise click.ClickException(f"Config file not found: {config_file}")
            file_config = load_config(config_file)
        
        # Check required parameters (skip if only creating config)
        if create_config_path is None:
            # Validate required parameters
            if iso_path is None and (file_config is None or file_config.get("iso_path") is None):
                raise click.ClickException("--iso-path is required (via CLI or --config file)")
            
            if ref_name is None and (file_config is None or file_config.get("ref_name") is None):
                raise click.ClickException("--ref-name is required (via CLI or --config file)")
            
            # Validate paths exist for actual run
            if iso_path is not None and not iso_path.exists():
                raise click.ClickException(f"Isolate path does not exist: {iso_path}")
            
            if anc_path is not None and not anc_path.exists():
                raise click.ClickException(f"Ancestor path does not exist: {anc_path}")

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
            "create_plots": create_plots,
        }

        # Merge configurations
        merged = merge_config(cli_args, file_config)
        config = Config(**merged)
        
        # If --create-config, save and exit
        if create_config_path is not None:
            config.save_to_file(create_config_path, log=False)
            
            click.echo(f"Config file created: {create_config_path}")
            click.echo("\nRun with: amplifinder --config " + str(create_config_path))
            return
        
        # CLI startup message
        click.echo(f"AmpliFinder v{__version__}")

        # Run pipeline (pipeline handles all logger setup)
        pipeline = Pipeline(config, verbose=verbose)
        if breseq_only:
            click.echo("Running breseq-only mode")
            pipeline.run_breseq_only()
        else:
            pipeline.run()

    except Exception as e:
        raise click.ClickException(str(e))

    click.echo("Done")


if __name__ == "__main__":
    main()

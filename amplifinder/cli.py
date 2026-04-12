"""Command-line interface for AmpliFinder."""

from __future__ import annotations

import asyncio
import csv
from concurrent.futures import Executor, ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import click

from amplifinder import __version__
from amplifinder.config import Config, ISDetectionMethod, load_config, merge_config
from amplifinder.pipeline import Pipeline
from amplifinder.utils.dataclass_utils import convert_csv_row_types
from amplifinder.utils.file_lock import cleanup_lock_files
from amplifinder.utils.file_utils import fmt_separator


MAX_PARALLEL_DEFAULT = 2


def _lock_dirs_from_config(config: Config) -> List[Path]:
    """Collect all directories that may contain lock files for a given config."""
    dirs = [
        config.iso_run_dir,
        config.anc_run_dir if config.has_ancestor else None,
        config.ref_path,
        config.iso_breseq_path,
        config.anc_breseq_path,
    ]
    if config.is_detection_method == ISDetectionMethod.ISFINDER:
        from amplifinder.data import get_builtin_isfinder_db_path
        dirs.append(get_builtin_isfinder_db_path().parent)
    # Note: ISEScan doesn't need a separate lock dir - it writes to tn_loc_dir (under ref_path)
    return [d for d in dirs if d is not None]


@dataclass
class RunResult:
    row_idx: int
    run_id: str
    config: Config
    start_time: datetime
    end_time: datetime
    error: Optional[str] = None

    @property
    def exit_code(self) -> int:
        return 0 if self.error is None else 1

    def duration_seconds(self) -> int:
        """Calculate duration in seconds."""
        return int((self.end_time - self.start_time).total_seconds())


def _execute_pipeline(config: Config, breseq_only: bool, verbose: bool) -> Optional[str]:
    """Execute the pipeline and return error message or None on success."""
    try:
        pipeline = Pipeline(config, verbose=verbose)
        if breseq_only:
            pipeline.run_breseq_only()
        else:
            pipeline.run()
        return None
    except Exception as e:
        return str(e)


async def _run_one_batch(
    run_idx: int,
    run_id: str,
    config: Config,
    breseq_only: bool,
    verbose: bool,
    semaphore: asyncio.Semaphore,
    executor: Optional[Executor],
) -> RunResult:
    """Run one amplifinder job in batch mode."""
    async with semaphore:
        start_time = datetime.now()
        click.echo(f"[{run_id}] Starting at {start_time.strftime('%H:%M:%S')}")

        # Run in executor (thread pool or process pool)
        loop = asyncio.get_event_loop()
        error = await loop.run_in_executor(
            executor, _execute_pipeline, config, breseq_only, verbose
        )

        end_time = datetime.now()

        return RunResult(
            row_idx=run_idx,
            run_id=run_id,
            config=config,
            start_time=start_time,
            end_time=end_time,
            error=error,
        )


@click.command()
@click.option(
    "--batch-input",
    "batch_csv",
    type=click.Path(path_type=Path),
    default=None,
    help="Batch mode: input CSV with runs (columns = Config fields).",
)
@click.option(
    "--batch-output",
    type=click.Path(path_type=Path),
    default=None,
    help="Batch mode: output status CSV (default: run_status.csv).",
)
@click.option(
    "--max-parallel",
    type=click.IntRange(min=1),
    default=None,
    help="Batch mode: max parallel runs (default: 2).",
)
@click.option(
    "--use-processes",
    is_flag=True,
    default=False,
    help="Batch mode: use process pool instead of thread pool for true parallelism.",
)
@click.option(
    "-i", "--iso-path", "--iso-fastq-path", "iso_fastq_path",
    type=click.Path(path_type=Path),
    required=False,
    default=None,
    help="Directory containing isolate FASTQ files (*.fastq*). All matching files are used.",
)
@click.option(
    "-r", "--ref-name",
    type=str,
    required=False,
    default=None,
    help="Reference genome name (e.g., U00096 for E. coli K-12).",
)
@click.option(
    "-a", "--anc-path", "--anc-fastq-path", "anc_fastq_path",
    type=click.Path(path_type=Path),
    default=None,
    help="Directory containing ancestor FASTQ files (*.fastq*). All matching files are used.",
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
    default=None,
    help="Path to reference genome files (default: genomesDB).",
)
@click.option(
    "-o", "--output-dir",
    type=click.Path(path_type=Path),
    default=None,
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
    default=None,
    help="Fetch reference from NCBI (default: True).",
)
@click.option(
    "--is-detection-method",
    "is_detection_method",
    type=click.Choice(["genbank", "isfinder", "isescan"], case_sensitive=False),
    default=None,
    help="IS detection method to USE: genbank (default), isfinder, or isescan.",
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
    "--breseq-only",
    is_flag=True,
    default=None,
    help="Only run through breseq step, then exit.",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Report which step is running and its output files.",
)
@click.option(
    "--debug",
    is_flag=True,
    default=False,
    help="Enable debug-level logging.",
)
@click.option(
    "--create-plots/--no-create-plots",
    default=None,
    help="Create junction and amplicon coverage plots (default: True).",
)
@click.version_option(version=__version__)
def main(
    batch_csv: Optional[Path],
    batch_output: Optional[Path],
    max_parallel: Optional[int],
    use_processes: bool,
    iso_fastq_path: Optional[Path],
    ref_name: Optional[str],
    anc_fastq_path: Optional[Path],
    iso_name: Optional[str],
    anc_name: Optional[str],
    ref_path: Optional[Path],
    output_dir: Optional[Path],
    iso_breseq_path: Optional[Path],
    anc_breseq_path: Optional[Path],
    ncbi: Optional[bool],
    is_detection_method: Optional[str],
    config_file: Optional[Path],
    create_config_path: Optional[Path],
    breseq_only: Optional[bool],
    verbose: bool,
    debug: bool,
    create_plots: Optional[bool],
) -> None:
    """AmpliFinder: Detect IS-mediated gene amplifications from WGS data."""

    if debug:
        from amplifinder.env import DEBUG
        DEBUG.set(True)

    click.echo(f"AmpliFinder v{__version__}")
    if breseq_only:
        click.echo("Running breseq-only mode")

    try:
        # Load config file if provided
        file_config_kwargs = None
        if config_file is not None:
            if not config_file.exists():
                raise FileNotFoundError(f"Config file not found: {config_file}")
            file_config_kwargs = load_config(config_file)

        # Collect CLI arguments
        cli_config_kwargs = {
            "iso_fastq_path": iso_fastq_path,
            "ref_name": ref_name,
            "anc_fastq_path": anc_fastq_path,
            "iso_name": iso_name,
            "anc_name": anc_name,
            "ref_path": ref_path,
            "output_dir": output_dir,
            "iso_breseq_path": iso_breseq_path,
            "anc_breseq_path": anc_breseq_path,
            "ncbi": ncbi,
            "is_detection_method": is_detection_method,
            "create_plots": create_plots,
        }

        # Merge CLI args and config file once (priority: cli_args > file_config > defaults)
        config_kwargs = merge_config(file_config_kwargs, cli_config_kwargs)

        if batch_csv is None:
            # single mode
            if use_processes:
                raise ValueError("--use-processes requires --batch-csv")
            if max_parallel is not None:
                raise ValueError("--max-parallel requires --batch-csv")
            if batch_output is not None:
                raise ValueError("--batch-output requires --batch-csv")
            _run_single(
                config_kwargs=config_kwargs,
                create_config_path=create_config_path,
                breseq_only=breseq_only,
                verbose=verbose,
            )
        else:
            # batch mode
            if create_config_path is not None:
                raise ValueError("--create-config is not supported with --batch-csv")
            _run_batch(
                batch_csv=batch_csv,
                base_config_kwargs=config_kwargs,
                max_parallel=max_parallel,
                use_processes=use_processes,
                batch_output=batch_output,
                breseq_only=breseq_only,
                verbose=verbose,
            )
    except Exception as e:
        raise click.ClickException(str(e))


def _run_single(
    config_kwargs: Dict[str, Any],
    create_config_path: Optional[Path],
    breseq_only: bool,
    verbose: bool,
) -> None:
    # Create config from merged base (dataclass validates required fields)
    config = Config(**config_kwargs)

    # If --create-config, save and exit (skip validation)
    if create_config_path is not None:
        header = [
            "AmpliFinder run configuration template",
            f"Generated with: amplifinder --create-config {create_config_path.name}",
            "All fields shown with their default values.",
            "Override only the fields you need; the rest will use these defaults.",
        ]
        config.save_to_file(create_config_path, log=False, header=header)
        click.echo(f"Config file created: {create_config_path}")
        return

    # Validate required arguments before running
    arg_errors = config.validate_args()
    if arg_errors:
        joined = "\n  - ".join(arg_errors)
        raise ValueError(f"Missing required arguments:\n  - {joined}")

    # Validate input paths before running
    path_errors = config.validate_paths()
    if path_errors:
        joined = "\n  - ".join(path_errors)
        raise ValueError(f"Invalid paths:\n  - {joined}")

    # Run pipeline
    error = _execute_pipeline(config, breseq_only, verbose)
    if error is not None:
        raise RuntimeError(f"Pipeline failed:\n{error}")

    # Clean up lock files after run
    cleanup_lock_files(*_lock_dirs_from_config(config))

    click.echo("Done")


def _run_batch(
    batch_csv: Path,
    base_config_kwargs: Dict[str, Any],
    max_parallel: Optional[int],
    use_processes: bool,
    batch_output: Optional[Path],
    breseq_only: bool,
    verbose: bool,
) -> None:
    """Run AmpliFinder for multiple inputs from CSV."""
    if not batch_csv.exists():
        raise FileNotFoundError(f"CSV file not found: {batch_csv}")

    # Set default max_parallel
    max_parallel = max_parallel or MAX_PARALLEL_DEFAULT

    # Read CSV rows
    with batch_csv.open(newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)

    if not rows:
        raise ValueError("No runs found in CSV.")

    # Build configs for each row
    configs: List[Tuple[int, str, Config]] = []
    all_errors: List[str] = []

    for idx, row in enumerate(rows, start=1):
        try:
            # Extract and convert row values (Config.__post_init__ handles Path conversions)
            row_args = {k: v for k, v in row.items() if v and v.strip()}
            convert_csv_row_types(row_args, Config)

            # Merge: row_args > base_config (which is already cli_args > file_config > defaults)
            merged = merge_config(base_config_kwargs, row_args)

            config = Config(**merged)
            arg_errors = config.validate_args()
            path_errors = config.validate_paths()
            all_config_errors = arg_errors + path_errors
            if all_config_errors:
                all_errors.append(f"Row {idx}: {'; '.join(all_config_errors)}")
            else:
                run_id = row.get("run_id") or config.iso_name or f"row_{idx}"
                configs.append((idx, run_id, config))
        except Exception as exc:
            all_errors.append(f"Row {idx}: {exc}")

    if all_errors:
        details = "\n  - ".join(all_errors)
        raise ValueError(f"Invalid rows:\n  - {details}")

    # Print batch summary
    click.echo(fmt_separator("BATCH MODE STARTING"))
    click.echo(f"Number of runs: {len(configs)}")
    click.echo(f"Input CSV: {batch_csv}")
    click.echo(f"Max parallel runs: {max_parallel}")
    click.echo(f"Executor: {'ProcessPool' if use_processes else 'ThreadPool'}")
    click.echo(f"Verbose: {verbose}")
    if breseq_only:
        click.echo("Mode: breseq-only")
    click.echo("=" * 80)

    # Create executor (process pool or thread pool)
    executor: Optional[Executor] = None
    if use_processes:
        executor = ProcessPoolExecutor(max_workers=max_parallel)
    else:
        pass  # ThreadPool not explicitly created, using None for default

    # Run all jobs concurrently
    batch_start_time = datetime.now()

    async def run_all():
        semaphore = asyncio.Semaphore(max_parallel)
        tasks = []
        for run_idx, run_id, config in configs:
            tasks.append(
                asyncio.create_task(
                    _run_one_batch(run_idx, run_id, config, breseq_only, verbose, semaphore, executor)
                )
            )

        output_path = batch_output or Path("run_status.csv")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        results: List[RunResult] = []

        with output_path.open("w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(
                f,
                fieldnames=[
                    "run_id", "row_index", "iso_name", "anc_name", "ref_name",
                    "exit_code", "start_time", "end_time", "duration_sec", "output_dir", "error",
                ],
            )
            writer.writeheader()

            for task in asyncio.as_completed(tasks):
                result: RunResult = await task
                results.append(result)

                writer.writerow({
                    "run_id": result.run_id,
                    "row_index": result.row_idx,
                    "iso_name": result.config.iso_name or "",
                    "anc_name": result.config.anc_name or "",
                    "ref_name": result.config.ref_name,
                    "exit_code": result.exit_code,
                    "start_time": result.start_time.isoformat(),
                    "end_time": result.end_time.isoformat(),
                    "duration_sec": result.duration_seconds(),
                    "output_dir": str(result.config.output_dir),
                    "error": result.error or "",
                })
                f.flush()
                click.echo(f"[{result.run_id}] exit {result.exit_code} ({result.duration_seconds()}s)")

        return results

    try:
        results = asyncio.run(run_all())
        batch_end_time = datetime.now()
        total_elapsed = int((batch_end_time - batch_start_time).total_seconds())
        failures = [r for r in results if r.exit_code != 0]
        click.echo(fmt_separator("BATCH MODE COMPLETED"))
        click.echo(f"Completed {len(results)} run(s). Failures: {len(failures)}")
        click.echo(f"Total batch time: {total_elapsed}s ({total_elapsed // 60}m {total_elapsed % 60}s)")
        click.echo("=" * 80)
        raise SystemExit(1 if failures else 0)
    finally:
        # Clean up executor
        if executor is not None:
            executor.shutdown(wait=True)
        # Clean up lock files after all workers are done
        all_dirs: set[Path] = set()
        for cfg_tuple in configs:
            all_dirs.update(_lock_dirs_from_config(cfg_tuple[2]))
        cleanup_lock_files(*all_dirs)


if __name__ == "__main__":
    main()

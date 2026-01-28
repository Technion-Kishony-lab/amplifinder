#!/usr/bin/env python3
"""Run AmpliFinder for multiple inputs from a CSV file."""

from __future__ import annotations

import argparse
import asyncio
import csv
import shlex
import sys
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List, Optional


BOOL_TRUE = {"1", "true", "yes", "y", "t"}
BOOL_FALSE = {"0", "false", "no", "n", "f"}


@dataclass
class RunSpec:
    row_index: int
    run_id: str
    iso_path: Optional[str]
    ref_name: Optional[str]
    anc_path: Optional[str]
    iso_name: Optional[str]
    anc_name: Optional[str]
    ref_path: Optional[str]
    output_dir: Optional[str]
    iso_breseq_path: Optional[str]
    anc_breseq_path: Optional[str]
    ncbi: Optional[bool]
    use_isfinder: Optional[bool]
    create_plots: Optional[bool]
    breseq_only: Optional[bool]
    verbose: Optional[bool]
    log_level: Optional[str]
    config_path: Optional[str]
    extra_args: Optional[str]


@dataclass
class RunResult:
    spec: RunSpec
    command: List[str]
    exit_code: int
    status: str
    log_path: Path
    start_time: str
    end_time: str


def _parse_bool(value: Optional[str]) -> Optional[bool]:
    if value is None:
        return None
    raw = value.strip().lower()
    if raw == "":
        return None
    if raw in BOOL_TRUE:
        return True
    if raw in BOOL_FALSE:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def _get_field(row: Dict[str, str], name: str) -> Optional[str]:
    value = row.get(name)
    if value is None:
        return None
    value = value.strip()
    if value == "":
        return None
    return value


def _build_run_id(row: Dict[str, str], row_index: int) -> str:
    for key in ("run_id", "iso_name", "iso_path"):
        value = _get_field(row, key)
        if value:
            return value.replace("/", "_")
    return f"row_{row_index}"


def _read_csv(csv_path: Path) -> List[RunSpec]:
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)

    specs: List[RunSpec] = []
    for idx, row in enumerate(rows, start=1):
        spec = RunSpec(
            row_index=idx,
            run_id=_build_run_id(row, idx),
            iso_path=_get_field(row, "iso_path"),
            ref_name=_get_field(row, "ref_name"),
            anc_path=_get_field(row, "anc_path"),
            iso_name=_get_field(row, "iso_name"),
            anc_name=_get_field(row, "anc_name"),
            ref_path=_get_field(row, "ref_path"),
            output_dir=_get_field(row, "output_dir"),
            iso_breseq_path=_get_field(row, "iso_breseq_path"),
            anc_breseq_path=_get_field(row, "anc_breseq_path"),
            ncbi=_parse_bool(_get_field(row, "ncbi")),
            use_isfinder=_parse_bool(_get_field(row, "use_isfinder")),
            create_plots=_parse_bool(_get_field(row, "create_plots")),
            breseq_only=_parse_bool(_get_field(row, "breseq_only")),
            verbose=_parse_bool(_get_field(row, "verbose")),
            log_level=_get_field(row, "log_level"),
            config_path=_get_field(row, "config_path"),
            extra_args=_get_field(row, "extra_args"),
        )
        specs.append(spec)
    return specs


def _validate_spec(spec: RunSpec, validate_paths: bool) -> List[str]:
    errors: List[str] = []
    if spec.config_path is None:
        if not spec.iso_path:
            errors.append("iso_path is required (or provide config_path)")
        if not spec.ref_name:
            errors.append("ref_name is required (or provide config_path)")
    if validate_paths:
        for label, value in (
            ("iso_path", spec.iso_path),
            ("anc_path", spec.anc_path),
            ("ref_path", spec.ref_path),
            ("iso_breseq_path", spec.iso_breseq_path),
            ("anc_breseq_path", spec.anc_breseq_path),
            ("config_path", spec.config_path),
        ):
            if value and not Path(value).exists():
                errors.append(f"{label} does not exist: {value}")
    return errors


def _build_command(spec: RunSpec, python_exe: str) -> List[str]:
    cmd = [python_exe, "-m", "amplifinder"]
    if spec.iso_path:
        cmd += ["--iso-path", spec.iso_path]
    if spec.ref_name:
        cmd += ["--ref-name", spec.ref_name]
    if spec.anc_path:
        cmd += ["--anc-path", spec.anc_path]
    if spec.iso_name:
        cmd += ["--iso-name", spec.iso_name]
    if spec.anc_name:
        cmd += ["--anc-name", spec.anc_name]
    if spec.ref_path:
        cmd += ["--ref-path", spec.ref_path]
    if spec.output_dir:
        cmd += ["--output-dir", spec.output_dir]
    if spec.iso_breseq_path:
        cmd += ["--iso-breseq-path", spec.iso_breseq_path]
    if spec.anc_breseq_path:
        cmd += ["--anc-breseq-path", spec.anc_breseq_path]
    if spec.ncbi is True:
        cmd.append("--ncbi")
    if spec.ncbi is False:
        cmd.append("--no-ncbi")
    if spec.use_isfinder is True:
        cmd.append("--use-isfinder")
    if spec.use_isfinder is False:
        cmd.append("--no-use-isfinder")
    if spec.create_plots is True:
        cmd.append("--create-plots")
    if spec.create_plots is False:
        cmd.append("--no-create-plots")
    if spec.breseq_only is True:
        cmd.append("--breseq-only")
    if spec.verbose is True:
        cmd.append("--verbose")
    if spec.log_level:
        cmd += ["--log-level", spec.log_level]
    if spec.config_path:
        cmd += ["--config", spec.config_path]
    if spec.extra_args:
        cmd += shlex.split(spec.extra_args)
    return cmd


async def _run_one(
    spec: RunSpec,
    python_exe: str,
    amplifinder_cwd: Path,
    log_dir: Path,
    semaphore: asyncio.Semaphore,
) -> RunResult:
    async with semaphore:
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / f"{spec.run_id}.log"
        command = _build_command(spec, python_exe)
        start_time = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        with log_path.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"Run ID: {spec.run_id}\n")
            log_handle.write(f"Start: {start_time}\n")
            log_handle.write(f"Command: {' '.join(command)}\n\n")
            process = await asyncio.create_subprocess_exec(
                *command,
                cwd=str(amplifinder_cwd),
                stdout=log_handle,
                stderr=log_handle,
            )
            exit_code = await process.wait()
        end_time = datetime.utcnow().isoformat(timespec="seconds") + "Z"
        status = "success" if exit_code == 0 else "error"
        return RunResult(
            spec=spec,
            command=command,
            exit_code=exit_code,
            status=status,
            log_path=log_path,
            start_time=start_time,
            end_time=end_time,
        )


async def _run_all(
    specs: Iterable[RunSpec],
    python_exe: str,
    amplifinder_cwd: Path,
    log_dir: Path,
    max_parallel: int,
    status_csv: Path,
) -> List[RunResult]:
    semaphore = asyncio.Semaphore(max_parallel)
    tasks = [
        asyncio.create_task(_run_one(spec, python_exe, amplifinder_cwd, log_dir, semaphore))
        for spec in specs
    ]

    status_csv.parent.mkdir(parents=True, exist_ok=True)
    results: List[RunResult] = []
    with status_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "run_id",
                "row_index",
                "iso_name",
                "anc_name",
                "ref_name",
                "status",
                "exit_code",
                "start_time",
                "end_time",
                "log_path",
                "output_dir",
                "command",
            ],
        )
        writer.writeheader()
        for task in asyncio.as_completed(tasks):
            result = await task
            results.append(result)
            writer.writerow(
                {
                    "run_id": result.spec.run_id,
                    "row_index": result.spec.row_index,
                    "iso_name": result.spec.iso_name or "",
                    "anc_name": result.spec.anc_name or "",
                    "ref_name": result.spec.ref_name or "",
                    "status": result.status,
                    "exit_code": result.exit_code,
                    "start_time": result.start_time,
                    "end_time": result.end_time,
                    "log_path": str(result.log_path),
                    "output_dir": result.spec.output_dir or "",
                    "command": " ".join(result.command),
                }
            )
            handle.flush()
            print(
                f"[{result.spec.run_id}] {result.status} (exit {result.exit_code})",
                flush=True,
            )
    return results


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run AmpliFinder for multiple inputs defined in a CSV file."
    )
    parser.add_argument("--csv", required=True, type=Path, help="CSV file with runs.")
    parser.add_argument(
        "--max-parallel",
        type=int,
        default=2,
        help="Maximum number of runs in parallel (default: 2).",
    )
    parser.add_argument(
        "--amplifinder-cwd",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "amplifinder",
        help="Working directory for running amplifinder (default: repo/amplifinder).",
    )
    parser.add_argument(
        "--python",
        dest="python_exe",
        default=sys.executable,
        help="Python executable to use (default: current interpreter).",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "logs",
        help="Directory for wrapper logs (default: runner/logs).",
    )
    parser.add_argument(
        "--status-csv",
        type=Path,
        default=Path(__file__).resolve().parent / "run_status.csv",
        help="Status CSV output path (default: runner/run_status.csv).",
    )
    parser.add_argument(
        "--no-validate-paths",
        action="store_true",
        help="Skip pre-checking that input paths exist.",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    specs = _read_csv(args.csv)
    if not specs:
        print("No runs found in CSV.", file=sys.stderr)
        return 1

    all_errors: List[str] = []
    valid_specs: List[RunSpec] = []
    for spec in specs:
        errors = _validate_spec(spec, validate_paths=not args.no_validate_paths)
        if errors:
            joined = "; ".join(errors)
            all_errors.append(f"Row {spec.row_index} ({spec.run_id}): {joined}")
        else:
            valid_specs.append(spec)

    if all_errors:
        print("Invalid rows detected:", file=sys.stderr)
        for err in all_errors:
            print(f"  - {err}", file=sys.stderr)
        print("Fix the CSV or use --no-validate-paths to bypass checks.", file=sys.stderr)
        return 1

    results = asyncio.run(
        _run_all(
            valid_specs,
            python_exe=args.python_exe,
            amplifinder_cwd=args.amplifinder_cwd,
            log_dir=args.log_dir,
            max_parallel=max(1, args.max_parallel),
            status_csv=args.status_csv,
        )
    )

    failures = [result for result in results if result.exit_code != 0]
    print(f"Completed {len(results)} run(s). Failures: {len(failures)}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())

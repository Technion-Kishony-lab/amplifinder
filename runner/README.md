# AmpliFinder Runner

This wrapper runs multiple AmpliFinder jobs from a CSV, with a cap on the
number of concurrent runs.

## Usage

```bash
python runner/run_amplifinder_batch.py --csv runner/example_runs.csv --max-parallel 2
```

By default it runs `python -m amplifinder` from the `amplifinder/` directory,
writes wrapper logs to `runner/logs/`, and writes a run summary to
`runner/run_status.csv`.

## CSV format

Each row represents a single AmpliFinder run. These columns are supported:

- `run_id`: Optional. Used for log file names and reporting.
- `iso_path`: Required if `config_path` is not provided.
- `ref_name`: Required if `config_path` is not provided.
- `anc_path`: Optional.
- `iso_name`: Optional.
- `anc_name`: Optional.
- `ref_path`: Optional.
- `output_dir`: Optional.
- `iso_breseq_path`: Optional.
- `anc_breseq_path`: Optional.
- `ncbi`: Optional boolean (`true/false`, `1/0`, `yes/no`).
- `use_isfinder`: Optional boolean.
- `create_plots`: Optional boolean.
- `breseq_only`: Optional boolean.
- `verbose`: Optional boolean.
- `log_level`: Optional (`DEBUG`, `INFO`, `WARNING`, `ERROR`).
- `config_path`: Optional path to a YAML/JSON config file.
- `extra_args`: Optional extra CLI args passed as-is.

If `config_path` is provided, `iso_path` and `ref_name` are not required.

## Notes

- The wrapper validates that input paths exist by default.
- Disable path validation with `--no-validate-paths` if needed.

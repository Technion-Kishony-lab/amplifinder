"""Summarize MATLAB xlsx files into a single markdown file."""
import pandas as pd
from pathlib import Path


def summarize_matlab_xlsx(matlab_output_dir: Path, output_md: Path, n_rows: int = 3) -> None:
    """Create markdown summary of all xlsx files in matlab output directory.

    Args:
        matlab_output_dir: Path to MATLAB output directory (e.g., output/SRR25242877)
        output_md: Path to output markdown file
        n_rows: Number of rows to include (default 3)
    """
    xlsx_files = sorted(set(
        list(matlab_output_dir.glob('*.xlsx')) +
        list(matlab_output_dir.glob('**/*.xlsx'))
    ))

    lines = ["# MATLAB xlsx files summary\n", f"Source: `{matlab_output_dir}`\n"]

    for f in xlsx_files:
        rel = f.relative_to(matlab_output_dir)
        lines.append(f"\n## {rel}\n")

        try:
            df = pd.read_excel(f)
            lines.append(f"**Shape:** {df.shape[0]} rows × {df.shape[1]} columns\n")
            lines.append(f"\n**Columns:** `{list(df.columns)}`\n")
            lines.append(f"\n**First {min(n_rows, len(df))} rows:**\n")
            lines.append(df.head(n_rows).to_markdown(index=False))
            lines.append("\n")
        except Exception as e:
            lines.append(f"**Error:** {e}\n")

    output_md.write_text('\n'.join(lines))
    print(f"Written to {output_md}")


if __name__ == "__main__":
    matlab_dir = Path("/zdata/user-data/rkishony/AmpliFinder_test/AmpliFinderWorkspace/output/SRR25242877")
    output_file = Path("/zdata/user-data/rkishony/amplifinder/docs/matlab_xlsx_summary.md")
    output_file.parent.mkdir(exist_ok=True)
    summarize_matlab_xlsx(matlab_dir, output_file)

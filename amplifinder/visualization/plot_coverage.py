"""Plotting utilities for coverage visualization."""
import numpy as np
from amplifinder.optional_deps import plt
from matplotlib.axes import Axes
from pathlib import Path
from typing import Optional, Dict

from amplifinder.data_types.record_types import AnalyzedTnJc2


def _format_amplicon_axes(ax: Axes, left: int, right: int, left_pos: int, right_pos: int, scaf: str) -> None:
    """Apply common formatting to amplicon plot axes.

    Args:
        ax: Matplotlib axis to format
        left: Left boundary genomic position
        right: Right boundary genomic position
        left_pos: Left boundary position in plot coordinates
        right_pos: Right boundary position in plot coordinates
        scaf: Scaffold name for x-axis label
    """
    # Add amplicon boundaries (dashed red lines)
    ax.axvline(left_pos, color='red', linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axvline(right_pos, color='red', linestyle='--', linewidth=1.5, alpha=0.7)

    # Add tick labels at amplicon boundaries with actual genomic positions
    ax.set_xticks([left_pos, right_pos])
    ax.set_xticklabels([str(left), str(right)])

    ax.set_xlabel(f"Genomic position, '{scaf}'", fontsize=10)
    ax.grid(True, alpha=0.3)


def plot_amplicon_coverage(
    tnjc2: AnalyzedTnJc2,
    iso_scafs_to_covs: Dict[str, np.ndarray],
    anc_scafs_to_covs: Optional[Dict[str, np.ndarray]],
    output_path: Path,
    flank_fraction: float = 0.2
) -> None:
    """Plot isolate and ancestor coverage along an amplicon region.

    Args:
        iso_coverage: Dictionary mapping scaffold names to isolate coverage arrays
        anc_coverage: Dictionary mapping scaffold names to ancestor coverage arrays (optional)
        tnjc2: AnalyzedTnJc2 record containing amplicon information
        output_path: Path to save plot
    """
    scaf_obj = tnjc2.scaffold

    # Calculate plot range (scaffold-relative, 1-based)
    flank = int(tnjc2.amplicon_length * flank_fraction)
    plot_start = tnjc2.left - flank
    plot_end = tnjc2.right + flank
    if not scaf_obj.is_circular:
        plot_start = max(plot_start, 1)
        plot_end = min(plot_end, scaf_obj.length)

    # Get coverage for the amplicon region
    iso_scaf_cov = iso_scafs_to_covs[tnjc2.scaf]
    iso_plot_cov = scaf_obj.slice(plot_start, plot_end, seq=iso_scaf_cov)
    iso_plot_cov = iso_plot_cov / tnjc2.iso_scaf_avg

    anc_plot_cov = None
    if anc_scafs_to_covs is not None and tnjc2.anc_scaf_avg is not None:
        anc_scaf_cov = anc_scafs_to_covs[tnjc2.scaf]
        anc_plot_cov = scaf_obj.slice(plot_start, plot_end, seq=anc_scaf_cov)
        anc_plot_cov = anc_plot_cov / tnjc2.anc_scaf_avg

    positions = np.arange(0, len(iso_plot_cov))

    # Calculate amplicon boundaries in plot coordinates
    left_pos = tnjc2.left - plot_start
    right_pos = tnjc2.right - plot_start

    # Create plot with two subplots
    fig, axes = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    ax1: Axes = axes[0]
    ax2: Axes = axes[1]

    # Top subplot: Plot isolate and ancestor coverage
    ax1.axhline(1.0, color='black', linestyle='--', linewidth=0.5)
    ax1.plot(positions, iso_plot_cov, 'b-', linewidth=1, label='Isolate', alpha=0.7)
    ax1.hlines(
        tnjc2.iso_scaf_norm_copy_number, left_pos, right_pos, colors='blue', linestyles='--',
        linewidth=1.5, alpha=0.8
    )

    if anc_plot_cov is not None:
        ax1.plot(positions, anc_plot_cov, 'gray', linewidth=1, label='Ancestor', alpha=0.7)
        ax1.hlines(
            tnjc2.anc_scaf_norm_copy_number, left_pos, right_pos, colors='gray', linestyles='--',
            linewidth=1.5, alpha=0.8
        )

    title = f"{tnjc2.raw_event.value} | {tnjc2.left}-{tnjc2.right} | copy_number: {tnjc2.copy_number:.1f}x"
    ax1.set_title(title, fontsize=12, fontweight='bold')
    ax1.set_ylabel('Coverage, scaffold-normalized', fontsize=10)
    ax1.set_yscale('log')
    ax1.legend(loc='upper right', frameon=False)
    _format_amplicon_axes(ax1, tnjc2.left, tnjc2.right, left_pos, right_pos, tnjc2.scaf)

    # Bottom subplot: Plot ratio iso_cov / anc_cov
    if anc_plot_cov is not None:
        ax2.axhline(1.0, color='black', linestyle='--', linewidth=0.5)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = iso_plot_cov / anc_plot_cov
        ax2.plot(positions, ratio, 'g-', linewidth=1, label='Isolate/Ancestor', alpha=0.7)
        ax2.hlines(
            tnjc2.scaf_norm_copy_number_ratio, left_pos, right_pos, colors='green', linestyles='--',
            linewidth=1.5, alpha=0.8
        )
        ax2.legend(loc='upper right', frameon=False)
    else:
        ax2.text(0.5, 0.5, 'Ancestor coverage not available',
                 ha='center', va='center', transform=ax2.transAxes, fontsize=12)

    ax2.set_yscale('log')
    _format_amplicon_axes(ax2, tnjc2.left, tnjc2.right, left_pos, right_pos, tnjc2.scaf)
    ax2.set_ylabel('Isolate/Ancestor coverage ratio', fontsize=10)

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

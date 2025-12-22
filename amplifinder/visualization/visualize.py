"""Step 15: Visualize candidate amplicons with coverage plots."""

from pathlib import Path
from typing import Optional, List
import numpy as np

try:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from amplifinder.data_types import AnalyzedTnJc2, RecordTypedDF
from amplifinder.config import load_config_from_run
from amplifinder.tools.breseq import load_breseq_coverage
from amplifinder.utils.coverage import get_coverage_in_range
from amplifinder.logger import info, warning


def plot_candidate_coverage(
    candidate: AnalyzedTnJc2,
    iso_coverage: np.ndarray,
    anc_coverage: Optional[np.ndarray],
    genome_length: int,
    output_path: Optional[Path] = None,
    flank_fraction: float = 0.2,
    show: bool = True,
) -> None:
    """Plot coverage for a single candidate amplicon.
    
    Args:
        candidate: Analyzed candidate record
        iso_coverage: Full isolate coverage array
        anc_coverage: Full ancestor coverage array (optional)
        genome_length: Genome length for circular handling
        output_path: Path to save plot (optional)
        flank_fraction: Fraction of amplicon length to show as flanking region
        show: Whether to display plot interactively
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for visualization. Install with: pip install matplotlib")
    
    # Calculate plot range
    amplicon_length = candidate.amplicon_length
    flank = int(amplicon_length * flank_fraction)
    plot_start = max(1, candidate.pos_chr_L - flank)
    plot_end = min(genome_length, candidate.pos_chr_R + flank)
    
    # Handle origin spanning
    if candidate.span_origin:
        # For origin-spanning, we need to handle circular coordinates
        # Simplified: just plot the visible range
        pass
    
    # Get coverage in plot range
    iso_plot_cov = get_coverage_in_range(
        iso_coverage, plot_start, plot_end, candidate.span_origin, genome_length
    )
    positions = np.arange(plot_start, plot_start + len(iso_plot_cov))
    
    anc_plot_cov = None
    if anc_coverage is not None:
        anc_plot_cov = get_coverage_in_range(
            anc_coverage, plot_start, plot_end, candidate.span_origin, genome_length
        )
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot coverage lines
    ax.plot(positions, iso_plot_cov, 'b-', linewidth=1.5, label='Isolate', alpha=0.7)
    if anc_plot_cov is not None:
        ax.plot(positions, anc_plot_cov, 'gray', linewidth=1.5, label='Ancestor', alpha=0.7)
    
    # Add amplicon boundaries (dashed red lines)
    ax.axvline(candidate.pos_chr_L, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax.axvline(candidate.pos_chr_R, color='red', linestyle='--', linewidth=2, alpha=0.7)
    
    # Shade amplicon region
    ax.axvspan(candidate.pos_chr_L, candidate.pos_chr_R, alpha=0.2, color='red', label='Amplicon')
    
    # Formatting
    copy_num = candidate.copy_number if hasattr(candidate, 'copy_number') else 0.0
    title = (
        f"{candidate.raw_event.value} | "
        f"{candidate.pos_chr_L}-{candidate.pos_chr_R} | "
        f"copy_number: {copy_num:.1f}x"
    )
    ax.set_title(title, fontsize=12, fontweight='bold')
    ax.set_xlabel('Genomic position', fontsize=10)
    ax.set_ylabel('Coverage depth', fontsize=10)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save or show
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        info(f"Saved plot to {output_path}")
    
    if show:
        plt.show()
    else:
        plt.close()


def visualize_candidates(
    run_dir: Path,
    save_plots: bool = False,
    interactive: bool = True,
    flank_fraction: float = 0.2,
) -> None:
    """Visualize all candidates from a completed run.
    
    Args:
        run_dir: Path to run directory (contains run_config.yaml and results)
        save_plots: If True, save plots to PNG files instead of displaying
        interactive: If True and not save_plots, show interactive plots with navigation
        flank_fraction: Fraction of amplicon length to show as flanking region
    """
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib is required for visualization. Install with: pip install matplotlib")
    
    run_dir = Path(run_dir)
    
    # Load config
    try:
        config = load_config_from_run(run_dir)
    except FileNotFoundError:
        raise FileNotFoundError(f"Config file not found in {run_dir}. Run pipeline first.")
    
    # Load analyzed candidates
    analyzed_file = run_dir / "tn_jc2_classified.csv"
    if not analyzed_file.exists():
        analyzed_file = run_dir / "tn_jc2_analyzed.csv"
    
    if not analyzed_file.exists():
        raise FileNotFoundError(f"Analyzed candidates file not found in {run_dir}")
    
    from amplifinder.data_types import RecordTypedDF
    analyzed = RecordTypedDF.from_csv(analyzed_file, AnalyzedTnJc2)
    
    if len(analyzed) == 0:
        warning("No candidates to visualize")
        return
    
    # Load coverage data
    iso_breseq_path = config.iso_breseq_path or (run_dir / "breseq")
    iso_coverage = load_breseq_coverage(iso_breseq_path, config.ref_name)
    
    anc_coverage = None
    if config.has_ancestor:
        anc_breseq_path = config.anc_breseq_path or (run_dir / "breseq")
        try:
            anc_coverage = load_breseq_coverage(anc_breseq_path, config.ref_name)
        except FileNotFoundError:
            warning(f"Ancestor coverage not found at {anc_breseq_path}")
    
    # Get genome length (simplified - would need to load genome)
    genome_length = len(iso_coverage)
    
    # Visualize each candidate
    if save_plots:
        # Save all plots to files
        for candidate in analyzed:
            plot_path = run_dir / candidate.analysis_dir / "coverage_plot.png"
            plot_path.parent.mkdir(parents=True, exist_ok=True)
            plot_candidate_coverage(
                candidate=candidate,
                iso_coverage=iso_coverage,
                anc_coverage=anc_coverage,
                genome_length=genome_length,
                output_path=plot_path,
                flank_fraction=flank_fraction,
                show=False,
            )
        info(f"Saved {len(analyzed)} coverage plots")
    elif interactive:
        # Interactive mode with navigation
        _interactive_visualization(
            analyzed, iso_coverage, anc_coverage, genome_length, flank_fraction
        )
    else:
        # Show all plots sequentially
        for candidate in analyzed:
            plot_candidate_coverage(
                candidate=candidate,
                iso_coverage=iso_coverage,
                anc_coverage=anc_coverage,
                genome_length=genome_length,
                flank_fraction=flank_fraction,
                show=True,
            )


def _interactive_visualization(
    candidates: RecordTypedDF[AnalyzedTnJc2],
    iso_coverage: np.ndarray,
    anc_coverage: Optional[np.ndarray],
    genome_length: int,
    flank_fraction: float,
) -> None:
    """Interactive visualization with navigation buttons."""
    if not MATPLOTLIB_AVAILABLE:
        return
    
    from matplotlib.widgets import Button
    
    candidates_list = list(candidates)
    current_idx = [0]  # Use list to allow modification in nested function
    
    fig, ax = plt.subplots(figsize=(12, 6))
    plt.subplots_adjust(bottom=0.2)
    
    def update_plot():
        ax.clear()
        candidate = candidates_list[current_idx[0]]
        plot_candidate_coverage(
            candidate=candidate,
            iso_coverage=iso_coverage,
            anc_coverage=anc_coverage,
            genome_length=genome_length,
            flank_fraction=flank_fraction,
            show=False,
        )
        plt.draw()
    
    def next_candidate(event):
        if current_idx[0] < len(candidates_list) - 1:
            current_idx[0] += 1
            update_plot()
    
    def prev_candidate(event):
        if current_idx[0] > 0:
            current_idx[0] -= 1
            update_plot()
    
    # Add buttons
    ax_prev = plt.axes([0.3, 0.05, 0.15, 0.075])
    ax_next = plt.axes([0.55, 0.05, 0.15, 0.075])
    btn_prev = Button(ax_prev, 'Previous')
    btn_next = Button(ax_next, 'Next')
    
    btn_prev.on_clicked(prev_candidate)
    btn_next.on_clicked(next_candidate)
    
    # Show first candidate
    update_plot()
    plt.show()

"""Plotting utilities for read alignment visualization."""
import numpy as np

from amplifinder.optional_deps import plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from pathlib import Path
from typing import Dict, List, Tuple

from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.data_types.enums import Element, Side, ReadType, JcCall
from amplifinder.visualization.genetic_elements import draw_genetic_element, draw_amplicon_structure


READ_TYPES_TO_COLORS = {
    Side.LEFT: 'red',
    Side.MIDDLE: 'green',
    Side.RIGHT: 'red',
    None: 'lightgrey'
}


JC_CALLS_TO_COLORS = {
    True: 'green',
    False: 'red',
    None: 'grey'
}


def format_read_info(jc_cov: JunctionReadCounts, scales: JunctionReadCounts) -> str:
    """Format downsampling info text."""
    lines = []
    lines.append(f"L={jc_cov.left}, 1:{scales.left}")
    lines.append(f"R={jc_cov.right}, 1:{scales.right}")
    lines.append(f"S={jc_cov.spanning}, 1:{scales.spanning}")
    lines.append(f"U={jc_cov.undetermined}, 1:{scales.undetermined}")
    return "\n".join(lines)


def _down_sample_alignments(alignments: List[Tuple[int, int, ReadType]], jc_cov: JunctionReadCounts, max_reads_per_plot: int) -> tuple[List[Tuple[int, int, ReadType]], JunctionReadCounts]:
    """Downsample alignments to a maximum number of reads per plot.
    
    Left and right alignments are scaled independently to max_reads_per_plot.
    Spanning and undetermined alignments are scaled by the minimum of left/right scaling factors.
    
    Returns:
        Tuple of (downsampled_alignments, scales)
        where scales is a JunctionReadCounts with left/right/spanning/undetermined as scaling factors
    """
    # Calculate scaling factors
    max_reads_per_side = max_reads_per_plot // 2
    if jc_cov.left > jc_cov.right:
        left_scale = max(1, jc_cov.left // max_reads_per_side)
        right_scale = max(1, (jc_cov.right + jc_cov.spanning + jc_cov.undetermined) // max_reads_per_side)
        spanning_scale = right_scale
    else:
        right_scale = max(1, jc_cov.right // max_reads_per_side)
        left_scale = max(1, (jc_cov.left + jc_cov.spanning + jc_cov.undetermined) // max_reads_per_side)
        spanning_scale = left_scale
    
    scales = JunctionReadCounts(left=left_scale, right=right_scale, spanning=spanning_scale, undetermined=spanning_scale)
    
    counters = JunctionReadCounts()
    
    downsampled = []
    for start, end, read_type in alignments:
        count = counters.counts[read_type]
        scale = scales.counts[read_type]
        if count % scale == 0:
            downsampled.append((start, end, read_type))
        counters.increment(read_type)
    
    return downsampled, scales


def plot_junctions_coverage(
    jc_lengths: Dict[JunctionType, int],
    alignment_data: Dict[JunctionType, List[Tuple[int, int, ReadType]]],
    alignment_data_anc: Dict[JunctionType, List[Tuple[int, int, ReadType]]] | None = None,
    jc_covs: Dict[JunctionType, JunctionReadCounts] | None = None,
    jc_covs_anc: Dict[JunctionType, JunctionReadCounts] | None = None,
    jc_calls: Dict[JunctionType, JcCall] | None = None,
    jc_calls_anc: Dict[JunctionType, JcCall] | None = None,
    title: str | None = None,
    output_path: Path | str = 'junctions_coverage.png',
    max_reads_per_plot: int = 200,
) -> None:
    """Plot junctions coverage by reads as horizontal lines (coverage plot).
    
    Args:
        alignment_data: Dict mapping JunctionType to list of (start, end, ReadType) tuples for isolate
        jc_lengths: Dict mapping JunctionType to junction length
        title: Plot title
        output_path: Path to save the PNG file
        max_reads_per_plot: Maximum reads to show per plot (downsampling)
        alignment_data_anc: Optional dict for ancestor reads (plotted below x-axis with negative y)
        jc_calls: Optional dict mapping JunctionType to coverage call (True/False/None) for isolate
        jc_calls_anc: Optional dict mapping JunctionType to coverage call for ancestor
    """
    with_calls = jc_calls is not None
    if with_calls:
        if alignment_data_anc and not jc_calls_anc:
            raise ValueError("Cannot plot with junction calls for isolate and no ancestor calls")

    h_pixels = 18

    # Create figure with 7 subplots (each with genetic element axes above)
    fig = plt.figure(figsize=(15, 11))
    
    # Junction plots grid - adjust to leave room for top axes
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.1, wspace=0.3,
                           top=0.95, bottom=0.05, left=0.08, right=0.92)
    
    for jt in JunctionType:
        if jt.order == 0:
            row, col = 1, 0
        else:
            row, col = (jt.side.value + 1) // 2, abs(jt.order)
        
        # Create sub-gridspec for this position: [genetic_element_axes, alignment_axes]
        # Make genetic element axes ~10% of height
        inner_gs = gridspec.GridSpecFromSubplotSpec(2, 1, 
                                                    subplot_spec=gs[row - 1, col - 1],
                                                    height_ratios=[0.08, 1],
                                                    hspace=0.01)
        
        # Genetic element axes (narrow, on top)
        ax_gene = fig.add_subplot(inner_gs[0])
        ax_gene.axis('off')
        
        # Main alignment axes
        ax = fig.add_subplot(inner_gs[1])
        
        alignments = alignment_data[jt]
        jc_length = jc_lengths[jt]  # TODO: what happen when anc/iso have different lengths?
        jc_cov = jc_covs[jt]
        midpoint = jc_length // 2
        
        # Sort alignments by start position
        sorted_alignments = sorted(alignments, key=lambda x: x[0])
        
        # Downsample isolate
        downsampled, scales_iso = _down_sample_alignments(sorted_alignments, jc_cov, max_reads_per_plot)
        
        # Group by read type for efficient plotting (isolate - positive y)
        segments_by_type = {read_type: {'x': [], 'y': []} for read_type in READ_TYPES_TO_COLORS}
        
        # start, end are 0-based end-exclusive coordinates.
        #  |---|---|---|---|---|---|---|---|
        # -4  -3  -2  -1   0   1   2   3   4  x-axis
        #    0   1   2   3   4   5   6   7    seq index
        # 
        # junction_length = 8, midpoint = 4
        # read with start=3, end=7 will be plotted as (len = 4)
        #              |---|---|---|---|
        #             -1               3      x-axis
        # x_start = start - midpoint = 3 - 4 = -1
        # x_end = end - midpoint = 7 - 4 = 3

        for y_pos, (start, end, read_type) in enumerate(downsampled, start=1):
            x_start = start - midpoint
            x_end = end - midpoint
            # Collect as [start, end] pairs
            segments_by_type[read_type]['x'].append([x_start, x_end])
            segments_by_type[read_type]['y'].append([y_pos, y_pos])
        
        # Plot each read type using 2D numpy arrays (isolate - above x-axis)
        for read_type, segments in segments_by_type.items():
            if segments['x']:
                x_arr = np.array(segments['x']).T  # shape (2, n_segments)
                y_arr = np.array(segments['y']).T  # shape (2, n_segments)
                ax.plot(x_arr, y_arr, color=READ_TYPES_TO_COLORS[read_type], 
                        linewidth=1.5, solid_capstyle='butt')
        
        # Plot ancestor reads below x-axis (negative y) if available
        if alignment_data_anc:
            alignments_anc = alignment_data_anc[jt]
            jc_cov_anc = jc_covs_anc[jt]
            sorted_alignments_anc = sorted(alignments_anc, key=lambda x: x[0])
            
            downsampled_anc, scales_anc = _down_sample_alignments(sorted_alignments_anc, jc_cov_anc, max_reads_per_plot)
            
            segments_by_type_anc = {read_type: {'x': [], 'y': []} for read_type in READ_TYPES_TO_COLORS}
            
            for y_pos, (start, end, read_type) in enumerate(downsampled_anc, start=1):
                x_start = start - midpoint
                x_end = end - midpoint
                # Negative y position for ancestor
                segments_by_type_anc[read_type]['x'].append([x_start, x_end])
                segments_by_type_anc[read_type]['y'].append([-y_pos, -y_pos])
            
            # Plot ancestor reads with negative y values
            for read_type, segments in segments_by_type_anc.items():
                if segments['x']:
                    x_arr = np.array(segments['x']).T
                    y_arr = np.array(segments['y']).T
                    ax.plot(x_arr, y_arr, color=READ_TYPES_TO_COLORS[read_type], 
                            linewidth=1.5, solid_capstyle='butt', alpha=0.7)
        
        # Set consistent y-axis (extend to negative if ancestor data exists)
        y_min = -(max_reads_per_plot + 1) if alignment_data_anc else 0
        y_max = max_reads_per_plot + 1
        
        # Add vertical line at junction point with color based on jc_calls (isolate) - from y=0 to y_max
        jc_call_color = JC_CALLS_TO_COLORS[jc_calls[jt]] if with_calls else 'black'
        ax.plot([0, 0], [0, y_max], color=jc_call_color, linestyle='--', linewidth=2.5)
        
        # Add vertical line for jc_calls_anc (ancestor) if available - from y=0 to y_min
        if alignment_data_anc:
            jc_call_anc_color = JC_CALLS_TO_COLORS[jc_calls_anc[jt]] if with_calls else 'black'
            ax.plot([0, 0], [y_min, 0], color=jc_call_anc_color, linestyle='--', linewidth=2.5)
        
        # Add horizontal line at y=0 if ancestor data exists
        if alignment_data_anc:
            ax.axhline(0, color='black', linestyle='-', linewidth=0.5)
        ax.set_ylim(y_min, y_max)
        ax.set_ylabel('Read number (iso +, anc -)' if alignment_data_anc else 'Read number')
        
        # Set x-axis limits centered at junction
        ax.set_xlim(-midpoint, midpoint)
        ax.set_xlabel('Position relative to junction (bp)')
        ax.grid(True, alpha=0.3, axis='x')
        
        # Draw genetic elements in the narrow axes above
        left_elem_type, right_elem_type = jt.elements
        
        # Set up genetic element axes
        ax_gene.set_xlim(-midpoint, midpoint)
        ax_gene.set_ylim(0, 1)
        ax_gene.set_aspect('auto')
        
        # Add title on top of genetic element axes
        # ax_gene.set_title(f'{jt.name}

        if alignment_data_anc:
            iso_text = "iso:\n" + format_read_info(jc_cov, scales_iso)
            anc_text = "anc:\n" + format_read_info(jc_cov_anc, scales_anc)
            ax.text(0.05, 0.95, iso_text, fontsize=7, transform=ax.transAxes, ha='left', va='top', verticalalignment='top')
            ax.text(0.05, 0.05, anc_text, fontsize=7, transform=ax.transAxes, ha='left', va='bottom', verticalalignment='bottom')
        else:
            iso_text = format_read_info(jc_cov, scales_iso)
            ax.text(0.05, 0.95, iso_text, fontsize=7, transform=ax.transAxes, ha='left', va='top', verticalalignment='top')

        # Draw left element (ending at x=0)
        draw_genetic_element(ax_gene, 0.1, -midpoint, 0, left_elem_type, h_pixels=h_pixels, wave_tail=True)
        
        # Draw right element (starting at x=0)
        draw_genetic_element(ax_gene, 0.1, 0, midpoint, right_elem_type, h_pixels=h_pixels, wave_head=True)
    
    # Add legend in empty subplot position (2, 4)
    ax_legend = fig.add_subplot(gs[1, 3])
    ax_legend.axis('off')
    ax_legend.set_xlim(0, 1)
    ax_legend.set_ylim(0, 1)
    
    # Draw genetic element legend at top - each as two parts with waves
    y_start = 0.85
    element_height = 0.08
    left_start = 0
    left_end = 0.2
    right_start = 0.22
    right_end = 0.42
    label_x = 0.5
    
    # Chromosome
    draw_genetic_element(ax_legend, y_start, left_start, left_end, 
                        Element.CHR, h_pixels=h_pixels, wave_tail=True)
    draw_genetic_element(ax_legend, y_start, right_start, right_end, 
                        Element.CHR, h_pixels=h_pixels, wave_head=True)
    ax_legend.text(label_x, y_start, 'chromosome', fontsize=10, va='bottom')
    
    # Amplicon
    y_amp = y_start - element_height - 0.05
    draw_genetic_element(ax_legend, y_amp, left_start, left_end, 
                        Element.AMP, h_pixels=h_pixels, wave_head=True)
    draw_genetic_element(ax_legend, y_amp, right_start, right_end, 
                        Element.AMP, h_pixels=h_pixels, wave_tail=True)
    ax_legend.text(label_x, y_amp, 'amplicon', fontsize=10, va='bottom')
    
    # IS element
    y_is = y_amp - element_height - 0.05
    draw_genetic_element(ax_legend, y_is, left_start, left_end, 
                        Element.TN, h_pixels=h_pixels, wave_head=True)
    draw_genetic_element(ax_legend, y_is, right_start, right_end, 
                        Element.TN, h_pixels=h_pixels, wave_tail=True)
    ax_legend.text(label_x, y_is, 'IS', fontsize=10, va='bottom')
    
    # Read type legend at bottom (combine left/right since they share the same color)
    legend_elements = [
        Line2D([0], [0], color=READ_TYPES_TO_COLORS[Side.MIDDLE], linewidth=2, label='spanning'),
        Line2D([0], [0], color=READ_TYPES_TO_COLORS[Side.LEFT], linewidth=2, label='left/right'),
        Line2D([0], [0], color=READ_TYPES_TO_COLORS[None], linewidth=2, label='undetermined')
    ]
    
    # Add jc_calls color legend if available (applies to both jc_calls and jc_calls_anc)
    if with_calls:
        jc_legend_elements = [
            Line2D([0], [0], color='green', linestyle='--', linewidth=2.5, label='Positive'),
            Line2D([0], [0], color='red', linestyle='--', linewidth=2.5, label='Negative'),
            Line2D([0], [0], color='grey', linestyle='--', linewidth=2.5, label='Undetermined'),
        ]
        legend_elements.extend(jc_legend_elements)
    
    # Set legend title based on downsampling
    legend_title = None
    
    ax_legend.legend(handles=legend_elements, loc='lower center', fontsize=11, frameon=True, title=legend_title)
    
    # Add overall title
    fig.suptitle(title, fontsize=14, fontweight='bold')
    
    # Save figure
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    # Demo plot with no reads - just to see layout with genetic elements
    alignment_data = {jt: [] for jt in JunctionType.sorted()}
    jct_lengths = {jt: 600 for jt in JunctionType.sorted()}
    
    plot_junctions_coverage(
        alignment_data=alignment_data,
        jc_lengths=jct_lengths,
        title='Alignment Coverage Demo (No Reads)',
        output_path='alignment_coverage_demo.png'
    )

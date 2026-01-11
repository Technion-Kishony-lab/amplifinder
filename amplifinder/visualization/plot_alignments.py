"""Plotting utilities for read alignment visualization."""
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from pathlib import Path
from typing import Dict, List, Tuple

from amplifinder.data_types import JunctionType
from amplifinder.visualization.genetic_elements import draw_genetic_element, draw_amplicon_structure, GeneticElementType


READ_TYPES_TO_COLORS = {
    'spanning': 'green',
    'left': 'red',
    'right': 'red',
    'undetermined': 'lightgrey'
}


JUNCTION_TYPES_TO_PLOT_POSITIONS = {
    JunctionType.CHR_TO_AMP_LEFT: (1, 1),
    JunctionType.CHR_TO_TN_LEFT: (2, 1),
    JunctionType.AMP_RIGHT_TO_TN_LEFT: (1, 2),
    JunctionType.AMP_RIGHT_TO_AMP_LEFT: (1, 4),
    JunctionType.TN_RIGHT_TO_AMP_LEFT: (2, 2),
    JunctionType.TN_RIGHT_TO_CHR: (2, 3),
    JunctionType.AMP_RIGHT_TO_CHR: (1, 3),
}


JUNCTION_TYPES_TO_LEFT_RIGHT_ELEMENTS = {
    JunctionType.CHR_TO_AMP_LEFT: (GeneticElementType.CHROMOSOME, GeneticElementType.AMPLICON),
    JunctionType.CHR_TO_TN_LEFT: (GeneticElementType.CHROMOSOME, GeneticElementType.TN_ELEMENT),
    JunctionType.AMP_RIGHT_TO_TN_LEFT: (GeneticElementType.AMPLICON, GeneticElementType.TN_ELEMENT),
    JunctionType.AMP_RIGHT_TO_AMP_LEFT: (GeneticElementType.AMPLICON, GeneticElementType.AMPLICON),
    JunctionType.TN_RIGHT_TO_AMP_LEFT: (GeneticElementType.TN_ELEMENT, GeneticElementType.AMPLICON),
    JunctionType.TN_RIGHT_TO_CHR: (GeneticElementType.TN_ELEMENT, GeneticElementType.CHROMOSOME),
    JunctionType.AMP_RIGHT_TO_CHR: (GeneticElementType.AMPLICON, GeneticElementType.CHROMOSOME),
}


def plot_alignment_coverage(
    alignment_data: Dict[JunctionType, List[Tuple[int, int, str]]],
    jct_lengths: Dict[JunctionType, int],
    title: str,
    output_path: Path | str ,
    max_reads_per_plot: int = 200,
    alignment_data_anc: Dict[JunctionType, List[Tuple[int, int, str]]] | None = None,
) -> None:
    """Plot read alignments as horizontal lines (coverage plot).
    
    Args:
        alignment_data: Dict mapping JunctionType to list of (start, end, read_type) tuples for isolate
        jct_lengths: Dict mapping JunctionType to junction length
        title: Plot title
        output_path: Path to save the PNG file
        max_reads_per_plot: Maximum reads to show per plot (downsampling)
        alignment_data_anc: Optional dict for ancestor reads (plotted below x-axis with negative y)
    """

    h_pixels = 18

    # Find max number of reads across all junction types for isolate and ancestor
    max_reads_iso = max(len(alignments) for alignments in alignment_data.values())
    if max_reads_iso == 0:
        max_reads_iso = 1
    
    max_reads_anc = 0
    if alignment_data_anc:
        max_reads_anc = max(len(alignments) for alignments in alignment_data_anc.values())
        if max_reads_anc == 0:
            max_reads_anc = 1
    
    # Determine separate downsampling factors for isolate and ancestor
    downsample_factor_iso = max(1, max_reads_iso // max_reads_per_plot)
    downsample_factor_anc = max(1, max_reads_anc // max_reads_per_plot) if max_reads_anc > 0 else 1
    max_reads_iso_after_downsample = max_reads_iso // downsample_factor_iso
    max_reads_anc_after_downsample = max_reads_anc // downsample_factor_anc if max_reads_anc > 0 else 0
    
    # Create figure with 7 subplots (each with genetic element axes above)
    fig = plt.figure(figsize=(15, 11))
    
    # Add narrow amplicon structure axes at top using explicit positioning
    # [left, bottom, width, height] in figure coordinates
    # ax_amplicon = fig.add_axes([0.08, 0.85, 0.84, 0.04])
    # ax_amplicon.axis('off')
    # draw_amplicon_structure(ax_amplicon, h_pixels=20, y=0.1)
    
    # Junction plots grid - adjust to leave room for top axes
    gs = gridspec.GridSpec(2, 4, figure=fig, hspace=0.1, wspace=0.3,
                           top=0.95, bottom=0.05, left=0.08, right=0.92)
    
    junction_types = JunctionType.sorted()
    for jt in junction_types:
        row, col = JUNCTION_TYPES_TO_PLOT_POSITIONS[jt]
        
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
        jct_length = jct_lengths[jt]
        midpoint = jct_length // 2
        
        # Sort alignments by start position
        sorted_alignments = sorted(alignments, key=lambda x: x[0])
        
        # Downsample isolate
        downsampled = sorted_alignments[::downsample_factor_iso]
        
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
            sorted_alignments_anc = sorted(alignments_anc, key=lambda x: x[0])
            downsampled_anc = sorted_alignments_anc[::downsample_factor_anc]
            
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
        
        # Add vertical line at junction point
        ax.axvline(0, color='black', linestyle='--', linewidth=2, alpha=0.5)
        
        # Add horizontal line at y=0 if ancestor data exists
        if alignment_data_anc:
            ax.axhline(0, color='black', linestyle='-', linewidth=1, alpha=0.3)
        
        # Set consistent y-axis (extend to negative if ancestor data exists)
        y_min = -max_reads_anc_after_downsample - 1 if alignment_data_anc else 0
        y_max = max_reads_iso_after_downsample + 1
        ax.set_ylim(y_min, y_max)
        ax.set_ylabel('Read number (iso +, anc -)' if alignment_data_anc else 'Read number')
        
        # Set x-axis limits centered at junction
        ax.set_xlim(-midpoint, midpoint)
        ax.set_xlabel('Position relative to junction (bp)')
        ax.grid(True, alpha=0.3, axis='x')
        
        # Draw genetic elements in the narrow axes above
        left_elem_type, right_elem_type = JUNCTION_TYPES_TO_LEFT_RIGHT_ELEMENTS[jt]
        
        # Set up genetic element axes
        ax_gene.set_xlim(-midpoint, midpoint)
        ax_gene.set_ylim(0, 1)
        ax_gene.set_aspect('auto')
        
        # Add title on top of genetic element axes
        # ax_gene.set_title(f'{jt.name}

        if alignment_data_anc:
            ax.text(0.05, 0.95, f'iso (n={len(alignments)})', fontsize=10, transform=ax.transAxes, ha='left', va='top')
            ax.text(0.05, 0.05, f'anc (n={len(alignment_data_anc[jt])})', fontsize=10, transform=ax.transAxes, ha='left', va='bottom')
        else:
            ax.text(0.05, 0.95, f'n={len(alignments)}', fontsize=10, transform=ax.transAxes, ha='left', va='top')

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
                        GeneticElementType.CHROMOSOME, h_pixels=h_pixels, wave_tail=True)
    draw_genetic_element(ax_legend, y_start, right_start, right_end, 
                        GeneticElementType.CHROMOSOME, h_pixels=h_pixels, wave_head=True)
    ax_legend.text(label_x, y_start, 'chromosome', fontsize=10, va='bottom')
    
    # Amplicon
    y_amp = y_start - element_height - 0.05
    draw_genetic_element(ax_legend, y_amp, left_start, left_end, 
                        GeneticElementType.AMPLICON, h_pixels=h_pixels, wave_head=True)
    draw_genetic_element(ax_legend, y_amp, right_start, right_end, 
                        GeneticElementType.AMPLICON, h_pixels=h_pixels, wave_tail=True)
    ax_legend.text(label_x, y_amp, 'amplicon', fontsize=10, va='bottom')
    
    # IS element
    y_is = y_amp - element_height - 0.05
    draw_genetic_element(ax_legend, y_is, left_start, left_end, 
                        GeneticElementType.TN_ELEMENT, h_pixels=h_pixels, wave_head=True)
    draw_genetic_element(ax_legend, y_is, right_start, right_end, 
                        GeneticElementType.TN_ELEMENT, h_pixels=h_pixels, wave_tail=True)
    ax_legend.text(label_x, y_is, 'IS', fontsize=10, va='bottom')
    
    # Read type legend at bottom (combine left/right since they share the same color)
    legend_elements = [
        Line2D([0], [0], color=READ_TYPES_TO_COLORS['spanning'], linewidth=2, label='spanning'),
        Line2D([0], [0], color=READ_TYPES_TO_COLORS['left'], linewidth=2, label='left/right'),
        Line2D([0], [0], color=READ_TYPES_TO_COLORS['undetermined'], linewidth=2, label='undetermined')
    ]
    
    # Set legend title based on downsampling
    legend_title = None
    if downsample_factor_iso > 1 or downsample_factor_anc > 1:
        if downsample_factor_iso == downsample_factor_anc:
            legend_title = f'downsampled 1:{downsample_factor_iso}'
        else:
            legend_title = f'iso 1:{downsample_factor_iso}, anc 1:{downsample_factor_anc}'
    
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
    
    plot_alignment_coverage(
        alignment_data=alignment_data,
        jct_lengths=jct_lengths,
        title='Alignment Coverage Demo (No Reads)',
        output_path='alignment_coverage_demo.png'
    )

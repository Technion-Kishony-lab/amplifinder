"""Plotting utilities for read alignment visualization."""
import numpy as np

from amplifinder.optional_deps import plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D

from pathlib import Path
from typing import Dict, List, Tuple

from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.data_types.enums import Element, ReadType, JcCall
from amplifinder.steps.jct_coverage.read_type import get_expected_counts
from amplifinder.visualization.genetic_elements import draw_genetic_element


READ_TYPES_TO_COLORS = {
    ReadType.LEFT: '#d00000',  # vibrant red
    ReadType.LEFT_MARGINAL: '#b06060',  # soft coral
    ReadType.MIDDLE: '#606060',  # teal green
    ReadType.RIGHT_MARGINAL: '#6060b0',  # warm amber
    ReadType.RIGHT: '#0000d0'  # bright orange
}


JC_CALLS_TO_COLORS = {
    True: 'green',
    False: 'red',
    None: 'grey'
}


def format_read_info(jc_cov: JunctionReadCounts, scales: JunctionReadCounts) -> list[Line2D]:
    """Create legend elements with colored markers for read counts."""
    return [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=READ_TYPES_TO_COLORS[rt],
               markersize=6, label=f"{jc_cov[rt]} :{scales[rt]}")
        for rt in ReadType
    ]


def add_pie_chart(ax, jc_cov: JunctionReadCounts, position: str = 'top'):
    """Add a small pie chart showing read distribution.

    Args:
        ax: The axes to add the pie chart to
        jc_cov: JunctionReadCounts with left, left_marginal, spanning, right_marginal, right
        position: 'top' for top-left or 'bottom' for bottom-right
    """
    # Calculate fractions (order: LEFT, LEFT_MARGINAL, MIDDLE, RIGHT_MARGINAL, RIGHT)
    if jc_cov.total == 0:
        return

    sizes = list(jc_cov.counts.values())
    colors = list(READ_TYPES_TO_COLORS.values())

    # Create inset axes for pie chart (bounds in axes coordinates: [left, bottom, width, height])
    y_offset = 0.32 if position == 'top' else 0.5
    inset_ax = ax.inset_axes([0.72, y_offset, 0.28, 0.18])

    inset_ax.pie(sizes, colors=colors, startangle=90, counterclock=False,
                 wedgeprops=dict(edgecolor='white', linewidth=0.5))
    inset_ax.set_aspect('equal')


def _down_sample_alignments(
    alignments: List[Tuple[int, int, ReadType]], jc_cov: JunctionReadCounts,
    max_reads_per_plot: int, arm_len: int, min_overlap_len: int, read_len: int
) -> tuple[List[Tuple[int, int, ReadType]], JunctionReadCounts]:
    """Downsample alignments to a maximum number of reads per plot.

    Left and right alignments are scaled independently to max_reads_per_plot.
    Spanning and marginal alignments are scaled by the minimum of left/right scaling factors.

    Returns:
        Tuple of (downsampled_alignments, scales)
        where scales is a JunctionReadCounts with scaling factors for each read type
    """
    # Sort alignments by start position
    sorted_alignments = sorted(alignments, key=lambda x: x[0])

    expected = get_expected_counts(arm_len, min_overlap_len, read_len)
    max_reads = expected * max_reads_per_plot // expected.total

    scales = jc_cov // max_reads

    counters = JunctionReadCounts()

    downsampled = []
    for start, end, read_type in sorted_alignments:
        count = counters[read_type]
        scale = scales[read_type]
        if scale == 0 or count % scale == 0:
            downsampled.append((start, end, read_type))
        counters.increment(read_type)

    return downsampled, scales


def _plot_alignments(ax, reads: List[Tuple[int, int, ReadType]], arm_len: int,
                     y_sign: int = 1, alpha: float = 1.0, use_nan: bool = True) -> None:
    # Plot isolate reads
    # start, end are 0-based end-exclusive coordinates.
    #
    #  |---|---|---|---|---|---|---|---|
    # -4  -3  -2  -1   0   1   2   3   4  x-axis
    #    0   1   2   3   4   5   6   7    seq index
    #
    # junction_length = 8, arm_len = 4
    # read with start=3, end=7 will be plotted as (len = 4)
    #              |---|---|---|---|
    #             -1               3      x-axis
    # x_start = start - arm_len = 3 - 4 = -1
    # x_end = end - arm_len = 7 - 4 = 3

    segments_by_type = {read_type: {'x': [], 'y': []} for read_type in ReadType}

    for y_pos, (start, end, read_type) in enumerate(reads, start=1):
        x_start = start - arm_len
        x_end = end - arm_len
        segments_by_type[read_type]['x'].append([x_start, x_end])
        segments_by_type[read_type]['y'].append([y_sign * y_pos, y_sign * y_pos])

    for read_type, segments in segments_by_type.items():
        x_arr = np.array(segments['x'])
        y_arr = np.array(segments['y'])
        if use_nan:  # Seems a bit faster to plot with nan separators
            x_arr = np.column_stack([x_arr, np.full(len(x_arr), np.nan)]).ravel()
            y_arr = np.column_stack([y_arr, np.full(len(y_arr), np.nan)]).ravel()
        else:
            x_arr = x_arr.T
            y_arr = y_arr.T
        ax.plot(x_arr, y_arr, color=READ_TYPES_TO_COLORS[read_type],
                linewidth=1, solid_capstyle='butt', alpha=alpha)


def _draw_genetic_element_legend(ax_legend, h_pixels: int = 18) -> None:
    """Draw genetic element legend showing chromosome, amplicon, and IS elements."""
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
    min_overlap_len: int = 10,
    read_len: int = 150,
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

        jc_length = jc_lengths[jt]  # TODO: what happen when anc/iso have different lengths?
        arm_len = jc_length // 2

        # Downsample and plot isolate
        alignments = alignment_data[jt]
        jc_cov = jc_covs[jt]
        downsampled, scales_iso = _down_sample_alignments(
            alignments, jc_cov, max_reads_per_plot, arm_len, min_overlap_len, read_len
        )
        _plot_alignments(ax, downsampled, arm_len, y_sign=1)

        # Downsample and plot ancestor reads below x-axis (negative y) if available
        if alignment_data_anc:
            alignments_anc = alignment_data_anc[jt]
            jc_cov_anc = jc_covs_anc[jt]
            downsampled_anc, scales_anc = _down_sample_alignments(
                alignments_anc, jc_cov_anc, max_reads_per_plot, arm_len, min_overlap_len, read_len
            )
            _plot_alignments(ax, downsampled_anc, arm_len, y_sign=-1, alpha=0.7)

        # Set consistent y-axis (extend to negative if ancestor data exists)
        y_min = -(max_reads_per_plot + 1) if alignment_data_anc else 0
        y_max = max_reads_per_plot + 1
        ax.set_ylim(y_min, y_max)
        ax.set_xlim(-arm_len, arm_len)

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

        # Set x-axis limits centered at junction
        ax.set_xlabel('Position relative to junction (bp)')
        ax.set_ylabel('Read number')
        ax.grid(True, alpha=0.3, axis='x')

        # Draw genetic elements in the narrow axes above
        left_elem_type, right_elem_type = jt.elements

        # Set up genetic element axes
        ax_gene.set_xlim(-arm_len, arm_len)
        ax_gene.set_ylim(0, 1)
        ax_gene.set_aspect('auto')

        fs = 7
        if alignment_data_anc:
            iso_legend = format_read_info(jc_cov, scales_iso)
            anc_legend = format_read_info(jc_cov_anc, scales_anc)
            leg1 = ax.legend(handles=iso_legend, loc='upper left', fontsize=fs,
                             framealpha=0.8, title='iso', title_fontsize=fs)
            ax.add_artist(leg1)
            ax.legend(handles=anc_legend, loc='lower left', fontsize=fs,
                      framealpha=0.8, title='anc', title_fontsize=fs)

            add_pie_chart(ax, jc_cov, position='top')
            add_pie_chart(ax, jc_cov_anc, position='bottom')
        else:
            iso_legend = format_read_info(jc_cov, scales_iso)
            ax.legend(handles=iso_legend, loc='upper left', fontsize=fs,
                      framealpha=0.8)
            # Add pie chart for isolate only
            add_pie_chart(ax, jc_cov, position='top')

        # Draw left element (ending at x=0)
        draw_genetic_element(ax_gene, 0.1, -arm_len, 0, left_elem_type, h_pixels=h_pixels, wave_tail=True)

        # Draw right element (starting at x=0)
        draw_genetic_element(ax_gene, 0.1, 0, arm_len, right_elem_type, h_pixels=h_pixels, wave_head=True)

    # Add legend in empty subplot position (2, 4)
    ax_legend = fig.add_subplot(gs[1, 3])
    _draw_genetic_element_legend(ax_legend, h_pixels=h_pixels)

    # Read-type legend (top)
    read_type_legend_elements = [
        Line2D([0], [0], color=READ_TYPES_TO_COLORS[rt], linewidth=2, label=rt.value)
        for rt in ReadType
    ]
    leg_read = ax_legend.legend(handles=read_type_legend_elements, bbox_to_anchor=(0.5, 0.35),
                                loc='center', fontsize=11, frameon=True, title='Read type')
    ax_legend.add_artist(leg_read)

    # Junction call legend (bottom)
    if with_calls:
        jc_legend_elements = [
            Line2D([0], [0], color='green', linestyle='--', linewidth=2.5, label='Positive'),
            Line2D([0], [0], color='red', linestyle='--', linewidth=2.5, label='Negative'),
            Line2D([0], [0], color='grey', linestyle='--', linewidth=2.5, label='Undetermined'),
        ]
        ax_legend.legend(handles=jc_legend_elements, bbox_to_anchor=(0.5, 0.05),
                         loc='center', fontsize=11, frameon=True, title='Junction call')

    # Add overall title
    fig.suptitle(title, fontsize=14, fontweight='bold')

    # Save figure
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)

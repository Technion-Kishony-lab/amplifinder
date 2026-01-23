"""Plot junction coverage figures combining genetic elements and alignments."""

from amplifinder.optional_deps import plt
if plt is not None:
    from matplotlib import gridspec
    from matplotlib.lines import Line2D
    from matplotlib.axes import Axes


from collections import defaultdict
from typing import Dict, List, Optional
from pathlib import Path

from amplifinder.env import PLOT_ALIGNMENT_SNP_INDELS

from amplifinder.config import AlignmentClassifyParams
from amplifinder.data_types import JunctionReadCounts, JunctionType
from amplifinder.data_types.enums import Element, ReadType, JcCall

from amplifinder.steps.jct_coverage.alignment_data import AlignmentData
from amplifinder.steps.jct_coverage.read_type import get_expected_counts
from amplifinder.steps.jct_coverage.alignment_segments import AlignmentElements, CoordsAlignmentElements, \
    convert_coords_to_nan_separated_arrays, get_alignment_segments

from amplifinder.visualization.genetic_elements import draw_genetic_element


READ_TYPES_TO_COLORS = {
    ReadType.LEFT: '#d00000',  # vibrant red
    ReadType.LEFT_MARGINAL: '#b06060',  # soft coral
    ReadType.MIDDLE: '#606060',  # teal green
    ReadType.PAIRED: '#606060',  # green for paired-end reads
    ReadType.RIGHT_MARGINAL: '#6060b0',  # warm amber
    ReadType.RIGHT: '#0000d0'  # bright orange
}


JC_CALLS_TO_COLORS: dict[Optional[bool], str] = {
    True: 'green',
    False: 'red',
    None: 'grey'
}


ALIGNMENT_ELEMENTS_PLOT_KWARGS = AlignmentElements(
    match={"color": "#1f77b4", "linewidth": 1},
    mismatch={"color": "#1f77b4", "linewidth": 1},
    deletion={"color": "#d62728", "linewidth": 3},
    insertion={"color": "#2ca02c", "marker": "v", "markersize": 4},
    snp={"color": "#1f77b4", "linewidth": 0.25},
)


COMMON_PLOT_KWARGS = {
    "solid_capstyle": "butt",
    "marker": "None",
    "color": "#1f77b4",
    "linewidth": 1,
}


def plot_alignments(ax: Axes, alignments: list[AlignmentData], arm_len: int,
                    show_events: bool = False,
                    y_step: float = 1.0,
                    read_type_to_color: dict = None,
                    alignment_elements_plot_kwargs: AlignmentElements = None,
                    common_plot_kwargs: dict = None) -> None:
    """Plot multiple alignments stacked vertically, colored by read_type."""

    read_type_to_color = read_type_to_color or {}
    alignment_elements_plot_kwargs = alignment_elements_plot_kwargs or ALIGNMENT_ELEMENTS_PLOT_KWARGS
    common_plot_kwargs = common_plot_kwargs or {}

    # Group alignments by read_type, keeping original index for y-position
    alignments_by_type = defaultdict(list)
    for i, alignment in enumerate(alignments):
        alignments_by_type[alignment.read_type].append((i, alignment))

    # Plot each read_type group
    for read_type, indexed_alignments in alignments_by_type.items():
        all_segments = CoordsAlignmentElements.create()

        # Collect segments from all alignments of this type
        for i, alignment in indexed_alignments:
            y = i * y_step
            segments = get_alignment_segments(alignment, show_events, -arm_len, y)
            all_segments.extend(segments)

        # Determine 'match' color for this read_type
        match_color = read_type_to_color.get(read_type)

        # Plot all segments for this read_type
        for key, segment in all_segments.items():
            x, y = convert_coords_to_nan_separated_arrays(segment)

            plot_kwargs_combined = {
                **COMMON_PLOT_KWARGS,
                **common_plot_kwargs,
                **ALIGNMENT_ELEMENTS_PLOT_KWARGS[key],
                **alignment_elements_plot_kwargs[key],
            }
            if match_color and key == 'match':
                plot_kwargs_combined['color'] = match_color
            ax.plot(x, y, **plot_kwargs_combined)


def add_hit_legend_with_info(ax: Axes, jc_cov: JunctionReadCounts, scales: JunctionReadCounts,
                             loc: str = 'upper left', title: str = None, fontsize: int = 7):
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=READ_TYPES_TO_COLORS[rt],
               markersize=6, label=f"{jc_cov[rt]} :{scales[rt]}")
        for rt in ReadType
    ]
    return ax.legend(handles=legend_elements, loc=loc, fontsize=fontsize,
                     framealpha=0.8, title=title, title_fontsize=fontsize)


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

    sizes_to_colors = {jc_cov.counts[rt]: READ_TYPES_TO_COLORS[rt] for
                       rt in ReadType}

    # Create inset axes for pie chart (bounds in axes coordinates: [left, bottom, width, height])
    y_offset = 0.32 if position == 'top' else 0.5
    inset_ax = ax.inset_axes([0.72, y_offset, 0.28, 0.18])

    inset_ax.pie(list(sizes_to_colors), colors=list(sizes_to_colors.values()), startangle=90, counterclock=False,
                 wedgeprops=dict(edgecolor='white', linewidth=0.5))
    inset_ax.set_aspect('equal')


def _down_sample_alignments(
    alignments: List[AlignmentData],
    jc_cov: JunctionReadCounts,
    max_reads_per_plot: int, expected_counts: JunctionReadCounts
) -> tuple[List[AlignmentData], JunctionReadCounts]:
    """Downsample alignments to a maximum number of reads per plot.

    Left and right alignments are scaled independently to max_reads_per_plot.
    Spanning and marginal alignments are scaled by the minimum of left/right scaling factors.

    Args:
        alignments: List of SingleAlignment or PairedAlignment tuples
        jc_cov: JunctionReadCounts with actual counts
        max_reads_per_plot: Maximum reads to show per plot
        expected_counts: Expected counts for calculating scaling factors

    Returns:
        Tuple of (downsampled_alignments, scales)
        where scales is a JunctionReadCounts with scaling factors for each read type
    """
    sorted_alignments = sorted(alignments, key=lambda a: a.left)

    max_reads = expected_counts * max_reads_per_plot // expected_counts.total

    scales = jc_cov // max_reads

    counters = JunctionReadCounts()

    downsampled = []
    for alignment in sorted_alignments:
        read_type = alignment.read_type
        count = counters[read_type]
        scale = scales[read_type]
        if scale == 0 or count % scale == 0:
            downsampled.append(alignment)
        counters.increment(read_type)

    return downsampled, scales


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
    alignment_data: Dict[JunctionType, List[AlignmentData]],
    alignment_data_anc: Dict[JunctionType, List[AlignmentData]] | None = None,
    jc_covs: Dict[JunctionType, JunctionReadCounts] | None = None,
    jc_covs_anc: Dict[JunctionType, JunctionReadCounts] | None = None,
    jc_calls: Dict[JunctionType, JcCall] | None = None,
    jc_calls_anc: Dict[JunctionType, JcCall] | None = None,
    title: str | None = None,
    output_path: Path | str = 'junctions_coverage.png',
    max_reads_per_plot: int = 200,
    alignment_classify_params: AlignmentClassifyParams | None = None,
    iso_read_len: int = 150,
    anc_read_len: int = 150,
) -> None:
    """Plot junctions coverage by reads as horizontal lines (coverage plot)."""
    with_calls = jc_calls is not None
    if with_calls:
        if alignment_data_anc and not jc_calls_anc:
            raise ValueError("Cannot plot with junction calls for isolate and no ancestor calls")

    h_pixels = 18

    alignment_classify_params = alignment_classify_params or AlignmentClassifyParams()

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
        expected_iso = get_expected_counts(iso_read_len, arm_len, alignment_classify_params)
        downsampled, scales_iso = _down_sample_alignments(
            alignments, jc_cov, max_reads_per_plot, expected_iso
        )
        plot_alignments(
            ax, downsampled, arm_len,
            show_events=PLOT_ALIGNMENT_SNP_INDELS,
            y_step=1.0,
            read_type_to_color=READ_TYPES_TO_COLORS,
        )

        # Downsample and plot ancestor reads below x-axis (negative y) if available
        if alignment_data_anc:
            alignments_anc = alignment_data_anc[jt]
            jc_cov_anc = jc_covs_anc[jt]
            expected_anc = get_expected_counts(anc_read_len, arm_len, alignment_classify_params)
            downsampled_anc, scales_anc = _down_sample_alignments(
                alignments_anc, jc_cov_anc, max_reads_per_plot, expected_anc
            )
            plot_alignments(
                ax, downsampled_anc, arm_len,
                show_events=PLOT_ALIGNMENT_SNP_INDELS,
                y_step=-1.0,
                read_type_to_color=READ_TYPES_TO_COLORS,
            )

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
        max_read_len = max(iso_read_len, anc_read_len)
        x_lim = max_read_len + alignment_classify_params.max_dist_from_junction
        ax_gene.set_xlim(-x_lim, x_lim)
        ax_gene.set_ylim(0, 1)
        ax_gene.set_aspect('auto')

        fs = 7
        if alignment_data_anc:
            leg1 = add_hit_legend_with_info(ax, jc_cov, scales_iso, loc='upper left',
                                            title='iso', fontsize=fs)
            ax.add_artist(leg1)
            add_hit_legend_with_info(ax, jc_cov_anc, scales_anc, loc='lower left',
                                     title='anc', fontsize=fs)

            add_pie_chart(ax, jc_cov, position='top')
            add_pie_chart(ax, jc_cov_anc, position='bottom')
        else:
            add_hit_legend_with_info(ax, jc_cov, scales_iso, loc='upper left', fontsize=fs)
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

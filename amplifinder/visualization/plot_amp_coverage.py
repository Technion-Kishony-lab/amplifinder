"""Plotting utilities for coverage visualization."""
import numpy as np
from amplifinder.optional_deps import plt
from matplotlib.axes import Axes
from pathlib import Path
from typing import Optional, Dict

from amplifinder.data_types import AnalyzedTnJc2, JunctionType, JcCall
from amplifinder.visualization.genetic_elements import draw_genetic_element
from amplifinder.data_types import Element


VLINE_STYLE = {'color': 'red', 'linestyle': '--', 'linewidth': 1.5, 'alpha': 0.7}
CONNECTOR_STYLE = {'color': 'black', 'linestyle': ':', 'linewidth': 1}
JC_ARM_RELATIVE_WIDTH = 0.05
SHIFT_MULTIPLIER = 1.05

H_PIXELS_SCHEM = 20


def _create_figure_layout(has_ancestor: bool):
    """Create figure and axes layout based on whether ancestor data is present.

    Returns:
        tuple: (fig, ax_schem_iso, ax_cov_iso, ax_schem_anc, ax_cov_anc, ax_cov_ratio)
               Last three are None if no ancestor data.
    """
    if has_ancestor:
        # 5-panel layout with gaps only between groups
        fig = plt.figure(figsize=(12, 14))

        # Define heights and gaps
        h_schem = 0.07
        h_cov = 0.21
        gap = 0.01

        # Calculate positions from top down
        top1 = 0.95
        bottom1 = top1 - h_schem - h_cov

        top2 = bottom1 - gap
        bottom2 = top2 - h_schem - h_cov

        top3 = bottom2 - gap
        bottom3 = top3 - h_cov

        # Iso group (schematic + coverage) - no gap
        gs_iso = fig.add_gridspec(2, 1, height_ratios=[h_schem, h_cov],
                                  hspace=0, top=top1, bottom=bottom1)
        ax_schem_iso = fig.add_subplot(gs_iso[0])
        ax_cov_iso = fig.add_subplot(gs_iso[1], sharex=ax_schem_iso)

        # Anc group (schematic + coverage) - no gap
        gs_anc = fig.add_gridspec(2, 1, height_ratios=[h_schem, h_cov],
                                  hspace=0, top=top2, bottom=bottom2)
        ax_schem_anc = fig.add_subplot(gs_anc[0], sharex=ax_schem_iso)
        ax_cov_anc = fig.add_subplot(gs_anc[1], sharex=ax_schem_iso)

        # Ratio panel with gap above
        gs_ratio = fig.add_gridspec(1, 1, top=top3, bottom=bottom3)
        ax_cov_ratio = fig.add_subplot(gs_ratio[0], sharex=ax_schem_iso)

        return fig, ax_schem_iso, ax_cov_iso, ax_schem_anc, ax_cov_anc, ax_cov_ratio
    else:
        # 2-panel layout: schematic_iso, coverage_iso
        fig = plt.figure(figsize=(12, 8))
        gs = fig.add_gridspec(2, 1, height_ratios=[0.78, 2.5], hspace=0, top=0.96)
        ax_schem_iso = fig.add_subplot(gs[0])
        ax_cov_iso = fig.add_subplot(gs[1], sharex=ax_schem_iso)

        return fig, ax_schem_iso, ax_cov_iso, None, None, None


def _format_coverage_axes(
    ax: Axes,
    tnjc2: AnalyzedTnJc2,
    plot_start: int,
    scaf: str,
    show_xlabel: bool,
    ref_tns_pos: list[tuple[int, int]],
    ylabel: str = 'Coverage, scaffold-normalized',
) -> None:
    left_pos = tnjc2.left - plot_start
    right_pos = tnjc2.right - plot_start

    # Axes formatting
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_yscale('log')

    # Add reference line at 1.0
    ax.axhline(1.0, color='grey', linestyle='--', linewidth=0.5)

    # Add amplicon boundaries (dashed red lines)
    ax.axvline(left_pos, **VLINE_STYLE)
    ax.axvline(right_pos, **VLINE_STYLE)
    for tn_start, tn_end in ref_tns_pos:
        ax.axvline(tn_start - plot_start, **VLINE_STYLE)
        ax.axvline(tn_end - plot_start, **VLINE_STYLE)

    # Add tick marks at amplicon boundaries
    ax.set_xticks([left_pos, right_pos])
    ax.set_xticklabels([str(tnjc2.left), str(tnjc2.right)])

    # Conditionally hide tick labels (must use tick_params for sharex axes)
    if not show_xlabel:
        ax.tick_params(labelbottom=False)

    if show_xlabel:
        ax.set_xlabel(f"Genomic position, {scaf}", fontsize=10)


def _plot_coverage_with_reference(
    ax: Axes,
    positions: np.ndarray,
    coverage: np.ndarray,
    ylim: tuple,
    left_pos: int,
    right_pos: int,
    ref_copy_number: float,
) -> None:
    ax.set_ylim(ylim)
    ax.plot(positions, coverage, 'k-', linewidth=1, zorder=1)
    ax.hlines(ref_copy_number, left_pos, right_pos,
              colors='blue', linestyles='--', linewidth=2, alpha=0.8, zorder=2)


def _draw_single_junction(
    ax: Axes,
    jt: JunctionType,
    pos: int,
    plot_start: int, plot_end: int,
    y: float,
    h_pixels: int = H_PIXELS_SCHEM,
    y_target: Optional[float] = None, pos_target: Optional[int] = None,
) -> None:
    element_left, element_right = jt.element_pair
    plot_width = plot_end - plot_start
    element_width = plot_width * JC_ARM_RELATIVE_WIDTH

    # Draw genetic elements
    draw_genetic_element(ax, y, pos - element_width - plot_start, pos - plot_start,
                         element_left, h_pixels=h_pixels, wave_tail=True, zorder=2)
    draw_genetic_element(ax, y, pos - plot_start, pos + element_width - plot_start,
                         element_right, h_pixels=h_pixels, wave_head=True, zorder=2)

    # Draw diagonal connector line if shifted
    if y_target is not None and pos_target is not None:
        ax.plot([pos - plot_start, pos_target - plot_start], [y, y_target],
                **CONNECTOR_STYLE, zorder=1)


def _draw_reference_elements(
    ax: Axes,
    tnjc2: AnalyzedTnJc2,
    plot_start: int,
    plot_end: int,
    y_zero: float,
    h_pixels: int = H_PIXELS_SCHEM
) -> None:
    rjv = tnjc2.rudimentary_junction_values
    jcs_positions = rjv.get_junction_chr_positions()

    def draw(start, end, element, **kwargs):
        draw_genetic_element(ax, y_zero, start - plot_start, end - plot_start,
                             element, h_pixels=h_pixels, zorder=2, **kwargs)

    draw(plot_start, jcs_positions[JunctionType.CHR_AMP], Element.CHR, wave_tail=True)

    if tnjc2.tnjc_left.is_ref_tn_junction():
        draw(jcs_positions[JunctionType.CHR_TN], jcs_positions[JunctionType.TN_AMP], Element.TN)

    draw(jcs_positions[JunctionType.TN_AMP], jcs_positions[JunctionType.AMP_TN], Element.AMP)

    if tnjc2.tnjc_right.is_ref_tn_junction():
        draw(jcs_positions[JunctionType.AMP_TN], jcs_positions[JunctionType.TN_CHR], Element.TN)

    draw(plot_end, jcs_positions[JunctionType.AMP_CHR], Element.CHR, wave_tail=True)


def _draw_schematic_panel(
    ax: Axes,
    plot_start: int,
    plot_end: int,
    tnjc2: AnalyzedTnJc2,
    jc_calls: dict[JunctionType, JcCall],
) -> None:
    # Set up the axis - xlim is already set by sharex, only set ylim
    ax.set_ylim(-0.05, 1.3)
    ax.axis('off')

    # Define vertical positions for three tiers (ZERO at bottom to align with coverage plot)
    y_zero = 0.05    # ZERO tier for reference elements (at bottom)
    y_bottom = 0.50  # BOTTOM tier for CHR-AMP / AMP-CHR
    y_top = 0.95     # TOP tier for IS-related junctions

    plot_width = plot_end - plot_start
    element_width = plot_width * JC_ARM_RELATIVE_WIDTH
    shift = element_width * SHIFT_MULTIPLIER

    _draw_reference_elements(ax, tnjc2, plot_start, plot_end, y_zero)
    rjv = tnjc2.rudimentary_junction_values
    jcs_positions = rjv.get_junction_chr_positions()

    # Calculate element height in data coordinates
    inv_trans = ax.transData.inverted()
    origin = inv_trans.transform([[0, 0]])[0]
    h_in_y_data = abs(inv_trans.transform([[0, H_PIXELS_SCHEM]])[0][1] - origin[1])

    tn_jc_pairs = {
        JunctionType.TN_AMP: JunctionType.CHR_TN,
        JunctionType.CHR_TN: JunctionType.TN_AMP,
        JunctionType.TN_CHR: JunctionType.AMP_TN,
        JunctionType.AMP_TN: JunctionType.TN_CHR,
    }

    y_target = y_zero + h_in_y_data
    for jt, jc_call in jc_calls.items():
        if jc_call is not JcCall.POS:
            continue
        if jt is JunctionType.AMP_AMP:
            continue
        pos = jcs_positions[jt]
        if jt in [JunctionType.CHR_AMP, JunctionType.AMP_CHR]:
            y = y_bottom
            pos_shifted = pos
        else:
            tn_jc_pair = tn_jc_pairs[jt]
            pos_pair = jcs_positions[tn_jc_pair]
            mid_pos = (pos_pair + pos) / 2
            y = y_top
            shift_sign = 1 if jt in [JunctionType.TN_AMP, JunctionType.TN_CHR] else -1
            pos_shifted = mid_pos + shift_sign * shift
        _draw_single_junction(ax, jt, pos_shifted, plot_start, plot_end,
                              y=y, y_target=y_target, pos_target=pos)


def plot_amplicon_coverage(
    tnjc2: AnalyzedTnJc2,
    jc_calls: dict[JunctionType, JcCall],
    jc_calls_anc: dict[JunctionType, JcCall],
    iso_scafs_to_covs: Dict[str, np.ndarray],
    anc_scafs_to_covs: Optional[Dict[str, np.ndarray]],
    output_path: Path,
    flank_fraction: float = 0.2
) -> None:
    """Plot isolate and ancestor coverage along an amplicon region with junction schematics.

    Args:
        iso_coverage: Dictionary mapping scaffold names to isolate coverage arrays
        anc_coverage: Dictionary mapping scaffold names to ancestor coverage arrays (optional)
        tnjc2: AnalyzedTnJc2 record containing amplicon information
        output_path: Path to save plot
    """
    scaf_obj = tnjc2.scaffold

    # Get reference transposon positions
    ref_tns_pos = []
    if ref_tn := tnjc2.tnjc_left.ref_tn:
        ref_tns_pos.append((ref_tn.start, ref_tn.end))

    if ref_tn := tnjc2.tnjc_right.ref_tn:
        ref_tns_pos.append((ref_tn.start, ref_tn.end))

    all_tn_pos = [start for start, end in ref_tns_pos] + [end for start, end in ref_tns_pos]

    left_most = min([tnjc2.left, *all_tn_pos])
    right_most = max([tnjc2.right, *all_tn_pos])

    # Calculate plot range (scaffold-relative, 1-based)
    flank = int((right_most - left_most) * flank_fraction)
    plot_start = left_most - flank
    plot_end = right_most + flank
    if not scaf_obj.is_circular:
        plot_start = max(plot_start, 1)
        plot_end = min(plot_end, scaf_obj.length)
    plot_width = plot_end - plot_start

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

    # Create plot layout
    has_ancestor = anc_plot_cov is not None
    fig, ax_schem_iso, ax_cov_iso, ax_schem_anc, ax_cov_anc, ax_cov_ratio = \
        _create_figure_layout(has_ancestor)

    # Set xlim and draw schematics
    ax_schem_iso.set_xlim(0, plot_width)
    _draw_schematic_panel(ax_schem_iso, plot_start, plot_end, tnjc2, jc_calls)
    if ax_schem_anc is not None:
        _draw_schematic_panel(ax_schem_anc, plot_start, plot_end, tnjc2, jc_calls_anc)

    # Add figure title (above everything)
    title = (f"{tnjc2.raw_event.description} | {tnjc2.left}-{tnjc2.right} | "
             f"copy_number: {tnjc2.copy_number:.1f}x")
    fig.suptitle(title, fontsize=12, fontweight='bold')

    # Calculate shared y-axis limits for iso and anc using 99th percentile
    if anc_plot_cov is None:
        ymax = np.percentile(iso_plot_cov, 99)
        ymin = np.percentile(iso_plot_cov, 1)
    else:
        ymax = max(np.percentile(iso_plot_cov, 99), np.percentile(anc_plot_cov, 99))
        ymin = min(np.percentile(iso_plot_cov, 1), np.percentile(anc_plot_cov, 1))
    ymin = max(ymin / 1.1, 1e-1)
    ylim = (ymin, ymax * 1.1)

    # Calculate amplicon boundaries in plot coordinates
    left_pos = tnjc2.left - plot_start
    right_pos = tnjc2.right - plot_start

    # Plot isolate coverage
    _format_coverage_axes(
        ax_cov_iso, tnjc2, plot_start, tnjc2.scaf,
        show_xlabel=(ax_cov_anc is None and ax_cov_ratio is None),
        ref_tns_pos=ref_tns_pos)
    _plot_coverage_with_reference(ax_cov_iso, positions, iso_plot_cov, ylim,
                                  left_pos, right_pos, tnjc2.iso_scaf_norm_copy_number)

    # Plot ancestor coverage
    if ax_cov_anc is not None and anc_plot_cov is not None:
        _format_coverage_axes(ax_cov_anc, tnjc2, plot_start, tnjc2.scaf,
                              show_xlabel=False, ref_tns_pos=ref_tns_pos)
        _plot_coverage_with_reference(ax_cov_anc, positions, anc_plot_cov, ylim,
                                      left_pos, right_pos,
                                      tnjc2.anc_scaf_norm_copy_number)

    # Plot ratio iso_cov / anc_cov
    if has_ancestor:
        _format_coverage_axes(ax_cov_ratio, tnjc2, plot_start, tnjc2.scaf,
                              show_xlabel=True, ref_tns_pos=ref_tns_pos,
                              ylabel='Isolate/Ancestor coverage ratio')
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio = iso_plot_cov / anc_plot_cov
        ylim = (1e-1, np.max(ratio[np.isfinite(ratio)]))
        _plot_coverage_with_reference(
            ax_cov_ratio, positions, ratio, ylim,
            left_pos, right_pos, tnjc2.scaf_norm_copy_number_ratio)

    plt.tight_layout()
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

import numpy as np

import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Polygon, Rectangle

from enum import Enum


class ArrowHead(str, Enum):
    """Arrow head style."""
    BLUNT = "blunt"                      # ======|
    TRIANGLE = "triangle"                # ======>
    INNER_TRIANGLE = "inner_triangle"    # ======<
    ROUND = "round"                      # ======)
    INNER_ROUND = "inner_round"          # ======(
    WAVE = "wave"                        # ======s


y = np.linspace(0, 1, 15)

ARROW_HEAD_TO_VERTICES = {
    ArrowHead.BLUNT: np.array([[0, 0], [0, 1]]),
    ArrowHead.TRIANGLE: np.array([[-1, 0], [0, 0.5], [-1, 1]]),
    ArrowHead.INNER_TRIANGLE: np.array([[0, 0], [-1, 0.5], [0, 1]]),
    ArrowHead.ROUND: np.array([-1 + np.sin(y * np.pi), y]).T,
    ArrowHead.INNER_ROUND: np.array([-np.sin(y * np.pi), y]).T,
    ArrowHead.WAVE: np.array([-1 + np.sin(y * 5 * np.pi), y]).T,
}


class GeneticElementType(str, Enum):
    CHROMOSOME = "chromosome"
    AMPLICON = "amplicon"
    TN_ELEMENT = "tn_element"


GENETIC_ELEMENT_TYPE_TO_ARROW_PARAMS = {
    GeneticElementType.CHROMOSOME: {
        'color': 'lightgrey',
        'head': ArrowHead.BLUNT,
        'tail': ArrowHead.BLUNT,
    },
    GeneticElementType.AMPLICON: {
        'color': 'lightblue',
        'head': ArrowHead.ROUND,
        'tail': ArrowHead.BLUNT,
    },
    GeneticElementType.TN_ELEMENT: {
        'color': 'goldenrod',
        'head': ArrowHead.TRIANGLE,
        'tail': ArrowHead.BLUNT,
    },
}


def draw_horizontal_arrow(ax, y, x1, x2, h_pixels, color, 
        head: ArrowHead | str = ArrowHead.BLUNT, tail: ArrowHead | str = ArrowHead.BLUNT, 
        head_width_ratio: float = 0.5, tail_width_ratio: float = 0.5):
    """Draw a horizontal arrow.
    
    Args:
        ax: matplotlib axes object
        y: y-coordinate (bottom of element)
        x1: start x-coordinate (in sequence coordinates)
        x2: end x-coordinate (in sequence coordinates)
        h_pixels: height (in pixels)
        color: color of the element
        head: arrow head style at x2 end
        tail: arrow head style at x1 end (default BLUNT)
        arrow_width_ratio: arrow width as fraction of element height (default 0.5)
    """
    head = ArrowHead(head) if isinstance(head, str) else head
    tail = ArrowHead(tail) if isinstance(tail, str) else tail

    length = x2 - x1
    
    # Convert h_pixels from pixels to data coordinates
    inv_trans = ax.transData.inverted()
    y_bottom_display = ax.transData.transform([[0, y]])[0][1]
    y_top_display = y_bottom_display + h_pixels
    y_top_data = inv_trans.transform([[0, y_top_display]])[0][1]
    h = y_top_data - y
    
    # Calculate arrow widths in data coordinates
    x_origin = inv_trans.transform([[0, 0]])[0][0]
    h_in_x_data = abs(inv_trans.transform([[h_pixels, 0]])[0][0] - x_origin)
        
    # Build vertices: tail (bottom-up) -> head (bottom->top) -> close
    coords1 = ARROW_HEAD_TO_VERTICES[tail][::-1]  # reverse to allow concatenation with head vertices
    coords2 = ARROW_HEAD_TO_VERTICES[head]

    dir = np.sign(length)
    
    coords1 = coords1 * np.array([[-h_in_x_data * tail_width_ratio * dir, h]]) + np.array([[x1, y]])
    coords2 = coords2 * np.array([[+h_in_x_data * head_width_ratio * dir, h]]) + np.array([[x2, y]])

    vertices = np.vstack([coords1, coords2])
    polygon = Polygon(vertices, facecolor=color, edgecolor='black', linewidth=0.5)
    ax.add_patch(polygon)


def draw_genetic_element(ax, y, x1, x2, element_type: GeneticElementType, h_pixels=10, 
                         wave_tail: bool = False, wave_head: bool = False, **kwargs):
    """Draw a genetic element.
    
    Args:
        ax: matplotlib axes object
        y: y-coordinate (bottom of element)
        x1: start x-coordinate (in sequence coordinates)
        x2: end x-coordinate (in sequence coordinates)
        element_type: type of genetic element
        wave_tail: if True, add a wave to the tail
        wave_head: if True, add a wave to the head
    """
    params = GENETIC_ELEMENT_TYPE_TO_ARROW_PARAMS[element_type].copy()
    if wave_tail:
        params['tail'] = ArrowHead.WAVE
        params['tail_width_ratio'] = 0.2
    if wave_head:
        params['head'] = ArrowHead.WAVE
        params['head_width_ratio'] = 0.2
    params.update(kwargs)
    draw_horizontal_arrow(ax, y, x1, x2, h_pixels, **params)


def draw_amplicon_structure(
    ax,
    chr_left_len: float = 120,
    amp_len: float = 100,
    tn_len: float = 60,
    h_pixels: float = 40,
    y: float = 0,
    chr_right_len: float = 120,
    n_amplicon_copies: int = 2,
) -> float:
    """Draw amplicon structure into existing axes: ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~
    
    Args:
        ax: matplotlib axes to draw into
        chr_left_len: Length of left chromosome region (in sequence coordinates)
        amp_len: Length of amplicon (in sequence coordinates)
        tn_len: Length of transposon (in sequence coordinates)
        h_pixels: Height of elements in pixels
        y: Y-coordinate (bottom of element) in data units
        chr_right_len: Length of right chromosome region (in sequence coordinates)
        n_amplicon_copies: Number of amplicon copies (each contains TN)
    """
    # Calculate total length
    total_length = chr_left_len + tn_len + n_amplicon_copies * (amp_len + tn_len) + chr_right_len
    
    # Set axis limits BEFORE drawing (needed for arrow width calculation)
    ax.set_xlim(0, total_length)
    ax.set_ylim(0, 10)
    
    # Build structure from left to right
    x = 0
    
    # Left chromosome
    draw_genetic_element(ax, y, x, (x := x + chr_left_len), GeneticElementType.CHROMOSOME, h_pixels, wave_tail=True)
    
    # TN element
    draw_genetic_element(ax, y, x, (x := x + tn_len), GeneticElementType.TN_ELEMENT, h_pixels)

    # Amplicon copies (each with TN inside)
    for i in range(n_amplicon_copies):
        # Amplicon
        draw_genetic_element(ax, y, x, (x := x + amp_len), GeneticElementType.AMPLICON, h_pixels)
        
        # TN element
        draw_genetic_element(ax, y, x, (x := x + tn_len), GeneticElementType.TN_ELEMENT, h_pixels)
        
    # Right chromosome
    draw_genetic_element(ax, y, x, (x := x + chr_right_len), GeneticElementType.CHROMOSOME, h_pixels, wave_head=True)

    return total_length


def plot_amplicon_structure(
    output_path: Path | str,
    chr_left_len: float = 70,
    amp_len: float = 100,
    tn_len: float = 60,
    h_pixels: float = 40,
    y: float = 0,
    chr_right_len: float = 200,
    n_amplicon_copies: int = 3,
    figsize: tuple[float, float] = (12, 3)
) -> None:
    """Plot amplicon structure: ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~
    
    Args:
        output_path: Path to save the PNG file
        chr_left_len: Length of left chromosome region (in sequence coordinates)
        amp_len: Length of amplicon (in sequence coordinates)
        tn_len: Length of transposon (in sequence coordinates)
        h: Height of elements in pixels
        y: Y-coordinate (bottom of element) in data units
        chr_right_len: Length of right chromosome region (in sequence coordinates)
        n_amplicon_copies: Number of amplicon copies (each contains TN)
        figsize: Figure size (width, height)
    """
    fig, ax = plt.subplots(figsize=figsize)
    draw_amplicon_structure(ax, chr_left_len, amp_len, tn_len, h_pixels, y, chr_right_len, n_amplicon_copies)
    
    # Save figure
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    # Demo: plot amplicon structure
    plot_amplicon_structure(
        output_path='amplicon_structure_demo.png',
    )
    fig, ax = plt.subplots(figsize=(12, 3))
    plt.axis([0, 1000, 0, 10])
    draw_horizontal_arrow(ax, 2, 100, 800, h_pixels=10, color='lightblue', head=ArrowHead.TRIANGLE, tail=ArrowHead.INNER_TRIANGLE)
    draw_horizontal_arrow(ax, 3, 100, 800, h_pixels=10, color='lightblue', head=ArrowHead.INNER_TRIANGLE, tail=ArrowHead.TRIANGLE)
    draw_horizontal_arrow(ax, 4, 600, 200, h_pixels=10, color='lightblue', head=ArrowHead.ROUND, tail=ArrowHead.INNER_ROUND)
    draw_horizontal_arrow(ax, 5, 600, 200, h_pixels=10, color='lightblue', head=ArrowHead.INNER_TRIANGLE, tail=ArrowHead.TRIANGLE)
    plt.savefig('genetic_elements_demo.png', dpi=150, bbox_inches='tight')
    

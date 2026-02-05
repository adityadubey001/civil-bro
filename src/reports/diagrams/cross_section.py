"""
Cross-section diagram generator using Matplotlib.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Circle, Rectangle, FancyBboxPatch
import numpy as np
from typing import Optional
import io


def generate_cross_section(
    width: float,
    depth: float,
    cover: float,
    tension_bars: int,
    tension_bar_dia: float,
    stirrup_dia: float,
    compression_bars: Optional[int] = None,
    compression_bar_dia: Optional[float] = None,
    return_figure: bool = True,
):
    """
    Generate cross-section diagram of RC beam.

    Args:
        width: Beam width in mm
        depth: Overall depth in mm
        cover: Clear cover in mm
        tension_bars: Number of tension (bottom) bars
        tension_bar_dia: Tension bar diameter in mm
        stirrup_dia: Stirrup diameter in mm
        compression_bars: Number of compression (top) bars (optional)
        compression_bar_dia: Compression bar diameter in mm (optional)
        return_figure: If True, return figure; if False, return bytes

    Returns:
        Matplotlib figure or PNG bytes
    """
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(6, 8))

    # Scale for display
    scale = 1

    # Draw concrete section
    concrete = FancyBboxPatch(
        (0, 0), width, depth,
        boxstyle="round,pad=0,rounding_size=5",
        linewidth=2, edgecolor='#2c3e50', facecolor='#ecf0f1'
    )
    ax.add_patch(concrete)

    # Draw stirrup (dashed rectangle inside)
    stirrup_offset_x = cover
    stirrup_offset_y = cover
    stirrup_width = width - 2 * cover
    stirrup_height = depth - 2 * cover

    stirrup = Rectangle(
        (stirrup_offset_x, stirrup_offset_y),
        stirrup_width, stirrup_height,
        linewidth=1.5, edgecolor='#27ae60', facecolor='none',
        linestyle='--'
    )
    ax.add_patch(stirrup)

    # Draw tension bars (bottom)
    bar_y = cover + stirrup_dia + tension_bar_dia / 2

    if tension_bars == 1:
        # Center the single bar
        bar_positions = [width / 2]
    else:
        # Distribute bars evenly
        first_bar_x = cover + stirrup_dia + tension_bar_dia / 2
        last_bar_x = width - cover - stirrup_dia - tension_bar_dia / 2
        bar_positions = np.linspace(first_bar_x, last_bar_x, tension_bars)

    for bar_x in bar_positions:
        bar = Circle(
            (bar_x, bar_y),
            tension_bar_dia / 2,
            facecolor='#e74c3c',
            edgecolor='#c0392b',
            linewidth=1
        )
        ax.add_patch(bar)

    # Draw compression bars (top) if provided
    if compression_bars and compression_bar_dia:
        comp_y = depth - cover - stirrup_dia - compression_bar_dia / 2

        if compression_bars == 1:
            comp_positions = [width / 2]
        else:
            first_comp_x = cover + stirrup_dia + compression_bar_dia / 2
            last_comp_x = width - cover - stirrup_dia - compression_bar_dia / 2
            comp_positions = np.linspace(first_comp_x, last_comp_x, compression_bars)

        for comp_x in comp_positions:
            bar = Circle(
                (comp_x, comp_y),
                compression_bar_dia / 2,
                facecolor='#3498db',
                edgecolor='#2980b9',
                linewidth=1
            )
            ax.add_patch(bar)

    # Add dimensions
    dim_offset = 30

    # Width dimension (bottom)
    ax.annotate(
        '', xy=(0, -dim_offset), xytext=(width, -dim_offset),
        arrowprops=dict(arrowstyle='<->', color='black', lw=1)
    )
    ax.text(width / 2, -dim_offset - 15, f'{width:.0f} mm',
            ha='center', va='top', fontsize=10, fontweight='bold')

    # Depth dimension (right side)
    ax.annotate(
        '', xy=(width + dim_offset, 0), xytext=(width + dim_offset, depth),
        arrowprops=dict(arrowstyle='<->', color='black', lw=1)
    )
    ax.text(width + dim_offset + 10, depth / 2, f'{depth:.0f} mm',
            ha='left', va='center', fontsize=10, fontweight='bold', rotation=90)

    # Cover dimension (left side)
    ax.annotate(
        '', xy=(-dim_offset / 2, 0), xytext=(-dim_offset / 2, cover),
        arrowprops=dict(arrowstyle='<->', color='gray', lw=0.8)
    )
    ax.text(-dim_offset / 2 - 5, cover / 2, f'{cover:.0f}',
            ha='right', va='center', fontsize=8, color='gray')

    # Labels
    ax.text(width / 2, bar_y + tension_bar_dia + 10,
            f'{tension_bars}-{int(tension_bar_dia)}φ',
            ha='center', va='bottom', fontsize=9, color='#c0392b', fontweight='bold')

    if compression_bars and compression_bar_dia:
        comp_y = depth - cover - stirrup_dia - compression_bar_dia / 2
        ax.text(width / 2, comp_y - compression_bar_dia - 10,
                f'{compression_bars}-{int(compression_bar_dia)}φ',
                ha='center', va='top', fontsize=9, color='#2980b9', fontweight='bold')

    # Stirrup label
    ax.text(width - cover - 5, depth / 2,
            f'{int(stirrup_dia)}φ stir.',
            ha='right', va='center', fontsize=8, color='#27ae60', rotation=90)

    # Legend
    legend_y = depth + 40
    ax.add_patch(Circle((20, legend_y), 5, facecolor='#e74c3c', edgecolor='#c0392b'))
    ax.text(35, legend_y, 'Tension steel', va='center', fontsize=8)

    if compression_bars:
        ax.add_patch(Circle((120, legend_y), 5, facecolor='#3498db', edgecolor='#2980b9'))
        ax.text(135, legend_y, 'Compression steel', va='center', fontsize=8)

    # Set axis limits and aspect
    margin = 60
    ax.set_xlim(-margin, width + margin)
    ax.set_ylim(-margin, depth + margin + 30)
    ax.set_aspect('equal')
    ax.axis('off')

    # Title
    ax.set_title('Cross Section', fontsize=12, fontweight='bold', pad=20)

    plt.tight_layout()

    if return_figure:
        return fig
    else:
        # Return as bytes
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()


def generate_reinforcement_layout(
    width: float,
    depth: float,
    span: float,
    tension_arrangement: str,
    stirrup_dia: int,
    stirrup_spacing_support: float,
    stirrup_spacing_center: float,
    compression_arrangement: Optional[str] = None,
    return_figure: bool = True,
):
    """
    Generate longitudinal reinforcement layout diagram.

    Args:
        width: Beam width in mm
        depth: Overall depth in mm
        span: Beam span in mm
        tension_arrangement: e.g., "4-20φ"
        stirrup_dia: Stirrup diameter in mm
        stirrup_spacing_support: Stirrup spacing at support in mm
        stirrup_spacing_center: Stirrup spacing at center in mm
        compression_arrangement: e.g., "2-16φ" (optional)
        return_figure: If True, return figure; if False, return bytes

    Returns:
        Matplotlib figure or PNG bytes
    """
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))

    # Scale for display
    scale_x = 1000 / span  # Normalize span to 1000 units
    scale_y = 1

    span_scaled = span * scale_x
    depth_scaled = depth * scale_y * 0.3  # Reduce vertical scale

    # Draw beam outline
    beam = Rectangle(
        (0, 0), span_scaled, depth_scaled,
        linewidth=2, edgecolor='#2c3e50', facecolor='#ecf0f1'
    )
    ax.add_patch(beam)

    # Draw supports
    support_size = 30
    # Left support (triangle)
    left_support = plt.Polygon(
        [(0, 0), (-support_size/2, -support_size), (support_size/2, -support_size)],
        facecolor='#95a5a6', edgecolor='#7f8c8d', linewidth=1.5
    )
    ax.add_patch(left_support)

    # Right support (triangle with roller)
    right_support = plt.Polygon(
        [(span_scaled, 0), (span_scaled - support_size/2, -support_size),
         (span_scaled + support_size/2, -support_size)],
        facecolor='#95a5a6', edgecolor='#7f8c8d', linewidth=1.5
    )
    ax.add_patch(right_support)

    # Draw roller circle under right support
    roller = Circle((span_scaled, -support_size - 8), 8,
                   facecolor='#bdc3c7', edgecolor='#7f8c8d')
    ax.add_patch(roller)

    # Draw tension bars (red line at bottom)
    bar_offset = depth_scaled * 0.15
    ax.plot([0, span_scaled], [bar_offset, bar_offset],
            color='#e74c3c', linewidth=3, label='Tension steel')
    ax.text(span_scaled / 2, bar_offset - 15, tension_arrangement,
            ha='center', va='top', fontsize=10, color='#c0392b', fontweight='bold')

    # Draw compression bars if provided (blue line at top)
    if compression_arrangement:
        comp_offset = depth_scaled * 0.85
        ax.plot([0, span_scaled], [comp_offset, comp_offset],
                color='#3498db', linewidth=2, label='Compression steel')
        ax.text(span_scaled / 2, comp_offset + 15, compression_arrangement,
                ha='center', va='bottom', fontsize=10, color='#2980b9', fontweight='bold')

    # Draw stirrups (vertical lines)
    # At supports (closer spacing)
    support_zone = span_scaled * 0.25
    n_stirrups_support = int(support_zone / (stirrup_spacing_support * scale_x))

    for i in range(n_stirrups_support):
        x = i * stirrup_spacing_support * scale_x
        ax.plot([x, x], [bar_offset + 5, depth_scaled - 5],
                color='#27ae60', linewidth=1, alpha=0.7)
        # Right side
        x_right = span_scaled - i * stirrup_spacing_support * scale_x
        ax.plot([x_right, x_right], [bar_offset + 5, depth_scaled - 5],
                color='#27ae60', linewidth=1, alpha=0.7)

    # At center (wider spacing)
    center_start = support_zone
    center_end = span_scaled - support_zone
    n_stirrups_center = int((center_end - center_start) / (stirrup_spacing_center * scale_x))

    for i in range(n_stirrups_center):
        x = center_start + i * stirrup_spacing_center * scale_x
        ax.plot([x, x], [bar_offset + 5, depth_scaled - 5],
                color='#27ae60', linewidth=0.8, alpha=0.5)

    # Add annotations
    # Span dimension
    ax.annotate(
        '', xy=(0, -support_size - 30), xytext=(span_scaled, -support_size - 30),
        arrowprops=dict(arrowstyle='<->', color='black', lw=1)
    )
    ax.text(span_scaled / 2, -support_size - 45, f'Span = {span/1000:.1f} m',
            ha='center', va='top', fontsize=10, fontweight='bold')

    # Stirrup spacing annotations
    ax.text(support_zone / 2, depth_scaled + 15,
            f'{stirrup_dia}φ @ {int(stirrup_spacing_support)} c/c',
            ha='center', va='bottom', fontsize=8, color='#27ae60')
    ax.text(span_scaled / 2, depth_scaled + 15,
            f'{stirrup_dia}φ @ {int(stirrup_spacing_center)} c/c',
            ha='center', va='bottom', fontsize=8, color='#27ae60')
    ax.text(span_scaled - support_zone / 2, depth_scaled + 15,
            f'{stirrup_dia}φ @ {int(stirrup_spacing_support)} c/c',
            ha='center', va='bottom', fontsize=8, color='#27ae60')

    # Set limits
    margin = 50
    ax.set_xlim(-margin, span_scaled + margin)
    ax.set_ylim(-support_size - 60, depth_scaled + 40)
    ax.set_aspect('equal')
    ax.axis('off')

    # Title
    ax.set_title('Longitudinal Section', fontsize=12, fontweight='bold', pad=10)

    plt.tight_layout()

    if return_figure:
        return fig
    else:
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        buf.seek(0)
        return buf.getvalue()

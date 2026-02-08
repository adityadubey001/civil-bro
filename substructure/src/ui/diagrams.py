"""Engineering diagrams for the Substructure Design app.

Generates matplotlib figures for cross-sections, elevations, P-M interaction
diagrams, and plan views.  Each public function returns a ``matplotlib.figure.Figure``
that the caller can pass to ``st.pyplot()``.

Colour conventions (matching typical design-report style):
    - Concrete: light gray fill (#d9d9d9), black outline
    - Rebar: red filled circles (#e74c3c)
    - Stirrups / spiral: green dashed (#27ae60)
    - Punching perimeters: blue dashed (#2980b9)
    - Dimension lines: black, thin
"""

from __future__ import annotations

import math
from typing import Any

import matplotlib
matplotlib.use("Agg")          # non-interactive backend for Streamlit
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ── Colours ──────────────────────────────────────────────────────────────────

_CONCRETE = "#d9d9d9"
_REBAR    = "#e74c3c"
_SPIRAL   = "#27ae60"
_PUNCH    = "#2980b9"
_DIM      = "#333333"


# ── Helpers ──────────────────────────────────────────────────────────────────

def _safe(obj: Any, attr: str, default: Any = 0.0) -> Any:
    """Null-safe attribute access."""
    if obj is None:
        return default
    return getattr(obj, attr, default)


def _add_dim_h(ax, y, x0, x1, text, offset=0.06, fontsize=8):
    """Draw a horizontal dimension line with arrows and a label."""
    ax.annotate(
        "", xy=(x1, y), xytext=(x0, y),
        arrowprops=dict(arrowstyle="<->", color=_DIM, lw=0.8),
    )
    ax.text((x0 + x1) / 2, y + offset, text, ha="center", va="bottom",
            fontsize=fontsize, color=_DIM)


def _add_dim_v(ax, x, y0, y1, text, offset=0.06, fontsize=8):
    """Draw a vertical dimension line with arrows and a label."""
    ax.annotate(
        "", xy=(x, y1), xytext=(x, y0),
        arrowprops=dict(arrowstyle="<->", color=_DIM, lw=0.8),
    )
    ax.text(x + offset, (y0 + y1) / 2, text, ha="left", va="center",
            fontsize=fontsize, color=_DIM, rotation=90)


# ═══════════════════════════════════════════════════════════════════════════════
# 1. Pier Cap Tapered Elevation
# ═══════════════════════════════════════════════════════════════════════════════

def draw_pier_cap_elevation(pc_result, geom) -> plt.Figure:
    """Side elevation of tapered pier cap cantilever.

    Parameters
    ----------
    pc_result : PierCapDesignResult
    geom : GeometryResults
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    # Dimensions from geometry (metres → normalised drawing coords)
    piercap_sec = _safe(geom, "piercap_section")
    width_m  = _safe(piercap_sec, "width", 2.85)       # longitudinal width (m)
    depth_max_m = _safe(piercap_sec, "depth", 1.75)     # max depth at pier face (m)

    # From config / result — depths in mm, convert to m for drawing
    depth_pier_mm = _safe(pc_result, "depth_at_pier", depth_max_m * 1000)
    depth_curt_mm = _safe(pc_result, "depth_at_curtailment", 750)
    cap_width_mm  = _safe(pc_result, "width", width_m * 1000)

    # Pier section dimensions (m)
    pier_sec = _safe(geom, "pier_section")
    pier_w_long = _safe(pier_sec, "width_long", 2.0)  # m

    # Pier cap transverse half-length from config
    piercap_length_trans_m = getattr(geom, "piercap_section", None)
    # Best effort: use bearing max_y to estimate cantilever length
    btl = _safe(geom, "btl_lhs")
    max_y = _safe(btl, "max_y_bearing", 3.7)
    cant_length_m = max_y + 0.5  # approx edge distance

    # Convert to mm for annotations
    depth_pier = depth_pier_mm
    depth_curt = depth_curt_mm
    cap_width  = cap_width_mm

    # Drawing in metres for simplicity
    D_pier = depth_pier / 1000.0
    D_curt = depth_curt / 1000.0
    W = cap_width / 1000.0
    L = cant_length_m
    pier_hw = pier_w_long / 2.0  # half-width of pier in longitudinal dir

    # Top of pier cap at y=0, pier face at x=0
    # Draw right-hand cantilever (from pier face outward)

    # Pier cap outline (trapezoidal)
    # Top edge: flat from x=-pier_hw to x=L
    # Bottom edge: tapers from D_pier at x=0 to D_curt at x=L
    # Left side: pier face

    xs_top = [-pier_hw, L]
    ys_top = [0, 0]
    xs_bot = [-pier_hw, 0, L, L]
    ys_bot = [-D_pier, -D_pier, -D_curt, -D_curt]

    # Draw concrete outline
    outline_x = [-pier_hw, L, L, 0, -pier_hw, -pier_hw]
    outline_y = [0, 0, -D_curt, -D_pier, -D_pier, 0]
    ax.fill(outline_x, outline_y, color=_CONCRETE, edgecolor="black", linewidth=1.5)

    # Pier body (rectangle below pier cap)
    pier_h = min(_safe(pier_sec, "width_trans", 2.0), 1.5)  # show partial pier
    pier_rect = mpatches.Rectangle(
        (-pier_hw, -D_pier), pier_w_long, -pier_h,
        linewidth=1.5, edgecolor="black", facecolor="#c0c0c0", zorder=1,
    )
    ax.add_patch(pier_rect)
    ax.text(pier_w_long / 2 - pier_hw, -D_pier - pier_h / 2, "PIER",
            ha="center", va="center", fontsize=9, fontweight="bold", color="#555")

    # Reinforcement lines (schematic)
    cover_m = (_safe(pc_result, "cover", 50)) / 1000.0
    # Top bars (tension for hogging) — dashed
    ax.plot([-pier_hw + cover_m, L - cover_m], [-cover_m, -cover_m],
            color=_REBAR, linewidth=1.5, linestyle="--", label="Top bars (tension)")
    # Bottom bars
    bot_y_pier = -(D_pier - cover_m)
    bot_y_curt = -(D_curt - cover_m)
    ax.plot([0, L - cover_m], [bot_y_pier, bot_y_curt],
            color=_REBAR, linewidth=1.5, linestyle="-", label="Bottom bars")

    # Stirrups (vertical dashed lines at intervals)
    n_stirrups = 6
    for i in range(n_stirrups + 1):
        x = i * L / n_stirrups
        depth_here = D_pier + (D_curt - D_pier) * (x / L) if L > 0 else D_pier
        ax.plot([x, x], [-cover_m, -(depth_here - cover_m)],
                color=_SPIRAL, linewidth=0.8, linestyle="--")

    # Bearing position marker
    bearing_x = _safe(geom, "btl_lhs")
    b_max_y = _safe(bearing_x, "max_y_bearing", L * 0.8)
    b_x = min(b_max_y, L * 0.85)
    ax.plot(b_x, 0, marker="v", color="blue", markersize=10, zorder=5)
    ax.text(b_x, 0.08, "Bearing", ha="center", va="bottom", fontsize=7, color="blue")

    # Dimension annotations
    margin = 0.15
    _add_dim_v(ax, L + margin, -D_curt, 0,
               f"{depth_curt:.0f} mm", offset=margin * 0.5, fontsize=7)
    _add_dim_v(ax, -pier_hw - margin, -D_pier, 0,
               f"{depth_pier:.0f} mm", offset=-margin * 3, fontsize=7)
    _add_dim_h(ax, margin, 0, L,
               f"Cantilever ≈ {L * 1000:.0f} mm", offset=0.04, fontsize=7)
    _add_dim_h(ax, -D_pier - pier_h - margin, -pier_hw, -pier_hw + pier_w_long,
               f"Pier {pier_w_long * 1000:.0f} mm", offset=-0.08, fontsize=7)

    # Section labels
    ax.text(0, -D_pier - 0.05, "A", fontsize=10, fontweight="bold",
            ha="center", va="top", color="red")
    ax.text(L * 0.6, -D_curt * 0.5 - D_pier * 0.5, "B", fontsize=10,
            fontweight="bold", ha="center", va="center", color="red")

    ax.set_xlim(-pier_hw - 0.6, L + 0.5)
    ax.set_ylim(-D_pier - pier_h - 0.5, 0.4)
    ax.set_aspect("equal")
    ax.set_title("Pier Cap — Tapered Elevation", fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=7, framealpha=0.8)
    ax.axis("off")

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# 2. Circular Cross-Section with Rebar (generic — used for pier & pile)
# ═══════════════════════════════════════════════════════════════════════════════

def _draw_circular_section(
    ax,
    diameter_mm: float,
    n_bars: int,
    bar_dia_mm: float,
    cover_mm: float,
    spiral_dia_mm: float = 0,
    spiral_spacing_mm: float = 0,
    title: str = "",
) -> None:
    """Draw a circular RC cross-section on the given axes."""

    R = diameter_mm / 2.0  # radius in mm

    # Concrete circle
    concrete = mpatches.Circle((0, 0), R, facecolor=_CONCRETE,
                               edgecolor="black", linewidth=1.5, zorder=1)
    ax.add_patch(concrete)

    # Spiral / hoop circle
    if spiral_dia_mm > 0:
        r_spiral = R - cover_mm - spiral_dia_mm / 2.0
        spiral = mpatches.Circle((0, 0), r_spiral, facecolor="none",
                                 edgecolor=_SPIRAL, linewidth=1.2,
                                 linestyle="--", zorder=2)
        ax.add_patch(spiral)

    # Rebar dots
    if n_bars > 0:
        r_bar_center = R - cover_mm - (spiral_dia_mm if spiral_dia_mm > 0 else 0) - bar_dia_mm / 2.0
        r_bar_center = max(r_bar_center, R * 0.3)  # safety
        bar_r_draw = max(bar_dia_mm / 2.0, R * 0.025)  # visual minimum
        for i in range(n_bars):
            angle = 2 * math.pi * i / n_bars
            bx = r_bar_center * math.cos(angle)
            by = r_bar_center * math.sin(angle)
            bar = mpatches.Circle((bx, by), bar_r_draw, facecolor=_REBAR,
                                  edgecolor="black", linewidth=0.4, zorder=3)
            ax.add_patch(bar)

    # Annotations
    pad = R * 0.15
    # Diameter
    ax.annotate(
        "", xy=(R, 0), xytext=(-R, 0),
        arrowprops=dict(arrowstyle="<->", color=_DIM, lw=0.8),
    )
    ax.text(0, -pad * 0.3, f"⌀ {diameter_mm:.0f} mm", ha="center", va="top",
            fontsize=8, color=_DIM, fontweight="bold")

    # Cover
    ax.annotate(
        "", xy=(R, R * 0.5), xytext=(R - cover_mm, R * 0.5),
        arrowprops=dict(arrowstyle="<->", color=_DIM, lw=0.6),
    )
    ax.text(R - cover_mm / 2, R * 0.5 + pad * 0.4, f"c={cover_mm:.0f}",
            ha="center", fontsize=6.5, color=_DIM)

    # Rebar info
    info_lines = [f"{n_bars} nos × ⌀{bar_dia_mm:.0f} mm"]
    if spiral_dia_mm > 0 and spiral_spacing_mm > 0:
        info_lines.append(f"Spiral ⌀{spiral_dia_mm:.0f} @ {spiral_spacing_mm:.0f} mm")
    ax.text(0, -(R + pad * 2), "\n".join(info_lines),
            ha="center", va="top", fontsize=7.5, color=_DIM,
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#ccc"))

    margin = R * 1.35
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin - pad * 2.5, margin)
    ax.set_aspect("equal")
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold")
    ax.axis("off")


def draw_pier_section(pier_result) -> plt.Figure:
    """Circular cross-section of the pier with rebar layout.

    Parameters
    ----------
    pier_result : PierDesignResult
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    _draw_circular_section(
        ax,
        diameter_mm=_safe(pier_result, "diameter", 2000),
        n_bars=int(_safe(pier_result, "n_bars", 20)),
        bar_dia_mm=_safe(pier_result, "bar_dia", 20),
        cover_mm=_safe(pier_result, "cover", 50),
        spiral_dia_mm=_safe(pier_result, "spiral_dia", 12),
        spiral_spacing_mm=_safe(pier_result, "spiral_spacing", 150),
        title="Pier — Cross Section",
    )
    fig.tight_layout()
    return fig


def draw_pile_section(pile_result) -> plt.Figure:
    """Circular cross-section of the pile with rebar layout.

    Parameters
    ----------
    pile_result : PileDesignResult
    """
    fig, ax = plt.subplots(figsize=(6, 6))
    _draw_circular_section(
        ax,
        diameter_mm=_safe(pile_result, "diameter", 1200),
        n_bars=int(_safe(pile_result, "n_bars", 16)),
        bar_dia_mm=_safe(pile_result, "bar_dia", 25),
        cover_mm=_safe(pile_result, "cover", 75),
        title="Pile — Cross Section",
    )
    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# 3. P-M Interaction Diagram (generic — used for pier and pile)
# ═══════════════════════════════════════════════════════════════════════════════

def draw_pm_diagram(
    interaction,
    P_applied: float,
    M_applied: float,
    title: str = "P-M Interaction Diagram",
) -> plt.Figure:
    """Plot P-M interaction curve with applied load point.

    Parameters
    ----------
    interaction : InteractionDiagram
        Has ``.points`` (list of InteractionPoint with .P and .M), plus
        ``P_max``, ``P_min``, ``M_max``, ``P_balanced``.
    P_applied : float
        Applied axial load (kN).
    M_applied : float
        Applied moment (kN.m).
    title : str
    """
    fig, ax = plt.subplots(figsize=(8, 5))

    if interaction is None:
        ax.text(0.5, 0.5, "No interaction data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        ax.set_title(title)
        return fig

    points = _safe(interaction, "points", [])
    if not points:
        ax.text(0.5, 0.5, "No interaction points", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        ax.set_title(title)
        return fig

    # Capacity curve
    Ms = [pt.M for pt in points]
    Ps = [pt.P for pt in points]
    ax.plot(Ms, Ps, "b-", linewidth=1.8, label="Capacity envelope", zorder=2)
    ax.fill_betweenx(Ps, 0, Ms, alpha=0.08, color="blue")

    # Key points
    P_max = _safe(interaction, "P_max", 0)
    P_min = _safe(interaction, "P_min", 0)
    M_max = _safe(interaction, "M_max", 0)
    P_bal = _safe(interaction, "P_balanced", 0)

    ax.plot(0, P_max, "bs", markersize=6, zorder=4)
    ax.text(M_max * 0.02 + 10, P_max, f"P_max = {P_max:.0f} kN",
            fontsize=7, va="bottom")

    ax.plot(M_max, P_bal, "b^", markersize=6, zorder=4)
    ax.text(M_max + M_max * 0.02, P_bal, f"Balanced ({M_max:.0f}, {P_bal:.0f})",
            fontsize=7, va="center")

    # Applied load point
    M_abs = abs(M_applied) if M_applied is not None else 0
    P_abs = P_applied if P_applied is not None else 0
    ax.plot(M_abs, P_abs, "ro", markersize=9, zorder=5, label="Applied (P, M)")

    # Line from origin to demand, extended to show utilisation
    ax.plot([0, M_abs], [0, P_abs], "r--", linewidth=1.0, alpha=0.6)

    ax.text(M_abs + M_max * 0.02, P_abs,
            f"({M_abs:.0f}, {P_abs:.0f})",
            fontsize=7, color="red", va="center")

    ax.set_xlabel("Moment M (kN.m)", fontsize=9)
    ax.set_ylabel("Axial Force P (kN)", fontsize=9)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color="gray", linewidth=0.5)
    ax.axvline(0, color="gray", linewidth=0.5)

    # Set limits with padding
    all_M = Ms + [M_abs]
    all_P = Ps + [P_abs]
    m_range = max(all_M) - min(all_M) if all_M else 1
    p_range = max(all_P) - min(all_P) if all_P else 1
    ax.set_xlim(-m_range * 0.05, max(all_M) * 1.15 + 1)
    ax.set_ylim(min(all_P) - p_range * 0.1, max(all_P) * 1.1 + 1)

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# 4. Pile Arrangement — Plan View
# ═══════════════════════════════════════════════════════════════════════════════

def draw_pile_arrangement(geom) -> plt.Figure:
    """Top-down plan view showing pile positions within the pile cap.

    Parameters
    ----------
    geom : GeometryResults
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    if geom is None:
        ax.text(0.5, 0.5, "No geometry data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        return fig

    # Pile cap dims (m)
    L = _safe(geom, "pilecap_length_long", 6)
    W = _safe(geom, "pilecap_width_trans", 10)

    # Pile cap rectangle centered on origin
    cap = mpatches.Rectangle(
        (-L / 2, -W / 2), L, W,
        linewidth=2, edgecolor="black", facecolor=_CONCRETE, zorder=1,
    )
    ax.add_patch(cap)

    # Pier outline
    pier_sec = _safe(geom, "pier_section")
    pier_type = _safe(pier_sec, "pier_type", "Circular")
    pier_wl = _safe(pier_sec, "width_long", 2.0)
    pier_wt = _safe(pier_sec, "width_trans", 2.0)

    if pier_type == "Circular":
        pier = mpatches.Circle(
            (0, 0), pier_wl / 2,
            facecolor="#b0b0b0", edgecolor="black", linewidth=1.5,
            linestyle="-", zorder=3,
        )
    else:
        pier = mpatches.Rectangle(
            (-pier_wl / 2, -pier_wt / 2), pier_wl, pier_wt,
            facecolor="#b0b0b0", edgecolor="black", linewidth=1.5, zorder=3,
        )
    ax.add_patch(pier)
    ax.text(0, 0, "PIER", ha="center", va="center", fontsize=9,
            fontweight="bold", color="white", zorder=4)

    # Pile circles
    pile_coords = _safe(geom, "pile_coordinates", [])
    pile_props = _safe(geom, "pile_props")
    pile_dia_m = _safe(pile_props, "diameter", 1.2)
    pile_r = pile_dia_m / 2.0

    for i, pc in enumerate(pile_coords):
        x = _safe(pc, "x", 0)
        y = _safe(pc, "y", 0)
        circle = mpatches.Circle(
            (x, y), pile_r,
            facecolor="white", edgecolor="black", linewidth=1.2, zorder=2,
        )
        ax.add_patch(circle)
        ax.text(x, y, f"P{i + 1}", ha="center", va="center", fontsize=7,
                fontweight="bold", zorder=3)

    # Dimension annotations
    pad = max(L, W) * 0.08

    # Pile cap dimensions
    _add_dim_h(ax, -W / 2 - pad, -L / 2, L / 2,
               f"{L * 1000:.0f} mm", offset=-pad * 0.6, fontsize=7)
    _add_dim_v(ax, L / 2 + pad, -W / 2, W / 2,
               f"{W * 1000:.0f} mm", offset=pad * 0.3, fontsize=7)

    # Pile spacing (if ≥ 2 piles exist)
    if len(pile_coords) >= 2:
        sp_long = _safe(geom, "pile_spacing_long", 0)
        sp_trans = _safe(geom, "pile_spacing_trans", 0)
        if sp_long > 0:
            # Find two piles in same column but adjacent rows
            first = pile_coords[0]
            for pc2 in pile_coords[1:]:
                if _safe(pc2, "col", 0) == _safe(first, "col", 0):
                    x0, y0 = _safe(first, "x", 0), _safe(first, "y", 0)
                    x1, y1 = _safe(pc2, "x", 0), _safe(pc2, "y", 0)
                    mid_y = (y0 + y1) / 2
                    ax.annotate(
                        "", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="<->", color="blue", lw=0.8),
                    )
                    ax.text((x0 + x1) / 2 - pad * 0.8, mid_y,
                            f"{sp_long * 1000:.0f}", fontsize=6.5, color="blue",
                            ha="center")
                    break

        if sp_trans > 0:
            first = pile_coords[0]
            for pc2 in pile_coords[1:]:
                if _safe(pc2, "row", 0) == _safe(first, "row", 0):
                    x0, y0 = _safe(first, "x", 0), _safe(first, "y", 0)
                    x1, y1 = _safe(pc2, "x", 0), _safe(pc2, "y", 0)
                    mid_x = (x0 + x1) / 2
                    ax.annotate(
                        "", xy=(x1, y1), xytext=(x0, y0),
                        arrowprops=dict(arrowstyle="<->", color="blue", lw=0.8),
                    )
                    ax.text(mid_x, (y0 + y1) / 2 + pad * 0.8,
                            f"{sp_trans * 1000:.0f}", fontsize=6.5, color="blue",
                            ha="center")
                    break

    # Axes labels
    ax.text(L / 2 + pad * 2.5, 0, "X (Long.)\n→ Traffic",
            ha="left", va="center", fontsize=7, color="#555")
    ax.text(0, W / 2 + pad * 2, "Y (Trans.)", ha="center", va="bottom",
            fontsize=7, color="#555")

    margin = max(L, W) * 0.3
    ax.set_xlim(-L / 2 - margin, L / 2 + margin)
    ax.set_ylim(-W / 2 - margin, W / 2 + margin)
    ax.set_aspect("equal")
    ax.set_title("Pile Arrangement — Plan View", fontsize=11, fontweight="bold")
    ax.axis("off")

    fig.tight_layout()
    return fig


# ═══════════════════════════════════════════════════════════════════════════════
# 5. Pile Cap Punching Shear Perimeters
# ═══════════════════════════════════════════════════════════════════════════════

def draw_pilecap_punching(pilecap_result, geom) -> plt.Figure:
    """Plan view of pile cap with punching shear control perimeters.

    Parameters
    ----------
    pilecap_result : PilecapDesignResult
    geom : GeometryResults
    """
    fig, ax = plt.subplots(figsize=(8, 8))

    if geom is None or pilecap_result is None:
        ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                ha="center", va="center", fontsize=12)
        return fig

    # Pile cap dims (from result, mm → m for drawing)
    L_mm = _safe(pilecap_result, "length_long", 6000)
    W_mm = _safe(pilecap_result, "width_trans", 10000)
    d_eff_mm = _safe(pilecap_result, "d_eff", 1400)
    L = L_mm / 1000.0
    W = W_mm / 1000.0
    d_eff = d_eff_mm / 1000.0

    # Pile cap rectangle
    cap = mpatches.Rectangle(
        (-L / 2, -W / 2), L, W,
        linewidth=2, edgecolor="black", facecolor=_CONCRETE, zorder=1,
    )
    ax.add_patch(cap)

    # Pier
    pier_sec = _safe(geom, "pier_section")
    pier_type = _safe(pier_sec, "pier_type", "Circular")
    pier_wl = _safe(pier_sec, "width_long", 2.0)
    pier_wt = _safe(pier_sec, "width_trans", 2.0)

    if pier_type == "Circular":
        pier_r = pier_wl / 2
        pier_patch = mpatches.Circle(
            (0, 0), pier_r,
            facecolor="#b0b0b0", edgecolor="black", linewidth=1.5, zorder=3,
        )
        ax.add_patch(pier_patch)

        # Punching perimeter at face (u_0)
        face_r = pier_r
        face_p = mpatches.Circle(
            (0, 0), face_r,
            facecolor="none", edgecolor=_PUNCH, linewidth=1.2,
            linestyle="-", zorder=4, label="Face perimeter",
        )
        ax.add_patch(face_p)

        # Punching perimeter at 2d (u_1)
        ctrl_r = pier_r + 2 * d_eff
        ctrl_p = mpatches.Circle(
            (0, 0), ctrl_r,
            facecolor="none", edgecolor=_PUNCH, linewidth=1.5,
            linestyle="--", zorder=4, label=f"Control at 2d ({2 * d_eff * 1000:.0f} mm)",
        )
        ax.add_patch(ctrl_p)

    else:
        # Rectangular pier
        pier_patch = mpatches.Rectangle(
            (-pier_wl / 2, -pier_wt / 2), pier_wl, pier_wt,
            facecolor="#b0b0b0", edgecolor="black", linewidth=1.5, zorder=3,
        )
        ax.add_patch(pier_patch)

        # Face perimeter
        face_rect = mpatches.Rectangle(
            (-pier_wl / 2, -pier_wt / 2), pier_wl, pier_wt,
            facecolor="none", edgecolor=_PUNCH, linewidth=1.2, zorder=4,
            label="Face perimeter",
        )
        ax.add_patch(face_rect)

        # Control perimeter at 2d (rounded rectangle approx)
        ext = 2 * d_eff
        ctrl_rect = mpatches.FancyBboxPatch(
            (-pier_wl / 2 - ext, -pier_wt / 2 - ext),
            pier_wl + 2 * ext, pier_wt + 2 * ext,
            boxstyle=mpatches.BoxStyle.Round(pad=ext),
            facecolor="none", edgecolor=_PUNCH, linewidth=1.5,
            linestyle="--", zorder=4,
            label=f"Control at 2d ({2 * d_eff * 1000:.0f} mm)",
        )
        ax.add_patch(ctrl_rect)

    ax.text(0, 0, "PIER", ha="center", va="center", fontsize=9,
            fontweight="bold", color="white", zorder=5)

    # Piles
    pile_coords = _safe(geom, "pile_coordinates", [])
    pile_props = _safe(geom, "pile_props")
    pile_dia_m = _safe(pile_props, "diameter", 1.2)
    pile_r = pile_dia_m / 2.0

    for i, pc in enumerate(pile_coords):
        x = _safe(pc, "x", 0)
        y = _safe(pc, "y", 0)
        circle = mpatches.Circle(
            (x, y), pile_r,
            facecolor="white", edgecolor="black", linewidth=1.0, zorder=2,
        )
        ax.add_patch(circle)

    # Punching perimeter around a corner pile (illustrative)
    if pile_coords:
        # Pick the pile with largest |x| + |y| (corner)
        corner = max(pile_coords, key=lambda p: abs(_safe(p, "x", 0)) + abs(_safe(p, "y", 0)))
        cx = _safe(corner, "x", 0)
        cy = _safe(corner, "y", 0)

        # Pile face perimeter
        pile_face = mpatches.Circle(
            (cx, cy), pile_r,
            facecolor="none", edgecolor="#e67e22", linewidth=1.2,
            linestyle="-", zorder=4,
        )
        ax.add_patch(pile_face)

        # Control perimeter at 2d from pile
        pile_ctrl_r = pile_r + 2 * d_eff
        pile_ctrl = mpatches.Circle(
            (cx, cy), pile_ctrl_r,
            facecolor="none", edgecolor="#e67e22", linewidth=1.2,
            linestyle=":", zorder=4, label="Pile control at 2d",
        )
        ax.add_patch(pile_ctrl)

    # 2d dimension annotation
    if pier_type == "Circular":
        dim_angle = math.pi / 4
        r_start = pier_r
        r_end = pier_r + 2 * d_eff
        ax.annotate(
            "", xy=(r_end * math.cos(dim_angle), r_end * math.sin(dim_angle)),
            xytext=(r_start * math.cos(dim_angle), r_start * math.sin(dim_angle)),
            arrowprops=dict(arrowstyle="<->", color=_PUNCH, lw=0.8),
        )
        ax.text((r_start + r_end) / 2 * math.cos(dim_angle) + 0.15,
                (r_start + r_end) / 2 * math.sin(dim_angle) + 0.15,
                f"2d = {2 * d_eff * 1000:.0f} mm",
                fontsize=7, color=_PUNCH)

    # Dimensions
    pad = max(L, W) * 0.06
    _add_dim_h(ax, -W / 2 - pad * 2, -L / 2, L / 2,
               f"{L_mm:.0f} mm", offset=-pad, fontsize=7)
    _add_dim_v(ax, L / 2 + pad * 2, -W / 2, W / 2,
               f"{W_mm:.0f} mm", offset=pad, fontsize=7)

    ax.legend(loc="upper left", fontsize=7, framealpha=0.8)
    margin = max(L, W) * 0.25
    ax.set_xlim(-L / 2 - margin, L / 2 + margin)
    ax.set_ylim(-W / 2 - margin, W / 2 + margin)
    ax.set_aspect("equal")
    ax.set_title("Pile Cap — Punching Shear Perimeters", fontsize=11, fontweight="bold")
    ax.axis("off")

    fig.tight_layout()
    return fig

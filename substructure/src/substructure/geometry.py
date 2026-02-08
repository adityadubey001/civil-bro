"""Geometry calculations for bridge substructure design.

Computes all derived dimensions from the validated input configuration:
levels, bearing positions, pier/piercap geometry, pile arrangement, pile cap
sizing, section properties, and fixity depth.  Handles both **Circular** and
**Rectangular** pier cross-sections.

Reference formulas extracted from the Excel *Inputs* sheet.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any


# ---------------------------------------------------------------------------
# Concrete elastic modulus (IRC 112 Table 6.5 short-term Ecm)
# ---------------------------------------------------------------------------
# Mean compressive strength: fcm = fck + 8  (IRC 112 Cl. 6.4.2.1)
# Ecm = 22 * (fcm / 10)^0.3  GPa  ->  multiply by 1000 for MPa

def _ecm_mpa(fck: float) -> float:
    """Return short-term elastic modulus Ecm in MPa per IRC 112."""
    fcm = fck + 8.0
    return 22_000.0 * (fcm / 10.0) ** 0.3


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class PileCoordinate:
    """Position of a single pile relative to the pile-cap centre."""
    row: int          # 1-based row index (longitudinal direction)
    col: int          # 1-based column index (transverse direction)
    x: float          # m  -- longitudinal offset from cap centre
    y: float          # m  -- transverse offset from cap centre


@dataclass
class BearingLevels:
    """Derived bearing / support levels for one side (LHS or RHS)."""
    max_y_bearing: float          # m  -- max transverse coordinate of bearings
    bearing_top_level: float      # m  -- top-of-bearing level (BTL)


@dataclass
class PierSection:
    """Cross-section properties for the pier shaft."""
    pier_type: str       # "Circular" or "Rectangular"
    area: float          # m2
    inertia_xx: float    # m4  -- about longitudinal axis
    inertia_yy: float    # m4  -- about transverse axis
    perimeter: float     # m
    width_long: float    # m  -- dimension in traffic direction
    width_trans: float   # m  -- dimension perpendicular to traffic


@dataclass
class PierCapSection:
    """Conservative rectangular section properties for the pier cap."""
    area: float       # m2
    perimeter: float  # m
    width: float      # m  -- width in longitudinal direction
    depth: float      # m  -- max depth


@dataclass
class PileProperties:
    """Derived properties for a single pile."""
    diameter: float             # m
    area: float                 # m2
    inertia: float              # m4  -- I_pile = pi/64 * D^4
    ecm: float                  # MPa -- pile concrete Ecm
    pile_length: float          # m   -- from bottom of pilecap to founding level
    pile_capacity_bottom: float # tonnes -- adjusted capacity at pile bottom
    stiffness_factor_T: float   # m   -- T = (Ecm * I / eta_b)^0.2
    fixity_depth_Lf: float      # m   -- Lf = fixity_ratio * T


@dataclass
class GeometryResults:
    """Aggregated results from all geometry calculations.

    Every field is populated by :func:`calculate_geometry`.
    """
    # -- Levels ---------------------------------------------------------
    pilecap_top_level: float              # m (ptl)
    pilecap_bottom_level: float           # m
    piercap_top_level: float              # m
    piercap_bottom_level: float           # m
    pier_height: float                    # m (piercap bottom to pilecap top)

    # -- Bearing / span -------------------------------------------------
    span_lhs: float                       # m  EJ-to-EJ left span
    span_rhs: float                       # m  EJ-to-EJ right span
    btl_lhs: BearingLevels
    btl_rhs: BearingLevels

    # -- Pier section ---------------------------------------------------
    pier_section: PierSection

    # -- Pier cap section -----------------------------------------------
    piercap_section: PierCapSection

    # -- Pile arrangement -----------------------------------------------
    pile_spacing_long: float              # m
    pile_spacing_trans: float             # m
    pile_coordinates: list[PileCoordinate]
    num_piles: int

    # -- Pile cap dimensions --------------------------------------------
    pilecap_length_long: float            # m  (in traffic direction)
    pilecap_width_trans: float            # m  (perpendicular to traffic)
    pilecap_thickness: float              # m

    # -- Pile derived ---------------------------------------------------
    pile_props: PileProperties

    # -- Raw bearing coordinates (echoed for downstream modules) --------
    bearing_coords_lhs: list[list[float]] = field(default_factory=list)
    bearing_coords_rhs: list[list[float]] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Internal calculation helpers
# ---------------------------------------------------------------------------

def _compute_bearing_levels(
    coords: list[list[float]],
    frl: float,
    camber_pct: float,
    wearing_coat_mm: float,
    depth_incl_slab: float,
) -> BearingLevels:
    """Derive Bearing Top Level from the bearing coordinates for one side.

    The bearing with the *largest absolute transverse offset* governs
    the camber drop.
    """
    max_y = max(abs(c[1]) for c in coords)
    camber_drop = max_y * camber_pct / 100.0
    wc_m = wearing_coat_mm / 1000.0
    btl = frl - camber_drop - wc_m - depth_incl_slab
    return BearingLevels(max_y_bearing=max_y, bearing_top_level=btl)


def _generate_pile_grid(
    piles_long: int,
    piles_trans: int,
    spacing_long: float,
    spacing_trans: float,
) -> list[PileCoordinate]:
    """Generate pile coordinates on a regular rectangular grid centred at
    the pilecap centroid.

    The grid is symmetric about both axes.

    Row indices run along the longitudinal (traffic) direction.
    Column indices run along the transverse direction.
    """
    coords: list[PileCoordinate] = []
    for row in range(1, piles_long + 1):
        x = (piles_long - 1) / 2.0 * spacing_long - (row - 1) * spacing_long
        for col in range(1, piles_trans + 1):
            y = (piles_trans - 1) / 2.0 * spacing_trans - (col - 1) * spacing_trans
            coords.append(PileCoordinate(row=row, col=col, x=x, y=y))
    return coords


def _pier_section_circular(diameter: float) -> PierSection:
    d = diameter
    area = math.pi * d ** 2 / 4.0
    inertia = math.pi * d ** 4 / 64.0
    perimeter = math.pi * d
    return PierSection(
        pier_type="Circular",
        area=area,
        inertia_xx=inertia,
        inertia_yy=inertia,
        perimeter=perimeter,
        width_long=d,
        width_trans=d,
    )


def _pier_section_rectangular(width_long: float, length_trans: float) -> PierSection:
    """Rectangular pier: *width_long* is in the traffic direction,
    *length_trans* is perpendicular."""
    area = width_long * length_trans
    inertia_xx = width_long * length_trans ** 3 / 12.0   # about long axis
    inertia_yy = length_trans * width_long ** 3 / 12.0   # about trans axis
    perimeter = 2.0 * (width_long + length_trans)
    return PierSection(
        pier_type="Rectangular",
        area=area,
        inertia_xx=inertia_xx,
        inertia_yy=inertia_yy,
        perimeter=perimeter,
        width_long=width_long,
        width_trans=length_trans,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def calculate_geometry(config: dict[str, Any]) -> GeometryResults:
    """Compute all derived geometric quantities from a validated config.

    Parameters
    ----------
    config:
        A dictionary returned by :func:`input_parser.parse_input`.

    Returns
    -------
    GeometryResults
        Dataclass with every derived dimension populated.
    """
    # Unpack input sections for readability.
    sup   = config["superstructure"]
    pier  = config["pier"]
    pcap  = config["pier_cap"]
    brg   = config["bearings"]
    lvl   = config["levels"]
    mat   = config["materials"]
    fnd   = config["foundation"]

    frl           = lvl["frl"]
    gl            = lvl["gl"]
    camber        = sup["camber_superelevation"]
    wc_mm         = sup["wearing_coat_thickness"]
    depth_slab    = sup["depth_incl_slab"]
    ped_depth     = sup["bearing_pedestal_depth"]
    offset_long   = pcap["bearing_offset_long"]

    # ----------------------------------------------------------------
    # 1. Levels
    # ----------------------------------------------------------------
    ptl = gl - fnd["pilecap_top_below_gl"] / 1000.0          # Pilecap Top Level
    pcb = ptl - fnd["pilecap_thickness"]                     # Pilecap Bottom Level

    # Spans EJ-to-EJ
    span_lhs = sup["left_span_cc_bearing"] + 2.0 * offset_long
    span_rhs = sup["right_span_cc_bearing"] + 2.0 * offset_long

    # Bearing top levels
    btl_lhs = _compute_bearing_levels(
        brg["coordinates_lhs"], frl, camber, wc_mm, depth_slab,
    )
    btl_rhs = _compute_bearing_levels(
        brg["coordinates_rhs"], frl, camber, wc_mm, depth_slab,
    )

    # Piercap levels
    piercap_top = min(btl_lhs.bearing_top_level,
                      btl_rhs.bearing_top_level) - ped_depth
    piercap_bottom = piercap_top - pcap["depth_max"]

    # Pier height (bottom of piercap to top of pilecap)
    pier_height = piercap_bottom - ptl

    # ----------------------------------------------------------------
    # 2. Pier section properties
    # ----------------------------------------------------------------
    if pier["type"] == "Circular":
        # For a circular pier the "diameter_bottom" is the main diameter.
        # If top != bottom we conservatively use the smaller for section
        # properties (most critical for capacity) but note both.
        d_design = min(pier["diameter_bottom"], pier["diameter_top"])
        pier_sec = _pier_section_circular(d_design)
    else:
        # Rectangular: diameter_bottom -> width_long, diameter_top -> length_trans
        pier_sec = _pier_section_rectangular(
            pier["diameter_bottom"], pier["diameter_top"]
        )

    # ----------------------------------------------------------------
    # 3. Pier cap section properties (conservative solid rectangle)
    # ----------------------------------------------------------------
    pcap_area = pcap["width_long"] * pcap["depth_max"]
    pcap_peri = 2.0 * (pcap["width_long"] + pcap["depth_max"])
    pcap_sec = PierCapSection(
        area=pcap_area,
        perimeter=pcap_peri,
        width=pcap["width_long"],
        depth=pcap["depth_max"],
    )

    # ----------------------------------------------------------------
    # 4. Pile arrangement
    # ----------------------------------------------------------------
    D = fnd["pile_diameter"]
    piles_long = fnd["piles_long"]
    piles_trans = fnd["piles_trans"]
    spacing_long = fnd["pile_spacing_factor"] * D
    spacing_trans = spacing_long  # typically identical

    pile_coords = _generate_pile_grid(
        piles_long, piles_trans, spacing_long, spacing_trans,
    )
    num_piles = len(pile_coords)

    # ----------------------------------------------------------------
    # 5. Pile cap dimensions (auto-sized from pile grid)
    # ----------------------------------------------------------------
    xs = [p.x for p in pile_coords]
    ys = [p.y for p in pile_coords]
    edge_clr_m = fnd["pilecap_edge_clearance"] / 1000.0

    pilecap_long = (max(xs) - min(xs)) + D + 2.0 * edge_clr_m
    pilecap_trans = (max(ys) - min(ys)) + D + 2.0 * edge_clr_m
    pilecap_thick = fnd["pilecap_thickness"]

    # ----------------------------------------------------------------
    # 6. Pile derived properties
    # ----------------------------------------------------------------
    pile_area = math.pi * D ** 2 / 4.0
    I_pile = math.pi * D ** 4 / 64.0
    ecm_pile = _ecm_mpa(mat["pile_pilecap"]["fck"])

    # Pile length: from bottom of pilecap to founding level
    pile_length = (
        fnd["total_pile_depth"]
        - fnd["pilecap_thickness"]
        - fnd["pilecap_top_below_gl"] / 1000.0
    )

    # Pile capacity at bottom (self-weight added, converted to tonnes)
    # concrete_density is in kN/m3; divide by 10 to get tonnes/m
    concrete_density_knm3 = mat["concrete_density"]
    selfweight_kn = pile_area * concrete_density_knm3 * pile_length  # kN
    selfweight_tonnes = selfweight_kn / 10.0  # 1 tonne-force ~ 10 kN
    pile_capacity_bottom = fnd["pile_capacity_top"] + selfweight_tonnes

    # Fixity depth
    # subgrade_modulus is in MN/m3; convert to kN/m3 for consistency
    # eta_b = k_s * D  (kN/m2 per m depth -> kN/m3 effectively)
    # However the formula T = (E*I / eta_b)^0.2 uses consistent units.
    # Here: Ecm in kPa = ecm_mpa * 1000, I in m4, eta_b in kN/m3.
    eta_b_knm3 = fnd["subgrade_modulus"] * 1000.0  # MN/m3 -> kN/m3
    ecm_kpa = ecm_pile * 1000.0                    # MPa -> kPa (= kN/m2)
    # T = (Ecm_kPa * I_pile / eta_b_kN/m3)^0.2  -> units: (kN/m2 * m4 / (kN/m3))^0.2 = m^0.2*5 = m
    T = (ecm_kpa * I_pile / eta_b_knm3) ** 0.2
    Lf = fnd["fixity_ratio"] * T

    pile_props = PileProperties(
        diameter=D,
        area=pile_area,
        inertia=I_pile,
        ecm=ecm_pile,
        pile_length=pile_length,
        pile_capacity_bottom=pile_capacity_bottom,
        stiffness_factor_T=T,
        fixity_depth_Lf=Lf,
    )

    # ----------------------------------------------------------------
    # Assemble results
    # ----------------------------------------------------------------
    return GeometryResults(
        # Levels
        pilecap_top_level=ptl,
        pilecap_bottom_level=pcb,
        piercap_top_level=piercap_top,
        piercap_bottom_level=piercap_bottom,
        pier_height=pier_height,
        # Spans
        span_lhs=span_lhs,
        span_rhs=span_rhs,
        btl_lhs=btl_lhs,
        btl_rhs=btl_rhs,
        # Pier
        pier_section=pier_sec,
        # Pier cap
        piercap_section=pcap_sec,
        # Piles
        pile_spacing_long=spacing_long,
        pile_spacing_trans=spacing_trans,
        pile_coordinates=pile_coords,
        num_piles=num_piles,
        # Pile cap
        pilecap_length_long=pilecap_long,
        pilecap_width_trans=pilecap_trans,
        pilecap_thickness=pilecap_thick,
        # Pile derived
        pile_props=pile_props,
        # Echoed bearing coords
        bearing_coords_lhs=brg["coordinates_lhs"],
        bearing_coords_rhs=brg["coordinates_rhs"],
    )

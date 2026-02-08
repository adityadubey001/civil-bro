"""Load calculations for bridge substructure design per IRC 6:2017.

Computes all design loads as 6-DOF force vectors (P, FL, FT, ML, MT, Mv) at
the pier base level (top of pile cap).  Results are structured for direct
consumption by :mod:`load_combinations`.

Key references
--------------
* IRC 6:2017  -- Standard Specifications and Code of Practice for Road
  Bridges, Section II: Loads and Load Combinations
* IRC SP-114:2018 -- Guidelines for Seismic Design of Road Bridges

Load types computed
-------------------
1. Dead Load (DL) -- pier, pier cap, superstructure self-weight
2. Superimposed Dead Load (SIDL) -- crash barriers, noise barriers, wearing coat
3. Live Load (LL) -- Class A, Class 70R with impact, lane reduction, congestion
4. Braking -- longitudinal force from vehicles (IRC 6 Cl.206)
5. Centrifugal -- transverse force on curved bridges (IRC 6 Cl.207)
6. Wind -- on superstructure, live load, and substructure (IRC 6 Cl.209)
7. Seismic -- elastic seismic coefficient method (IRC 6 Cl.211 / IRC SP-114)
8. Temperature -- uniform thermal expansion/contraction (IRC 6 Cl.215)

Validation reference values (from Excel workbook)
--------------------------------------------------
- Piercap weight: 951.625 kN -> with 10%: 1046.79 kN
- Pier weight: 575.81 kN
- DL per bearing: 1064.72 kN, total per span: 8517.75 kN
- Crash barrier: 8.08 kN/m each side
- WC: 10.89 kN/m
- Impact factor (RC, 60m): 0.088
- Braking (3xClassA): 312.8 kN
- Braking (1x70R+1xClassA): 485.8 kN
- Lane reduction factor (3 lanes): 0.9
- Congestion factor (60m span): 1.6
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from .geometry import GeometryResults
from .utils import load_irc_tables


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_G: float = 9.81             # m/s^2 -- standard gravity
_RHO_CONCRETE: float = 25.0  # kN/m^3 -- default, overridden by config

# IRC Class A single train axle loads (kN) and spacings (m) between axles
_CLASS_A_AXLES: list[float] = [27, 27, 114, 114, 68, 68, 68, 68]
_CLASS_A_SPACINGS: list[float] = [0, 1.1, 3.2, 1.2, 4.3, 3.0, 3.0, 3.0]
_CLASS_A_TOTAL: float = 554.0  # kN

# IRC 70R Wheeled axle loads (kN) and spacings (m) between axles
_70R_AXLES: list[float] = [80, 120, 120, 170, 170, 170, 170]
_70R_SPACINGS: list[float] = [0, 3.96, 1.52, 2.13, 1.37, 3.05, 1.37]
_70R_TOTAL: float = 1000.0  # kN

# Lane reduction factors (IRC 6 Table 5)
_LANE_REDUCTION: dict[int, float] = {
    1: 1.0, 2: 1.0, 3: 0.9, 4: 0.75, 5: 0.6, 6: 0.5,
}

# Congestion factor breakpoints (IRC 6 Cl.204.4)
# Interpolation table: span (m) -> congestion factor
_CONGESTION_SPANS: list[float] = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
_CONGESTION_FACTORS: list[float] = [1.15, 1.3, 1.3, 1.3, 1.5, 1.6, 1.7, 1.7, 1.7, 1.7]

# Approximate girder cross-section areas by superstructure type (m^2)
_GIRDER_AREA_APPROX: dict[str, float] = {
    "U Girder": 1.233,
    "Box Girder": 2.5,
    "I Girder": 0.6,
    "Slab": 0.0,  # slab is computed separately
}

# Default design speed for centrifugal force (km/h)
_DESIGN_SPEED_KMPH: float = 50.0

# Wind gust factor per IRC 6:2017 Cl.209.3.3.
# Pz is computed from the hourly mean wind speed (IRC 6:2017 Table 5 converts
# the IS 875 basic 3-second gust speed to hourly mean).  The gust factor G=2.0
# is then applied to account for gust effects in the static approach.
_GUST_FACTOR: float = 2.0

# Lift coefficient for deck plan area (IRC 6 Cl.209.3.4)
_CL_DECK: float = 0.75

# IRC 5 crash barrier profile: approximately 0.323 m^2 cross section
_CB_AREA_IRC5: float = 0.323  # m^2


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ForceVector:
    """6-DOF force/moment vector at a point.

    Sign conventions (viewed from traffic approach direction):
    - P   : vertical, positive downward (kN)
    - FL  : horizontal along traffic direction, positive in direction of
            travel (kN)
    - FT  : horizontal perpendicular to traffic, positive to the right (kN)
    - ML  : moment about the longitudinal (traffic) axis (kN.m)
    - MT  : moment about the transverse axis (kN.m)
    - Mv  : vertical torsion about the vertical axis (kN.m)
    """

    P: float = 0.0
    FL: float = 0.0
    FT: float = 0.0
    ML: float = 0.0
    MT: float = 0.0
    Mv: float = 0.0

    # -- Aliases for backward compatibility with load_combinations -----------
    # The combinations module references .HL and .HT on ForceVector instances.

    @property
    def HL(self) -> float:
        """Alias for FL (horizontal longitudinal)."""
        return self.FL

    @property
    def HT(self) -> float:
        """Alias for FT (horizontal transverse)."""
        return self.FT

    # -- Arithmetic ----------------------------------------------------------

    def __add__(self, other: ForceVector) -> ForceVector:
        if not isinstance(other, ForceVector):
            return NotImplemented
        return ForceVector(
            P=self.P + other.P,
            FL=self.FL + other.FL,
            FT=self.FT + other.FT,
            ML=self.ML + other.ML,
            MT=self.MT + other.MT,
            Mv=self.Mv + other.Mv,
        )

    def __sub__(self, other: ForceVector) -> ForceVector:
        if not isinstance(other, ForceVector):
            return NotImplemented
        return ForceVector(
            P=self.P - other.P,
            FL=self.FL - other.FL,
            FT=self.FT - other.FT,
            ML=self.ML - other.ML,
            MT=self.MT - other.MT,
            Mv=self.Mv - other.Mv,
        )

    def __mul__(self, factor: float) -> ForceVector:
        return ForceVector(
            P=self.P * factor,
            FL=self.FL * factor,
            FT=self.FT * factor,
            ML=self.ML * factor,
            MT=self.MT * factor,
            Mv=self.Mv * factor,
        )

    def __rmul__(self, factor: float) -> ForceVector:
        return self.__mul__(factor)

    def __neg__(self) -> ForceVector:
        return self.__mul__(-1.0)


@dataclass
class LLCase:
    """A single live-load case with per-side and combined forces."""
    name: str
    lhs: ForceVector              # per-side (left span)
    rhs: ForceVector              # per-side (right span)
    force_at_pier_base: ForceVector  # combined = lhs + rhs


@dataclass
class LoadResults:
    """All computed loads for the substructure.

    All force vectors are resolved at the **pier base** level (top of pile
    cap) unless explicitly noted.
    """

    # -- Dead loads ---------------------------------------------------------
    dl_piercap: ForceVector         # Self-weight of pier cap (incl 10% for pedestals)
    dl_pier: ForceVector            # Self-weight of pier
    dl_super_lhs: ForceVector       # Superstructure DL, left span total
    dl_super_rhs: ForceVector       # Superstructure DL, right span total

    # -- Superimposed dead loads -------------------------------------------
    sidl_lhs: ForceVector           # SIDL (crash barrier etc), left span total
    sidl_rhs: ForceVector           # SIDL (crash barrier etc), right span total
    wc_lhs: ForceVector             # Wearing coat, left span total
    wc_rhs: ForceVector             # Wearing coat, right span total

    # -- Live load cases ---------------------------------------------------
    ll_cases: list[LLCase]          # [Max P, Max ML, Max MT, Congestion]

    # -- Braking (one per LL case) -----------------------------------------
    braking_cases: list[dict[str, ForceVector]]

    # -- Centrifugal (one per LL case) -------------------------------------
    centrifugal_cases: list[dict[str, ForceVector]]

    # -- Wind loads --------------------------------------------------------
    wind_on_super_lhs: ForceVector   # wind on superstructure, left span
    wind_on_super_rhs: ForceVector   # wind on superstructure, right span
    wind_on_ll_lhs: ForceVector      # wind on live load, left span
    wind_on_ll_rhs: ForceVector      # wind on live load, right span
    wind_on_pier: ForceVector        # wind on pier
    wind_on_piercap: ForceVector     # wind on pier cap
    wind_vertical_lhs: ForceVector   # vertical wind uplift, left
    wind_vertical_rhs: ForceVector   # vertical wind uplift, right

    # -- Seismic loads -----------------------------------------------------
    seismic_super_long: ForceVector   # seismic on superstructure, longitudinal
    seismic_super_trans: ForceVector  # seismic on superstructure, transverse
    seismic_pier_long: ForceVector    # seismic on pier, longitudinal
    seismic_pier_trans: ForceVector   # seismic on pier, transverse
    seismic_piercap_long: ForceVector # seismic on pier cap, longitudinal
    seismic_piercap_trans: ForceVector  # seismic on pier cap, transverse
    seismic_hydrodynamic: ForceVector # hydrodynamic, if applicable
    Ah_long: float                    # seismic coefficient longitudinal
    Ah_trans: float                   # seismic coefficient transverse

    # -- Temperature loads -------------------------------------------------
    temp_long: ForceVector            # temperature force in longitudinal direction

    # -- Summary forces at pier base ---------------------------------------
    total_dl_at_pier_base: ForceVector
    total_sidl_at_pier_base: ForceVector
    total_wc_at_pier_base: ForceVector

    # -- Per-bearing reactions (keyed by load type) ------------------------
    bearing_reactions: dict[str, dict[str, list[float]]] = field(
        default_factory=dict,
    )

    # -- Convenience aliases used by downstream modules ----------------------
    @property
    def wind_on_super(self) -> ForceVector:
        """Total wind on superstructure (LHS + RHS)."""
        return self.wind_on_super_lhs + self.wind_on_super_rhs

    @property
    def wind_on_ll(self) -> ForceVector:
        """Total wind on live load (LHS + RHS)."""
        return self.wind_on_ll_lhs + self.wind_on_ll_rhs

    @property
    def wind_on_sub(self) -> ForceVector:
        """Total wind on substructure (pier + piercap)."""
        return self.wind_on_pier + self.wind_on_piercap

    @property
    def eq_super_long(self) -> ForceVector:
        """Alias for seismic_super_long."""
        return self.seismic_super_long

    @property
    def eq_super_trans(self) -> ForceVector:
        """Alias for seismic_super_trans."""
        return self.seismic_super_trans

    @property
    def eq_ll(self) -> ForceVector:
        """Seismic on live load -- combined longitudinal + transverse."""
        return ForceVector()

    @property
    def eq_sub(self) -> ForceVector:
        """Total seismic on substructure."""
        return (
            self.seismic_pier_long + self.seismic_pier_trans
            + self.seismic_piercap_long + self.seismic_piercap_trans
        )

    @property
    def temp_rise(self) -> ForceVector:
        """Alias for temperature longitudinal (rise case)."""
        return self.temp_long

    @property
    def temp_fall(self) -> ForceVector:
        """Alias for temperature longitudinal (fall case -- negative)."""
        return -self.temp_long


# ---------------------------------------------------------------------------
# Interpolation helpers
# ---------------------------------------------------------------------------

def _interp(x: float, xs: list[float], ys: list[float]) -> float:
    """Linear interpolation with clamping at boundaries."""
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    return float(np.interp(x, xs, ys))


def _get_height_factor(height_m: float, irc: dict[str, Any]) -> float:
    """Interpolate the wind height factor k2 from IRC 6 Table."""
    hf_table = irc["wind"]["height_factor"]
    heights = sorted(hf_table.keys())
    factors = [hf_table[h] for h in heights]
    return _interp(height_m, [float(h) for h in heights], factors)


def _congestion_factor(span_m: float) -> float:
    """Congestion factor per IRC 6 Cl.204.4, interpolated by span."""
    return _interp(span_m, _CONGESTION_SPANS, _CONGESTION_FACTORS)


def _lane_reduction_factor(num_lanes: int) -> float:
    """Lane reduction factor per IRC 6 Table 5."""
    if num_lanes <= 0:
        return 1.0
    if num_lanes >= 6:
        return _LANE_REDUCTION[6]
    return _LANE_REDUCTION.get(num_lanes, 0.5)


# ---------------------------------------------------------------------------
# Number of lanes per IRC 6 Table 5
# ---------------------------------------------------------------------------

def _num_lanes(clear_carriageway_width: float) -> int:
    """Determine number of design lanes from clear carriageway width.

    IRC 6 Table 5 lane widths: 3.5 m standard lane.
    Carriageway < 5.3 m : 1 lane
    5.3 <= cw < 9.6 m   : 2 lanes
    9.6 <= cw < 13.1 m  : 3 lanes
    13.1 <= cw < 16.6 m : 4 lanes
    Beyond : floor(cw / 3.5)
    """
    cw = clear_carriageway_width
    if cw < 5.3:
        return 1
    if cw < 9.6:
        return 2
    if cw < 13.1:
        return 3
    if cw < 16.6:
        return 4
    return max(1, int(cw / 3.5))


# ---------------------------------------------------------------------------
# IRC 6 Response Spectrum (Sa/g) per soil type
# ---------------------------------------------------------------------------

def _sa_over_g(T: float, soil_type: str) -> float:
    """Spectral acceleration coefficient Sa/g per IRC SP-114 / IRC 6 Fig.18.

    Parameters
    ----------
    T : float
        Natural period of vibration (seconds).
    soil_type : str
        One of "hard", "medium", "soft".
    """
    if T <= 0.0:
        return 2.5

    # Common initial ramp for all soil types
    if T <= 0.10:
        return 1.0 + 15.0 * T

    if soil_type == "hard":
        if T <= 0.40:
            return 2.5
        if T <= 4.0:
            return 1.0 / T
        return 0.25

    elif soil_type == "medium":
        if T <= 0.55:
            return 2.5
        if T <= 4.0:
            return 1.36 / T
        return 0.34

    else:  # soft
        if T <= 0.67:
            return 2.5
        if T <= 4.0:
            return 1.67 / T
        return 0.42


# ---------------------------------------------------------------------------
# Impact factor per IRC 6 Cl.204
# ---------------------------------------------------------------------------

def _impact_factor_rc(span_m: float, vehicle: str) -> float:
    """Impact factor for reinforced/prestressed concrete bridges.

    IF = 4.5 / (6 + L), capped per IRC 6 Table.
    """
    L = span_m
    IF = 4.5 / (6.0 + L)
    if vehicle in ("class_70R_wheeled", "70R"):
        return min(IF, 0.25)
    else:  # class_A
        return min(IF, 0.545)


# ---------------------------------------------------------------------------
# Influence-line reaction for a train of axle loads on a simply-supported span
# ---------------------------------------------------------------------------

def _max_reaction_simply_supported(
    span: float,
    axle_loads: list[float],
    spacings: list[float],
) -> float:
    """Maximum end reaction at the near support for a simply-supported span.

    Places the axle train at the position that maximises the left reaction
    by sliding the train across the span and summing contributions.

    The influence line ordinate for left reaction at position x from the
    left support is: eta(x) = (span - x) / span.

    The train is slid with the first axle from position 0 to span, and
    for each position all axles on span contribute.
    """
    n = len(axle_loads)
    # Compute absolute positions of each axle from the first axle
    positions = [0.0]
    for i in range(1, n):
        positions.append(positions[i - 1] + spacings[i])

    best_reaction = 0.0
    # Slide the first axle across the span in small increments
    num_steps = max(500, int(span * 20))
    for step in range(num_steps + 1):
        x_first = step * span / num_steps
        reaction = 0.0
        for k in range(n):
            x_axle = x_first + positions[k]
            if 0.0 <= x_axle <= span:
                # Reaction at left support = P * (span - x) / span
                reaction += axle_loads[k] * (span - x_axle) / span
        best_reaction = max(best_reaction, reaction)

    return best_reaction


# ---------------------------------------------------------------------------
# Bearing eccentricity distribution
# ---------------------------------------------------------------------------

def _distribute_to_bearings(
    total_P: float,
    total_ML: float,
    coords: list[list[float]],
) -> list[float]:
    """Distribute vertical force and transverse moment to bearings.

    Uses the elastic distribution formula:
        R_i = P/n + ML * y_i / sum(y_j^2)

    where y_i is the transverse coordinate of bearing i, ML is the moment
    about the longitudinal axis.
    """
    n = len(coords)
    if n == 0:
        return []
    ys = [c[1] for c in coords]
    sum_y2 = sum(y ** 2 for y in ys)
    reactions: list[float] = []
    for y in ys:
        r = total_P / n
        if sum_y2 > 1e-9:
            r += total_ML * y / sum_y2
        reactions.append(r)
    return reactions


# ---------------------------------------------------------------------------
# Ecm computation (IRC 112)
# ---------------------------------------------------------------------------

def _ecm_mpa(fck: float) -> float:
    """Short-term elastic modulus Ecm in MPa (IRC 112)."""
    fcm = fck + 8.0
    return 22_000.0 * (fcm / 10.0) ** 0.3


# ---------------------------------------------------------------------------
# 1. Dead Load
# ---------------------------------------------------------------------------

def _compute_dead_loads(
    config: dict[str, Any],
    geom: GeometryResults,
) -> tuple[ForceVector, ForceVector, ForceVector, ForceVector, dict]:
    """Compute dead-load force vectors at pier base.

    Returns (dl_piercap, dl_pier, dl_super_lhs, dl_super_rhs, bearing_dl).

    Pier cap volume formula from Excel workbook:
      Volume = (length_trans * width_long * avg(depth_max, depth_min))
             + (pier_dia + 2*0.15)^2 * (depth_max - avg(depth_max, depth_min)) / 2
      This accounts for the trapezoidal shape with full depth at pier location.

    A simpler but equivalent form:
      V_rect = width_long * avg_depth * length_trans
      V_haunch = width_long * (depth_max - depth_min)/2 * pier_dia
      V_total = V_rect + V_haunch + pedestal_area_adder
    Then multiply by 1.1 for pedestals.
    """
    sup = config["superstructure"]
    pcap = config["pier_cap"]
    pier_cfg = config["pier"]
    mat = config["materials"]
    density = float(mat.get("concrete_density", _RHO_CONCRETE))

    pier_base_level = geom.pilecap_top_level

    # -- Pier cap self-weight --
    # Using the Excel formula approach:
    # Component 1: Rectangular portion at average depth over full transverse
    # length excluding pier width
    pcap_w = float(pcap["width_long"])       # width in traffic direction
    pcap_lt = float(pcap["length_trans"])     # length transverse
    depth_max = float(pcap["depth_max"])
    depth_min = float(pcap["depth_min"])
    avg_depth = (depth_max + depth_min) / 2.0

    # Pier dimension at cap level (bottom of cap)
    if pier_cfg["type"] == "Circular":
        pier_dim = float(pier_cfg["diameter_bottom"])
    else:
        pier_dim = float(pier_cfg["diameter_top"])  # transverse dimension

    # Volume of pier cap (Excel-compatible formula):
    # The pier cap is trapezoidal in elevation:
    # - Full depth (depth_max) over the pier width
    # - Tapers to depth_min at the cantilever tips
    # Volume = cantilever portions + full-depth pier-width portion
    # Cantilever each side length = (pcap_lt - pier_dim) / 2
    # Cantilever volume = 2 * pcap_w * cant_length * avg_depth
    # Pier portion volume = pcap_w * pier_dim * depth_max
    cant_length = (pcap_lt - pier_dim) / 2.0
    if cant_length < 0:
        cant_length = 0.0

    vol_cant = 2.0 * pcap_w * cant_length * avg_depth
    vol_pier_portion = pcap_w * pier_dim * depth_max

    # Additional allowance: pedestal zone -- pier area with fillet
    # Excel adds (pier_dim + 2*0.15)^2 term for pedestal/haunch zone
    pedestal_adder = (pier_dim + 2.0 * 0.15) ** 2
    vol_pcap = vol_cant + vol_pier_portion + pedestal_adder

    wt_pcap_raw = vol_pcap * density
    # Add 10% for bearing pedestals, seating ledges, haunches
    wt_pcap = wt_pcap_raw * 1.1

    dl_piercap = ForceVector(P=wt_pcap)

    # -- Pier self-weight --
    area = geom.pier_section.area     # m^2
    h_pier = geom.pier_height         # m
    wt_pier = area * h_pier * density
    # CG at mid-height; concentric with pier axis, no moment
    dl_pier = ForceVector(P=wt_pier)

    # -- Superstructure DL per bearing --
    n_brg_lhs = len(geom.bearing_coords_lhs)
    n_brg_rhs = len(geom.bearing_coords_rhs)

    dl_per_brg_lhs = _super_dl_per_bearing(
        sup, density, geom.span_lhs, n_brg_lhs, "lhs",
    )
    dl_per_brg_rhs = _super_dl_per_bearing(
        sup, density, geom.span_rhs, n_brg_rhs, "rhs",
    )

    # Total superstructure DL per side
    total_P_lhs = dl_per_brg_lhs * n_brg_lhs
    total_P_rhs = dl_per_brg_rhs * n_brg_rhs

    # For DL, load is uniformly distributed so per-bearing is equal: no
    # eccentricity moment.
    dl_super_lhs = ForceVector(P=total_P_lhs)
    dl_super_rhs = ForceVector(P=total_P_rhs)

    # Per-bearing DL reactions
    bearing_dl = {
        "lhs": [dl_per_brg_lhs] * n_brg_lhs,
        "rhs": [dl_per_brg_rhs] * n_brg_rhs,
    }

    return dl_piercap, dl_pier, dl_super_lhs, dl_super_rhs, bearing_dl


def _super_dl_per_bearing(
    sup: dict[str, Any],
    density: float,
    span: float,
    n_bearings: int,
    side: str,
) -> float:
    """DL reaction per bearing for one side of the pier.

    If the user provided a non-zero dl_reaction, use that directly.
    Otherwise auto-calculate from geometry:
      DL = (girder_area * span * density * n_girders
            + slab_thickness * deck_width * span * density
            + diaphragm_allowance) / 2 / n_bearings
    """
    key = f"dl_reaction_{side}"
    user_val = float(sup.get(key, 0.0))
    if user_val > 0.0:
        return user_val

    # Auto-calculate
    deck_type = sup["type"]
    n_girders = int(sup["num_girders"])
    deck_width = float(sup["deck_width"])
    slab_t = float(sup["deck_slab_thickness"])

    # Girder weight per full span
    girder_area = _GIRDER_AREA_APPROX.get(deck_type, 1.0)
    girder_weight = girder_area * span * density * n_girders

    # Deck slab weight per full span
    slab_weight = slab_t * deck_width * span * density

    # Cross girder / diaphragm allowance: ~2% of girder + slab
    diaphragm_weight = 0.02 * (girder_weight + slab_weight)

    total_span_weight = girder_weight + slab_weight + diaphragm_weight

    # Each pier supports half the span from each side; DL per bearing
    half_span_weight = total_span_weight / 2.0
    if n_bearings <= 0:
        return half_span_weight
    return half_span_weight / n_bearings


# ---------------------------------------------------------------------------
# 2. Superimposed Dead Load (SIDL)
# ---------------------------------------------------------------------------

def _compute_sidl(
    config: dict[str, Any],
    geom: GeometryResults,
) -> tuple[ForceVector, ForceVector, ForceVector, ForceVector, dict]:
    """Compute SIDL forces at pier base.

    Returns (sidl_lhs, sidl_rhs, wc_lhs, wc_rhs, bearing_sidl).

    Crash barrier: IRC-5 profile, ~0.323 m^2 cross section * 25 kN/m^3
        = 8.075 kN/m per barrier (rounded to 8.08 kN/m in Excel)
    Noise barrier: ~1 kN/m per side (lightweight steel/acrylic panels)
    Wearing coat: thickness_mm/1000 * 22 * carriageway_width (kN/m)
        Example: 45/1000 * 22 * 11 = 10.89 kN/m
    """
    sup = config["superstructure"]
    mat = config["materials"]
    density = float(mat.get("concrete_density", _RHO_CONCRETE))

    n_brg_lhs = len(geom.bearing_coords_lhs)
    n_brg_rhs = len(geom.bearing_coords_rhs)

    # Crash barrier weight (kN/m): IRC-5 profile
    # Area ~0.323 m^2 * 25 kN/m^3 = 8.075 kN/m
    cb_line_load = _CB_AREA_IRC5 * density   # ~8.075 kN/m per barrier

    # Median barrier (if present)
    mb_width = float(sup.get("median_barrier_width", 0.0))
    mb_line_load = 0.0
    if mb_width > 0:
        # Median barrier: similar profile, scaled by width ratio
        cb_width = float(sup.get("crash_barrier_width", 0.5))
        if cb_width > 0:
            mb_line_load = cb_line_load * (mb_width / cb_width)
        else:
            mb_line_load = _CB_AREA_IRC5 * density

    # Noise barrier: approximate as lightweight panel, ~1 kN/m per side
    nb_height = float(sup.get("noise_barrier_height", 0.0))
    nb_line_load = 1.0 if nb_height > 0.0 else 0.0

    # Total SIDL line load (kN/m):
    # 2 crash barriers (one on each edge) + median + 2 noise barriers
    sidl_line_load = 2.0 * cb_line_load + mb_line_load + 2.0 * nb_line_load

    # Wearing coat (WC) line load
    deck_width = float(sup["deck_width"])
    wc_mm = float(sup["wearing_coat_thickness"])
    cb_width = float(sup.get("crash_barrier_width", 0.5))
    clear_cw = deck_width - 2.0 * cb_width  # clear carriageway width
    wc_density = 22.0  # kN/m^3 for bituminous wearing coat
    wc_line_load = wc_mm / 1000.0 * wc_density * clear_cw  # kN/m

    # SIDL and WC per side of pier
    def _sidl_for_side(
        span: float,
        n_brg: int,
        side: str,
        coords: list[list[float]],
    ) -> tuple[ForceVector, ForceVector, list[float], list[float]]:
        """Compute SIDL and WC force vectors for one side."""
        # Check user override
        user_sidl = float(sup.get(f"sidl_reaction_{side}", 0.0))

        if user_sidl > 0.0:
            total_sidl = user_sidl * n_brg
        else:
            # Half-span contribution to this pier
            total_sidl = sidl_line_load * span / 2.0

        total_wc = wc_line_load * span / 2.0

        # Symmetric barriers -> no net transverse eccentricity for DL-type SIDL
        sidl_per_brg = _distribute_to_bearings(total_sidl, 0.0, coords)
        wc_per_brg = _distribute_to_bearings(total_wc, 0.0, coords)

        fv_sidl = ForceVector(P=total_sidl)
        fv_wc = ForceVector(P=total_wc)

        return fv_sidl, fv_wc, sidl_per_brg, wc_per_brg

    sidl_lhs, wc_lhs, sidl_brg_lhs, wc_brg_lhs = _sidl_for_side(
        geom.span_lhs, n_brg_lhs, "lhs", geom.bearing_coords_lhs,
    )
    sidl_rhs, wc_rhs, sidl_brg_rhs, wc_brg_rhs = _sidl_for_side(
        geom.span_rhs, n_brg_rhs, "rhs", geom.bearing_coords_rhs,
    )

    bearing_sidl = {
        "sidl_lhs": sidl_brg_lhs,
        "sidl_rhs": sidl_brg_rhs,
        "wc_lhs": wc_brg_lhs,
        "wc_rhs": wc_brg_rhs,
    }

    return sidl_lhs, sidl_rhs, wc_lhs, wc_rhs, bearing_sidl


# ---------------------------------------------------------------------------
# 3. Live Load
# ---------------------------------------------------------------------------

def _compute_live_load(
    config: dict[str, Any],
    geom: GeometryResults,
) -> list[LLCase]:
    """Compute live-load cases (Max P, Max ML, Max MT, Congestion).

    Uses simplified influence-line approach.  Returns a list of LLCase objects.

    Each ll_case entry is a dict with 'lhs' and 'rhs' ForceVector.

    The user may provide manual LL reactions via the ``live_load`` section
    of the config YAML.  If provided, those are used directly.
    """
    sup = config["superstructure"]
    deck_width = float(sup["deck_width"])
    cb_width = float(sup.get("crash_barrier_width", 0.5))
    clear_cw = deck_width - 2.0 * cb_width

    num_lanes = _num_lanes(clear_cw)
    lane_red = _lane_reduction_factor(num_lanes)

    span_lhs = geom.span_lhs
    span_rhs = geom.span_rhs
    avg_span = (span_lhs + span_rhs) / 2.0

    continuity = sup.get("continuity", "Simply Supported")
    is_continuous = continuity == "Deck Continuity"
    n_cont_spans = (
        int(sup.get("num_continuous_spans", 1)) if is_continuous else 1
    )

    n_brg_lhs = len(geom.bearing_coords_lhs)
    n_brg_rhs = len(geom.bearing_coords_rhs)

    # Impact factors for RC bridge
    IF_A = _impact_factor_rc(avg_span, "class_A")
    IF_70R = _impact_factor_rc(avg_span, "70R")

    # Maximum reaction for one Class A train (half span loaded, one side)
    R_A_one = _max_reaction_simply_supported(
        avg_span, _CLASS_A_AXLES, _CLASS_A_SPACINGS,
    )
    R_A_one *= (1.0 + IF_A)

    # Maximum reaction for one 70R Wheeled (half span loaded, one side)
    R_70R_one = _max_reaction_simply_supported(
        avg_span, _70R_AXLES, _70R_SPACINGS,
    )
    R_70R_one *= (1.0 + IF_70R)

    # For continuous spans, the pier reaction is amplified by the continuity
    # support coefficient (approx 1.10 for intermediate pier)
    cont_factor = 1.10 if is_continuous else 1.0

    # Check for user-provided manual LL reactions
    ll_manual = config.get("live_load", None)

    ll_cases: list[LLCase] = []

    # Maximum transverse eccentricity for eccentric loading
    max_y_lhs = geom.btl_lhs.max_y_bearing
    max_y_rhs = geom.btl_rhs.max_y_bearing

    # Determine governing vehicle configuration for max P
    # (a) n lanes of Class A
    R_max_A_total = num_lanes * R_A_one * lane_red * cont_factor
    config_a = f"{num_lanes}xClass A"

    # (b) 1x70R + (n-2) x Class A (70R occupies 2 lanes)
    if num_lanes >= 2:
        n_extra_A = max(0, num_lanes - 2)
        R_max_70R_total = (
            R_70R_one + n_extra_A * R_A_one
        ) * lane_red * cont_factor
        config_b = f"1x70R + {n_extra_A}xClass A"
    else:
        R_max_70R_total = 0.0
        config_b = ""

    if R_max_70R_total > R_max_A_total:
        max_P_total = R_max_70R_total
        max_P_config = config_b
    else:
        max_P_total = R_max_A_total
        max_P_config = config_a

    # --- Case 1: Max P (maximum vertical load) ---
    # Both spans loaded symmetrically -> no ML, no MT
    if ll_manual and "max_P_lhs" in ll_manual and ll_manual["max_P_lhs"] > 0:
        P_lhs_c1 = float(ll_manual["max_P_lhs"])
        P_rhs_c1 = float(ll_manual.get("max_P_rhs", P_lhs_c1))
    else:
        P_lhs_c1 = max_P_total / 2.0
        P_rhs_c1 = max_P_total / 2.0

    lhs_c1 = ForceVector(P=P_lhs_c1)
    rhs_c1 = ForceVector(P=P_rhs_c1)
    ll_cases.append(LLCase(
        name="Max P", lhs=lhs_c1, rhs=rhs_c1,
        force_at_pier_base=lhs_c1 + rhs_c1,
    ))

    # --- Case 2: Max ML (maximum longitudinal moment) ---
    # One side loaded, other unloaded -> differential creates MT about trans axis
    # Longitudinal eccentricity = bearing offset from pier centre
    if ll_manual and "max_ML_lhs" in ll_manual and ll_manual["max_ML_lhs"] > 0:
        P_lhs_c2 = float(ll_manual["max_ML_lhs"])
        P_rhs_c2 = float(ll_manual.get("max_ML_rhs", 0.0))
    else:
        P_lhs_c2 = max_P_total
        P_rhs_c2 = 0.0

    # Longitudinal eccentricity from pier centre to bearing line
    ecc_long = abs(geom.bearing_coords_lhs[0][0]) if n_brg_lhs > 0 else 0.0
    MT_c2 = P_lhs_c2 * ecc_long  # moment about transverse axis

    lhs_c2 = ForceVector(P=P_lhs_c2, MT=MT_c2)
    rhs_c2 = ForceVector(P=P_rhs_c2)
    ll_cases.append(LLCase(
        name="Max ML", lhs=lhs_c2, rhs=rhs_c2,
        force_at_pier_base=lhs_c2 + rhs_c2,
    ))

    # --- Case 3: Max MT (maximum transverse moment) ---
    # Eccentric loading in transverse direction
    if ll_manual and "max_MT_lhs" in ll_manual and ll_manual["max_MT_lhs"] > 0:
        P_lhs_c3 = float(ll_manual["max_MT_lhs"])
        P_rhs_c3 = float(ll_manual.get("max_MT_rhs", P_lhs_c3))
        ML_c3 = float(ll_manual.get("max_MT_moment", 0.0))
    else:
        P_lhs_c3 = max_P_total / 2.0
        P_rhs_c3 = max_P_total / 2.0
        # Approximate eccentricity: CG of LL at ~50% of max bearing spread
        ecc_trans = max(max_y_lhs, max_y_rhs) * 0.5
        ML_c3 = max_P_total * ecc_trans  # moment about longitudinal axis

    lhs_c3 = ForceVector(P=P_lhs_c3, ML=ML_c3 / 2.0)
    rhs_c3 = ForceVector(P=P_rhs_c3, ML=ML_c3 / 2.0)
    ll_cases.append(LLCase(
        name="Max MT", lhs=lhs_c3, rhs=rhs_c3,
        force_at_pier_base=lhs_c3 + rhs_c3,
    ))

    # --- Case 4: Congestion ---
    # IRC 6 Cl.204.4 -- congestion loading using Class A only
    cf = _congestion_factor(avg_span)
    R_cong_per_lane = _max_reaction_simply_supported(
        avg_span, _CLASS_A_AXLES, _CLASS_A_SPACINGS,
    )  # without impact factor
    R_cong_total = num_lanes * R_cong_per_lane * lane_red * cf * cont_factor

    if ll_manual and "cong_lhs" in ll_manual and ll_manual["cong_lhs"] > 0:
        P_lhs_c4 = float(ll_manual["cong_lhs"])
        P_rhs_c4 = float(ll_manual.get("cong_rhs", P_lhs_c4))
    else:
        P_lhs_c4 = R_cong_total / 2.0
        P_rhs_c4 = R_cong_total / 2.0

    # Congestion: symmetric, no eccentricity, no braking
    lhs_c4 = ForceVector(P=P_lhs_c4)
    rhs_c4 = ForceVector(P=P_rhs_c4)
    ll_cases.append(LLCase(
        name="Congestion", lhs=lhs_c4, rhs=rhs_c4,
        force_at_pier_base=lhs_c4 + rhs_c4,
    ))

    return ll_cases


# ---------------------------------------------------------------------------
# 4. Braking Force
# ---------------------------------------------------------------------------

def _compute_braking(
    config: dict[str, Any],
    geom: GeometryResults,
    ll_cases: list[LLCase],
    irc: dict[str, Any],
) -> list[dict[str, ForceVector]]:
    """Compute braking force at pier base for each LL case.

    IRC 6 Cl.206:
    - Class A: 20% of first vehicle of each train + 5% of all vehicles
    - 70R: 20% of first 70R axle loads
    - Congestion: no braking (5% nominal if required)

    Braking lever arm = 1.2m above deck + depth of superstructure
    including slab + WC thickness = application height above bearing level.

    The braking force is shared among piers for continuous spans.
    Longitudinal moment at pier base = braking_force * lever_arm_to_base.
    Vertical reaction couple = +/- ML / span_cc_bearing / n_bearings_per_side.
    """
    sup = config["superstructure"]
    braking_tbl = irc.get("braking", {})

    depth_super = float(sup["depth_incl_slab"])
    wc_m = float(sup["wearing_coat_thickness"]) / 1000.0
    height_above_road = float(braking_tbl.get("height_above_road", 1.2))

    pier_base_level = geom.pilecap_top_level
    frl = config["levels"]["frl"]

    # Lever arm from road surface through superstructure to pier base
    # Application point is 1.2m above road surface.
    # Height of application above bearing = 1.2 + depth_super + wc_thickness
    lever_above_bearing = height_above_road + depth_super + wc_m

    # Bearing level to pier base
    btl_avg = (
        geom.btl_lhs.bearing_top_level + geom.btl_rhs.bearing_top_level
    ) / 2.0
    lever_bearing_to_base = btl_avg - pier_base_level

    # Total lever arm from braking application to pier base
    lever_arm = lever_above_bearing + lever_bearing_to_base

    continuity = sup.get("continuity", "Simply Supported")
    is_continuous = continuity == "Deck Continuity"
    n_cont_spans = (
        int(sup.get("num_continuous_spans", 1)) if is_continuous else 1
    )

    deck_width = float(sup["deck_width"])
    cb_width = float(sup.get("crash_barrier_width", 0.5))
    clear_cw = deck_width - 2.0 * cb_width
    num_lanes = _num_lanes(clear_cw)

    n_brg_lhs = len(geom.bearing_coords_lhs)
    n_brg_rhs = len(geom.bearing_coords_rhs)
    n_brg_total = n_brg_lhs + n_brg_rhs

    # Span centre-to-centre of bearings (longitudinal distance between
    # LHS and RHS bearing lines)
    if n_brg_lhs > 0 and n_brg_rhs > 0:
        span_cc_brg = abs(
            geom.bearing_coords_rhs[0][0] - geom.bearing_coords_lhs[0][0]
        )
    else:
        span_cc_brg = 1.65  # default from sample: 2 * 0.825

    braking_results: list[dict[str, ForceVector]] = []

    for idx, ll_case in enumerate(ll_cases):
        case_name = ll_case.name

        # Determine total braking force based on vehicle configuration
        if case_name == "Congestion":
            # IRC 6 Cl.206.4: 5% of all vehicles on congestion loading
            bf = 0.05 * num_lanes * _CLASS_A_TOTAL
        elif case_name == "Max P":
            # Governing config: whichever gave max P
            # For 3 lanes Class A: 20% of 1st vehicle per train + 5% rest
            # = 3 * (0.20*554) + (remaining not double-counted)
            # Simplified: 20% of first train + 5% of all trains
            # For nxClassA: 0.20 * CLASS_A_TOTAL + 0.05 * (n-1) * CLASS_A_TOTAL
            bf = (
                0.20 * _CLASS_A_TOTAL
                + 0.05 * max(0, num_lanes - 1) * _CLASS_A_TOTAL
            )
        elif case_name in ("Max ML", "Max MT"):
            # 1x70R + remaining Class A lanes (if 70R governs)
            # or nxClassA (if ClassA governs)
            if num_lanes >= 2:
                # Try 1x70R + ClassA combination
                bf_70r = 0.20 * _70R_TOTAL
                n_extra_A = max(0, num_lanes - 2)
                bf_classA = 0.05 * n_extra_A * _CLASS_A_TOTAL
                bf_combo = bf_70r + bf_classA

                # Pure Class A
                bf_pure_A = (
                    0.20 * _CLASS_A_TOTAL
                    + 0.05 * max(0, num_lanes - 1) * _CLASS_A_TOTAL
                )

                bf = max(bf_combo, bf_pure_A)
            else:
                bf = 0.20 * _CLASS_A_TOTAL
        else:
            bf = 0.20 * _CLASS_A_TOTAL

        # Share braking among piers for continuous spans
        bf_at_pier = bf / n_cont_spans

        # Longitudinal force per bearing
        fl_per_brg = bf_at_pier / n_brg_total if n_brg_total > 0 else 0.0

        # Longitudinal moment at pier base
        MT_braking = bf_at_pier * lever_arm

        # Vertical reaction couple from braking moment
        # P_couple = +/- MT / span_cc_bearing for each side
        if span_cc_brg > 1e-6 and n_brg_lhs > 0:
            P_couple_per_brg = MT_braking / span_cc_brg / n_brg_lhs
        else:
            P_couple_per_brg = 0.0

        # LHS bearings get +P (braking pushes down on leading side)
        # RHS bearings get -P (uplift on trailing side)
        lhs_fv = ForceVector(
            P=P_couple_per_brg * n_brg_lhs,
            FL=fl_per_brg * n_brg_lhs,
            MT=MT_braking / 2.0,
        )
        rhs_fv = ForceVector(
            P=-P_couple_per_brg * n_brg_rhs,
            FL=fl_per_brg * n_brg_rhs,
            MT=MT_braking / 2.0,
        )

        braking_results.append({
            "lhs": lhs_fv,
            "rhs": rhs_fv,
            "total": ForceVector(FL=bf_at_pier, MT=MT_braking),
        })

    return braking_results


# ---------------------------------------------------------------------------
# 5. Centrifugal Force
# ---------------------------------------------------------------------------

def _compute_centrifugal(
    config: dict[str, Any],
    geom: GeometryResults,
    ll_cases: list[LLCase],
) -> list[dict[str, ForceVector]]:
    """Centrifugal force per IRC 6 Cl.207.

    C = W * V^2 / (127 * R), applied at 1.2m above road surface.
    Not applicable for straight bridges or congestion case.

    The vertical reaction couple from the centrifugal moment:
        P_couple = C * lever_arm / (2 * max_bearing_spacing)
    """
    sup = config["superstructure"]
    R = float(sup.get("radius_of_curvature", 0.0))

    _zero_case: dict[str, ForceVector] = {
        "lhs": ForceVector(),
        "rhs": ForceVector(),
        "total": ForceVector(),
    }

    centrifugal_results: list[dict[str, ForceVector]] = []

    if R <= 0.0 or R > 50000.0:
        # Straight bridge or negligible curvature -- no centrifugal force
        for _ in ll_cases:
            centrifugal_results.append({
                "lhs": ForceVector(),
                "rhs": ForceVector(),
                "total": ForceVector(),
            })
        return centrifugal_results

    V = _DESIGN_SPEED_KMPH  # km/h
    pier_base_level = geom.pilecap_top_level
    frl = config["levels"]["frl"]
    depth_super = float(sup["depth_incl_slab"])
    wc_m = float(sup["wearing_coat_thickness"]) / 1000.0

    # Lever arm from centrifugal application (1.2m above road) to pier base
    # Application at 1.2m above road surface
    height_above_road = 1.2  # m
    btl_avg = (
        geom.btl_lhs.bearing_top_level + geom.btl_rhs.bearing_top_level
    ) / 2.0
    lever_above_bearing = height_above_road + depth_super + wc_m
    lever_bearing_to_base = btl_avg - pier_base_level
    lever_arm = lever_above_bearing + lever_bearing_to_base

    continuity = sup.get("continuity", "Simply Supported")
    is_continuous = continuity == "Deck Continuity"
    n_cont_spans = (
        int(sup.get("num_continuous_spans", 1)) if is_continuous else 1
    )

    # Maximum transverse bearing spacing for vertical couple
    max_y_lhs = geom.btl_lhs.max_y_bearing
    max_y_rhs = geom.btl_rhs.max_y_bearing
    max_spacing = max(max_y_lhs, max_y_rhs)

    n_brg_lhs = len(geom.bearing_coords_lhs)
    n_brg_rhs = len(geom.bearing_coords_rhs)
    n_brg_total = n_brg_lhs + n_brg_rhs

    for idx, ll_case in enumerate(ll_cases):
        if "Congestion" in ll_case.name:
            # No centrifugal for congestion case
            centrifugal_results.append({
                "lhs": ForceVector(),
                "rhs": ForceVector(),
                "total": ForceVector(),
            })
            continue

        W_ll = ll_case.lhs.P + ll_case.rhs.P  # total LL on the pier

        # Centrifugal force: C = W * V^2 / (127 * R)
        C = W_ll * V ** 2 / (127.0 * R)
        C_at_pier = C / n_cont_spans

        # Transverse force and moment at pier base
        ML_centrifugal = C_at_pier * lever_arm

        # Vertical reaction couple from centrifugal moment
        if max_spacing > 1e-6:
            P_couple = C_at_pier * lever_arm / (2.0 * max_spacing)
        else:
            P_couple = 0.0

        total_fv = ForceVector(
            FT=C_at_pier,
            ML=ML_centrifugal,
        )

        # Split equally between LHS and RHS
        lhs_fv = ForceVector(
            P=P_couple,
            FT=C_at_pier / 2.0,
            ML=ML_centrifugal / 2.0,
        )
        rhs_fv = ForceVector(
            P=-P_couple,
            FT=C_at_pier / 2.0,
            ML=ML_centrifugal / 2.0,
        )

        centrifugal_results.append({
            "lhs": lhs_fv,
            "rhs": rhs_fv,
            "total": total_fv,
        })

    return centrifugal_results


# ---------------------------------------------------------------------------
# 6. Wind Load
# ---------------------------------------------------------------------------

def _compute_wind(
    config: dict[str, Any],
    geom: GeometryResults,
    irc: dict[str, Any],
) -> tuple[
    ForceVector, ForceVector,   # wind_on_super_lhs, _rhs
    ForceVector, ForceVector,   # wind_on_ll_lhs, _rhs
    ForceVector, ForceVector,   # wind_on_pier, wind_on_piercap
    ForceVector, ForceVector,   # wind_vertical_lhs, _rhs
]:
    """Wind loads per IRC 6 Cl.209.

    Returns separate force vectors for each wind component:
      (wind_super_lhs, wind_super_rhs, wind_ll_lhs, wind_ll_rhs,
       wind_pier, wind_piercap, wind_vert_lhs, wind_vert_rhs)

    Wind pressure: Pz = 0.6 * Vz^2 (N/m^2), where Vz = k2 * Vb
    Transverse: Pz * area_elevation * CD * G / 1000 (kN)
    Longitudinal: 25% of transverse
    Vertical uplift: Pz * area_plan * CL * G / 1000 (kN)

    CD from b/d ratio interpolated from IRC 6 Table for superstructure.
    CD for pier: 0.5 smooth circular (from IRC table), or from b/d for rect.
    """
    sup = config["superstructure"]
    pcap_cfg = config["pier_cap"]
    wind_cfg = config["wind"]

    Vb = float(wind_cfg["basic_speed"])  # m/s
    pier_base_level = geom.pilecap_top_level
    frl = config["levels"]["frl"]
    gl = config["levels"]["gl"]

    # Average deck height above ground for wind height factor
    deck_level = frl  # top of deck approximately at FRL
    deck_height = deck_level - gl

    # Height factor at deck level
    k2_deck = _get_height_factor(max(deck_height, 10.0), irc)

    # Hourly mean wind pressure at deck level
    Vz_deck = k2_deck * Vb
    Pz_deck = 0.6 * Vz_deck ** 2  # N/m^2 = Pa

    deck_type = sup["type"]
    depth_super = float(sup["depth_incl_slab"])
    cb_height = float(sup.get("crash_barrier_height", 1.1))
    nb_height = float(sup.get("noise_barrier_height", 0.0))
    deck_width = float(sup["deck_width"])

    span_lhs = geom.span_lhs
    span_rhs = geom.span_rhs

    G = _GUST_FACTOR

    # --- Drag coefficient for superstructure from IRC 6 tables ---
    wind_tbl = irc["wind"]["drag_coefficients"]
    # Map deck type to table key
    deck_key_map = {
        "U Girder": "u_girder",
        "Box Girder": "box_girder",
        "I Girder": "i_girder",
        "Slab": "deck_slab",
    }
    deck_key = deck_key_map.get(deck_type, "u_girder")
    CD_super = float(wind_tbl.get(deck_key, 2.0))

    # --- Wind on superstructure (per span) ---
    # Elevation area (side view) per half-span at pier
    elevation_height = depth_super + cb_height + nb_height

    def _wind_super_for_span(span: float) -> tuple[ForceVector, ForceVector]:
        """Wind on superstructure for half-span at pier.

        Returns (transverse+longitudinal ForceVector, vertical ForceVector).
        """
        # Half-span contribution to this pier
        half_span = span / 2.0

        # Elevation area
        A_elev = half_span * elevation_height  # m^2

        # Transverse force on superstructure
        FT = Pz_deck * A_elev * CD_super * G / 1000.0  # kN

        # Longitudinal = 25% of transverse (IRC 6 Cl.209.3.5)
        FL = 0.25 * FT

        # CG of wind on superstructure (at bearing level + depth/2)
        btl_avg = (
            geom.btl_lhs.bearing_top_level + geom.btl_rhs.bearing_top_level
        ) / 2.0
        wind_super_cg = btl_avg + depth_super / 2.0
        lever = wind_super_cg - pier_base_level

        fv_ht = ForceVector(
            FL=FL,
            FT=FT,
            ML=FT * lever,
            MT=FL * lever,
        )

        # Plan area for vertical uplift
        A_plan = deck_width * half_span
        FV = Pz_deck * A_plan * _CL_DECK * G / 1000.0  # kN (uplift)

        fv_vert = ForceVector(P=-FV)  # negative P = uplift

        return fv_ht, fv_vert

    wind_super_lhs, wind_vert_lhs = _wind_super_for_span(span_lhs)
    wind_super_rhs, wind_vert_rhs = _wind_super_for_span(span_rhs)

    # --- Wind on live load (per span) ---
    # LL height above deck = 3m for wind calculation (minus barriers above deck)
    barrier_above_deck = cb_height  # crash barrier height above deck
    ll_exposed_height = max(0.0, 3.0 - barrier_above_deck)
    CD_ll = 1.2  # drag coefficient for vehicles

    # CG of wind on LL: at deck top + barrier + half exposed height
    ll_cg = frl + barrier_above_deck + ll_exposed_height / 2.0
    lever_ll = ll_cg - pier_base_level

    def _wind_ll_for_span(span: float) -> ForceVector:
        """Wind on live load for half-span at pier."""
        half_span = span / 2.0
        A_ll = half_span * ll_exposed_height
        FT = Pz_deck * A_ll * CD_ll * G / 1000.0
        FL = 0.25 * FT
        return ForceVector(
            FL=FL,
            FT=FT,
            ML=FT * lever_ll,
            MT=FL * lever_ll,
        )

    wind_ll_lhs = _wind_ll_for_span(span_lhs)
    wind_ll_rhs = _wind_ll_for_span(span_rhs)

    # --- Wind on pier ---
    h_pier = geom.pier_height
    pier_type = geom.pier_section.pier_type
    D_pier_trans = geom.pier_section.width_trans
    D_pier_long = geom.pier_section.width_long

    # Height factor at pier mid-height
    pier_mid_level = pier_base_level + h_pier / 2.0
    pier_mid_height_above_gl = pier_mid_level - gl
    k2_pier = _get_height_factor(max(pier_mid_height_above_gl, 10.0), irc)
    Vz_pier = k2_pier * Vb
    Pz_pier = 0.6 * Vz_pier ** 2

    # Pier drag coefficient
    if pier_type == "Circular":
        CD_pier = float(
            irc["wind"]["drag_coefficients"]["circular_pier"]["smooth"]
        )
    else:
        # Rectangular: use b/d ratio to pick coefficient
        b_d = D_pier_trans / D_pier_long if D_pier_long > 0 else 2.0
        if b_d >= 6:
            CD_pier = float(
                irc["wind"]["drag_coefficients"]["rectangular_pier"]["b_d_6"]
            )
        else:
            CD_pier = float(
                irc["wind"]["drag_coefficients"]["rectangular_pier"]["b_d_2"]
            )

    # Exposed pier height: from pilecap top to piercap bottom
    # If GL is above pilecap top, the pier below GL is underground
    pier_exposed_bottom = max(pier_base_level, gl)
    pier_exposed_top = geom.piercap_bottom_level
    pier_exposed_ht = max(0.0, pier_exposed_top - pier_exposed_bottom)

    # Transverse wind on pier (using trapezoidal pressure approximation)
    A_pier_trans = D_pier_long * pier_exposed_ht
    FT_pier = Pz_pier * A_pier_trans * CD_pier * G / 1000.0

    # Longitudinal wind on pier
    A_pier_long = D_pier_trans * pier_exposed_ht
    FL_pier = Pz_pier * A_pier_long * CD_pier * G / 1000.0

    # Lever arm: CG of exposed pier above pier base
    # For uniform pressure, CG at mid-height of exposed portion
    pier_cg_above_base = (
        (pier_exposed_bottom - pier_base_level)
        + pier_exposed_ht / 2.0
    )

    wind_pier = ForceVector(
        FL=FL_pier,
        FT=FT_pier,
        ML=FT_pier * pier_cg_above_base,
        MT=FL_pier * pier_cg_above_base,
    )

    # --- Wind on pier cap ---
    pcap_depth = float(pcap_cfg["depth_max"])
    pcap_length = float(pcap_cfg["length_trans"])
    pcap_width = float(pcap_cfg["width_long"])

    # Pier cap visible height: above GL
    pcap_visible_top = geom.piercap_top_level
    pcap_visible_bottom = max(geom.piercap_bottom_level, gl)
    pcap_visible_ht = max(0.0, pcap_visible_top - pcap_visible_bottom)

    # CD for pier cap: approximate as rectangular with t/b ratio
    t_b_pcap = pcap_length / pcap_width if pcap_width > 0 else 1.0
    if t_b_pcap >= 6:
        CD_pcap = 1.3
    elif t_b_pcap >= 2:
        # Linear interpolation between 1.5 at t/b=2 and 1.3 at t/b=6
        CD_pcap = 1.5 - (t_b_pcap - 2.0) * (1.5 - 1.3) / (6.0 - 2.0)
    else:
        CD_pcap = 1.5

    # Transverse wind on pier cap
    A_pcap_trans = pcap_width * pcap_visible_ht
    FT_pcap = Pz_deck * A_pcap_trans * CD_pcap * G / 1000.0

    # Longitudinal wind on pier cap
    A_pcap_long = pcap_length * pcap_visible_ht
    FL_pcap = Pz_deck * A_pcap_long * CD_pcap * G / 1000.0

    # CG of pier cap above pier base
    pcap_cg = geom.piercap_bottom_level + pcap_depth / 2.0
    lever_pcap = pcap_cg - pier_base_level

    wind_piercap = ForceVector(
        FL=FL_pcap,
        FT=FT_pcap,
        ML=FT_pcap * lever_pcap,
        MT=FL_pcap * lever_pcap,
    )

    return (
        wind_super_lhs, wind_super_rhs,
        wind_ll_lhs, wind_ll_rhs,
        wind_pier, wind_piercap,
        wind_vert_lhs, wind_vert_rhs,
    )


# ---------------------------------------------------------------------------
# 7. Seismic (Elastic Seismic Coefficient Method)
# ---------------------------------------------------------------------------

def _compute_seismic(
    config: dict[str, Any],
    geom: GeometryResults,
    irc: dict[str, Any],
    total_dl_super: float,
    total_sidl: float,
    total_wc: float,
    wt_piercap: float,
    wt_pier: float,
) -> tuple[
    ForceVector, ForceVector,   # super_long, super_trans
    ForceVector, ForceVector,   # pier_long, pier_trans
    ForceVector, ForceVector,   # piercap_long, piercap_trans
    ForceVector,                # hydrodynamic
    float, float,               # Ah_long, Ah_trans
]:
    """Seismic forces per IRC 6 Cl.211 / IRC SP-114.

    Returns individual seismic force vectors for each component and the
    horizontal seismic coefficients Ah_long, Ah_trans.

    Steps:
    1. Compute pier lateral stiffness: K = 3*E*I / H^3 (cantilever)
    2. Compute natural period: T = 2*pi * sqrt(W / (g*K))
    3. Get Sa/g from response spectrum per soil type
    4. Ah = Z/2 * I/R * Sa/g
    5. Compute forces on each component and moments at pier base
    """
    eq_cfg = config["seismic"]
    mat = config["materials"]
    fnd = config["foundation"]

    zone_str = eq_cfg["zone"]
    zone_factors = irc["seismic"]["zone_factors"]
    Z = float(zone_factors[zone_str])
    I_factor = float(eq_cfg["importance_factor"])
    R_factor = float(eq_cfg["response_reduction"])
    soil_type = eq_cfg.get("soil_type", "medium")

    density = float(mat.get("concrete_density", _RHO_CONCRETE))
    pier_base_level = geom.pilecap_top_level
    pcap_cfg = config["pier_cap"]

    # -- Pier stiffness --
    fck_pier = float(mat["pier"]["fck"])
    Ecm = _ecm_mpa(fck_pier)     # MPa
    Ecm_kpa = Ecm * 1000.0       # kPa = kN/m^2

    h_pier = geom.pier_height
    I_pier_xx = geom.pier_section.inertia_xx  # m^4 about longitudinal axis
    I_pier_yy = geom.pier_section.inertia_yy  # m^4 about transverse axis

    # Cracked stiffness: use 0.7 * I_gross (IRC SP-114 recommendation)
    cracked_factor = 0.70

    # Fixity depth contribution to effective height
    Lf = geom.pile_props.fixity_depth_Lf
    H_eff = h_pier + Lf  # effective cantilever height

    # Pier lateral stiffness (cantilever model): K = 3*E*I_cr / H_eff^3
    if H_eff > 0:
        K_trans = (
            3.0 * Ecm_kpa * cracked_factor * I_pier_xx / H_eff ** 3
        )  # kN/m
        K_long = (
            3.0 * Ecm_kpa * cracked_factor * I_pier_yy / H_eff ** 3
        )  # kN/m
    else:
        K_trans = 1e9
        K_long = 1e9

    # -- Seismic weight --
    # Weight above pier base: DL super + SIDL + WC + piercap + half pier
    W_total = (
        total_dl_super + total_sidl + total_wc
        + wt_piercap + 0.5 * wt_pier
    )

    # -- Natural period --
    # T = 2*pi * sqrt(W / (g * K))  where W is in kN, K in kN/m
    if K_trans > 0 and W_total > 0:
        T_trans = 2.0 * math.pi * math.sqrt(W_total / (_G * K_trans))
    else:
        T_trans = 0.5

    if K_long > 0 and W_total > 0:
        T_long = 2.0 * math.pi * math.sqrt(W_total / (_G * K_long))
    else:
        T_long = 0.5

    # -- Spectral acceleration --
    Sa_g_trans = _sa_over_g(T_trans, soil_type)
    Sa_g_long = _sa_over_g(T_long, soil_type)

    # -- Horizontal seismic coefficient --
    Ah_trans = Z * I_factor * Sa_g_trans / (2.0 * R_factor)
    Ah_long = Z * I_factor * Sa_g_long / (2.0 * R_factor)

    # -- Vertical seismic coefficient --
    Av = (2.0 / 3.0) * max(Ah_trans, Ah_long)

    # -- Lever arms --
    btl_avg = (
        geom.btl_lhs.bearing_top_level + geom.btl_rhs.bearing_top_level
    ) / 2.0
    lever_super = btl_avg - pier_base_level       # bearing to pier base
    lever_pcap = (
        geom.piercap_bottom_level + pcap_cfg["depth_max"] / 2.0
        - pier_base_level
    )  # piercap CG to pier base
    lever_pier_cg = h_pier / 2.0                  # pier CG above pier base

    # -- Weight of superstructure (DL + SIDL + WC) for seismic --
    W_super = total_dl_super + total_sidl + total_wc

    # -- Seismic on superstructure (longitudinal direction) --
    F_super_long = Ah_long * W_super
    seismic_super_long = ForceVector(
        P=Av * W_super,
        FL=F_super_long,
        MT=F_super_long * lever_super,
    )

    # -- Seismic on superstructure (transverse direction) --
    F_super_trans = Ah_trans * W_super
    seismic_super_trans = ForceVector(
        P=Av * W_super,
        FT=F_super_trans,
        ML=F_super_trans * lever_super,
    )

    # -- Seismic on pier (longitudinal direction) --
    F_pier_long = Ah_long * wt_pier
    seismic_pier_long = ForceVector(
        P=Av * wt_pier,
        FL=F_pier_long,
        MT=F_pier_long * lever_pier_cg,
    )

    # -- Seismic on pier (transverse direction) --
    F_pier_trans = Ah_trans * wt_pier
    seismic_pier_trans = ForceVector(
        P=Av * wt_pier,
        FT=F_pier_trans,
        ML=F_pier_trans * lever_pier_cg,
    )

    # -- Seismic on pier cap (longitudinal direction) --
    F_pcap_long = Ah_long * wt_piercap
    seismic_piercap_long = ForceVector(
        P=Av * wt_piercap,
        FL=F_pcap_long,
        MT=F_pcap_long * lever_pcap,
    )

    # -- Seismic on pier cap (transverse direction) --
    F_pcap_trans = Ah_trans * wt_piercap
    seismic_piercap_trans = ForceVector(
        P=Av * wt_piercap,
        FT=F_pcap_trans,
        ML=F_pcap_trans * lever_pcap,
    )

    # -- Hydrodynamic force (if applicable) --
    # Per IRC 6 Cl.219.5: applicable when pier is submerged
    hfl = config["levels"]["hfl"]
    seismic_hydro = ForceVector()

    if hfl > pier_base_level:
        # Submerged height of pier
        water_level = min(hfl, geom.piercap_bottom_level)
        submerged_ht = max(0.0, water_level - pier_base_level)

        if submerged_ht > 0.0:
            water_density = float(mat.get("water_density", 10.0))  # kN/m^3
            # Hydrodynamic pressure coefficient Ce depends on shape:
            # For circular pier: Ce envelope from Zangar's curves
            # Simplified: Ce_avg ~ 0.73 for circular, ~0.67 for rectangular
            if geom.pier_section.pier_type == "Circular":
                Ce_avg = 0.73
                D_hydro = geom.pier_section.width_trans
            else:
                Ce_avg = 0.67
                D_hydro = geom.pier_section.width_trans

            # Hydrodynamic force per IRC 6 Cl.219.5:
            # F_hydro = Ce * Ah * water_density * D * H^2 / 2
            # where H is submerged height, D is pier width perpendicular
            # to flow direction
            # The factor accounts for the pressure distribution
            Ah_max = max(Ah_long, Ah_trans)
            F_hydro = (
                Ce_avg * Ah_max * water_density * D_hydro
                * submerged_ht ** 2 / 2.0
            )

            # CG of hydrodynamic force: approximately at 0.4H from base
            lever_hydro = 0.4 * submerged_ht

            # Apply in both directions
            seismic_hydro = ForceVector(
                FL=F_hydro * Ah_long / Ah_max if Ah_max > 0 else 0.0,
                FT=F_hydro * Ah_trans / Ah_max if Ah_max > 0 else 0.0,
                ML=F_hydro * Ah_trans / Ah_max * lever_hydro
                if Ah_max > 0 else 0.0,
                MT=F_hydro * Ah_long / Ah_max * lever_hydro
                if Ah_max > 0 else 0.0,
            )

    return (
        seismic_super_long, seismic_super_trans,
        seismic_pier_long, seismic_pier_trans,
        seismic_piercap_long, seismic_piercap_trans,
        seismic_hydro,
        Ah_long, Ah_trans,
    )


# ---------------------------------------------------------------------------
# 8. Temperature
# ---------------------------------------------------------------------------

def _compute_temperature(
    config: dict[str, Any],
    geom: GeometryResults,
    irc: dict[str, Any],
) -> ForceVector:
    """Temperature-induced longitudinal force at pier base.

    For elastomeric bearings: force = bearing shear stiffness * thermal movement.
    Movement = alpha * delta_T * span/2 (expansion from pier location).

    Returns a single ForceVector (longitudinal).  For rise/fall cases the
    caller can negate.

    IRC 6 Cl.215: total temperature range = rise + fall.
    """
    sup = config["superstructure"]
    temp_tbl = irc.get("temperature", {})

    temp_rise = float(temp_tbl.get("rise", 30.0))
    temp_fall = float(temp_tbl.get("fall", 20.0))
    alpha = float(temp_tbl.get("coefficient", 0.000012))

    avg_span = (geom.span_lhs + geom.span_rhs) / 2.0

    # Maximum thermal movement at pier (governs for rise case)
    # Free expansion = alpha * dT * span/2
    max_dT = max(temp_rise, temp_fall)
    expansion = alpha * max_dT * avg_span / 2.0  # m

    # Bearing shear stiffness
    # Typical shear stiffness ~1500 kN/m per bearing for elastomeric
    n_brg_total = len(geom.bearing_coords_lhs) + len(geom.bearing_coords_rhs)
    bearing_type = config.get("bearings", {}).get("type", "Elastomeric")

    if bearing_type.lower() == "pot":
        K_per = 500.0  # kN/m equivalent for pot-PTFE
    else:
        K_per = 1500.0  # kN/m for elastomeric

    K_total = K_per * n_brg_total
    F_temp = K_total * expansion  # kN

    # Lever arm from bearing to pier base
    pier_base_level = geom.pilecap_top_level
    btl_avg = (
        geom.btl_lhs.bearing_top_level + geom.btl_rhs.bearing_top_level
    ) / 2.0
    lever = btl_avg - pier_base_level

    return ForceVector(
        FL=F_temp,
        MT=F_temp * lever,
    )


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def calculate_loads(
    config: dict[str, Any],
    geometry: GeometryResults,
) -> LoadResults:
    """Compute all design loads for the bridge substructure.

    This is the primary public API of the loads module.  It reads the
    validated project configuration and the pre-computed geometry, then
    calculates every load type specified by IRC 6:2017, resolving all
    forces to the pier base level.

    Parameters
    ----------
    config : dict
        Validated configuration dictionary from :func:`input_parser.parse_input`.
    geometry : GeometryResults
        Pre-computed geometry from :func:`geometry.calculate_geometry`.

    Returns
    -------
    LoadResults
        Dataclass containing all computed force vectors and per-bearing
        reactions, ready for consumption by :mod:`load_combinations`.
    """
    irc = load_irc_tables()

    # ------------------------------------------------------------------
    # 1. Dead Load
    # ------------------------------------------------------------------
    dl_piercap, dl_pier, dl_super_lhs, dl_super_rhs, brg_dl = (
        _compute_dead_loads(config, geometry)
    )

    # ------------------------------------------------------------------
    # 2. SIDL
    # ------------------------------------------------------------------
    sidl_lhs, sidl_rhs, wc_lhs, wc_rhs, brg_sidl = _compute_sidl(
        config, geometry,
    )

    # ------------------------------------------------------------------
    # 3. Live Load
    # ------------------------------------------------------------------
    ll_cases = _compute_live_load(config, geometry)

    # ------------------------------------------------------------------
    # 4. Braking
    # ------------------------------------------------------------------
    braking_cases = _compute_braking(
        config, geometry, ll_cases, irc,
    )

    # ------------------------------------------------------------------
    # 5. Centrifugal
    # ------------------------------------------------------------------
    centrifugal_cases = _compute_centrifugal(
        config, geometry, ll_cases,
    )

    # ------------------------------------------------------------------
    # 6. Wind
    # ------------------------------------------------------------------
    (
        wind_super_lhs, wind_super_rhs,
        wind_ll_lhs, wind_ll_rhs,
        wind_pier, wind_piercap,
        wind_vert_lhs, wind_vert_rhs,
    ) = _compute_wind(config, geometry, irc)

    # ------------------------------------------------------------------
    # 7. Seismic
    # ------------------------------------------------------------------
    total_dl_super = dl_super_lhs.P + dl_super_rhs.P
    total_sidl = sidl_lhs.P + sidl_rhs.P
    total_wc = wc_lhs.P + wc_rhs.P

    (
        seismic_super_long, seismic_super_trans,
        seismic_pier_long, seismic_pier_trans,
        seismic_piercap_long, seismic_piercap_trans,
        seismic_hydro,
        Ah_long, Ah_trans,
    ) = _compute_seismic(
        config, geometry, irc,
        total_dl_super=total_dl_super,
        total_sidl=total_sidl,
        total_wc=total_wc,
        wt_piercap=dl_piercap.P,
        wt_pier=dl_pier.P,
    )

    # ------------------------------------------------------------------
    # 8. Temperature
    # ------------------------------------------------------------------
    temp_long = _compute_temperature(config, geometry, irc)

    # ------------------------------------------------------------------
    # Per-bearing reactions (collated)
    # ------------------------------------------------------------------
    bearing_reactions: dict[str, dict[str, list[float]]] = {
        "DL": brg_dl,
        **brg_sidl,
    }
    # Add LL per-bearing for each case
    n_brg_lhs = len(geometry.bearing_coords_lhs)
    n_brg_rhs = len(geometry.bearing_coords_rhs)
    for ll_case in ll_cases:
        P_lhs = ll_case.lhs.P
        P_rhs = ll_case.rhs.P
        bearing_reactions[f"LL_{ll_case.name}"] = {
            "lhs": [P_lhs / n_brg_lhs] * n_brg_lhs if n_brg_lhs > 0 else [],
            "rhs": [P_rhs / n_brg_rhs] * n_brg_rhs if n_brg_rhs > 0 else [],
        }

    # ------------------------------------------------------------------
    # Summary totals
    # ------------------------------------------------------------------
    total_dl = dl_piercap + dl_pier + dl_super_lhs + dl_super_rhs
    total_sidl_fv = sidl_lhs + sidl_rhs
    total_wc_fv = wc_lhs + wc_rhs

    # ------------------------------------------------------------------
    # Assemble
    # ------------------------------------------------------------------
    return LoadResults(
        # Dead loads
        dl_piercap=dl_piercap,
        dl_pier=dl_pier,
        dl_super_lhs=dl_super_lhs,
        dl_super_rhs=dl_super_rhs,
        # SIDL
        sidl_lhs=sidl_lhs,
        sidl_rhs=sidl_rhs,
        wc_lhs=wc_lhs,
        wc_rhs=wc_rhs,
        # Live load
        ll_cases=ll_cases,
        # Braking
        braking_cases=braking_cases,
        # Centrifugal
        centrifugal_cases=centrifugal_cases,
        # Wind
        wind_on_super_lhs=wind_super_lhs,
        wind_on_super_rhs=wind_super_rhs,
        wind_on_ll_lhs=wind_ll_lhs,
        wind_on_ll_rhs=wind_ll_rhs,
        wind_on_pier=wind_pier,
        wind_on_piercap=wind_piercap,
        wind_vertical_lhs=wind_vert_lhs,
        wind_vertical_rhs=wind_vert_rhs,
        # Seismic
        seismic_super_long=seismic_super_long,
        seismic_super_trans=seismic_super_trans,
        seismic_pier_long=seismic_pier_long,
        seismic_pier_trans=seismic_pier_trans,
        seismic_piercap_long=seismic_piercap_long,
        seismic_piercap_trans=seismic_piercap_trans,
        seismic_hydrodynamic=seismic_hydro,
        Ah_long=Ah_long,
        Ah_trans=Ah_trans,
        # Temperature
        temp_long=temp_long,
        # Summaries
        total_dl_at_pier_base=total_dl,
        total_sidl_at_pier_base=total_sidl_fv,
        total_wc_at_pier_base=total_wc_fv,
        # Per-bearing
        bearing_reactions=bearing_reactions,
    )

"""Pier cap design (flexure, shear, crack width) per IRC 112:2020.

Designs the pier cap as a transverse cantilever beam projecting from
both faces of the pier.  Checks are performed at two critical sections:

* **Section A** -- at the face of the pier (maximum hogging moment).
* **Section B** -- at a curtailment location midway along the cantilever
  (or at the outermost bearing, whichever is more critical).

The tapered geometry is fully accounted for: the depth varies linearly
from ``depth_max`` (at the pier face) to ``depth_min`` (at the
cantilever tip).

Design checks
-------------
1. **ULS Flexure** -- IRC 112 Cl. 8.2 (rectangular stress-block).
2. **ULS Shear**   -- IRC 112 Cl. 10.3 (concrete capacity + stirrup
   design if required).
3. **SLS Crack Width** -- IRC 112 Cl. 12.3.4 (direct calculation).
4. **Minimum Reinforcement** -- IRC 112 Cl. 16.5.1.1.

Units convention
----------------
Forces in **kN**, moments in **kN.m**, dimensions in **mm** for
section-level calculations (converted from metres as needed).

Key references
--------------
* IRC 112:2020  -- Code of Practice for Concrete Road Bridges
* IRC 6:2017    -- Standard Specifications and Code of Practice for
  Road Bridges, Section II: Loads and Load Combinations
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from .loads import ForceVector
from .geometry import GeometryResults, PierCapSection
from .load_combinations import CombinationResult, CombinationResults
from .materials import (
    get_concrete_properties,
    get_steel_properties,
    ConcreteProperties,
    SteelProperties,
)
from .utils import load_irc_tables


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_XU_D_LIMIT_FE500: float = 0.46
"""Maximum permitted xu/d ratio for Fe 500 / Fe 550 steel (IRC 112)."""

_DEFAULT_BAR_DIA: float = 32.0
"""Default main bar diameter (mm) for reinforcement estimation."""

_DEFAULT_STIRRUP_DIA: float = 12.0
"""Default stirrup bar diameter (mm)."""

_DEFAULT_STIRRUP_LEGS: int = 4
"""Default number of stirrup legs across the pier cap width."""


# ---------------------------------------------------------------------------
# Result data class
# ---------------------------------------------------------------------------

@dataclass
class PierCapDesignResult:
    """Comprehensive results of the pier cap design.

    All dimensions are in **mm**, forces in **kN**, moments in **kN.m**,
    areas in **mm2**, stresses in **MPa**, and crack widths in **mm**
    unless stated otherwise.
    """

    # -- Section properties ------------------------------------------------
    width: float
    """Width of pier cap in traffic direction (mm)."""

    depth_at_pier: float
    """Overall depth at the pier face (mm)."""

    depth_at_curtailment: float
    """Overall depth at the curtailment / outer-bearing section (mm)."""

    cover: float
    """Clear cover to main reinforcement (mm)."""

    d_eff_pier: float
    """Effective depth at the pier face (mm)."""

    d_eff_curt: float
    """Effective depth at the curtailment section (mm)."""

    # -- ULS Flexure -------------------------------------------------------
    Mu_pier: float
    """Applied ULS moment at the pier face (kN.m)."""

    Mu_curt: float
    """Applied ULS moment at the curtailment section (kN.m)."""

    Mu_capacity_pier: float
    """Moment capacity at the pier face (kN.m)."""

    Mu_capacity_curt: float
    """Moment capacity at the curtailment section (kN.m)."""

    xu_d_pier: float
    """Neutral-axis depth ratio xu/d at the pier face."""

    xu_d_curt: float
    """Neutral-axis depth ratio xu/d at the curtailment section."""

    Ast_reqd_pier: float
    """Required tensile steel area at the pier face (mm2)."""

    Ast_reqd_curt: float
    """Required tensile steel area at the curtailment section (mm2)."""

    Ast_provided_pier: float
    """Provided tensile steel area at the pier face (mm2)."""

    Ast_provided_curt: float
    """Provided tensile steel area at the curtailment section (mm2)."""

    flexure_util_pier: float
    """Flexure utilisation ratio Mu / Mu_capacity at the pier face."""

    flexure_util_curt: float
    """Flexure utilisation ratio Mu / Mu_capacity at the curtailment."""

    # -- ULS Shear ---------------------------------------------------------
    VEd_pier: float
    """Applied ULS shear at the pier face (kN)."""

    VEd_curt: float
    """Applied ULS shear at the curtailment section (kN)."""

    VRdc_pier: float
    """Concrete shear resistance at the pier face (kN)."""

    VRdc_curt: float
    """Concrete shear resistance at the curtailment section (kN)."""

    shear_reinf_reqd: bool
    """True if shear reinforcement is needed at any section."""

    Asw_s_reqd: float
    """Required shear reinforcement ratio Asw/s (mm2/mm)."""

    Asw_s_provided: float
    """Provided shear reinforcement ratio Asw/s (mm2/mm)."""

    VRds: float
    """Shear resistance provided by stirrups (kN)."""

    VRd_max: float
    """Maximum shear resistance limited by strut crushing (kN)."""

    shear_util: float
    """Shear utilisation ratio max(VEd) / max(VRdc, VRds)."""

    # -- SLS Crack Width ---------------------------------------------------
    crack_width: float
    """Calculated crack width at the pier face under SLS (mm)."""

    crack_width_limit: float
    """Permissible crack width for the given exposure (mm)."""

    crack_ok: bool
    """True if calculated crack width <= limit."""

    # -- Minimum reinforcement ---------------------------------------------
    As_min: float
    """Minimum reinforcement area per IRC 112 Cl. 16.5.1.1 (mm2)."""

    min_reinf_ok: bool
    """True if provided steel >= As_min at all sections."""

    # -- Overall status ----------------------------------------------------
    status: str
    """``'OK'`` if all checks pass, otherwise ``'NOT OK'``."""


# ---------------------------------------------------------------------------
# Internal helpers -- section depth along cantilever
# ---------------------------------------------------------------------------

def _depth_at_x(
    x: float,
    L_cant: float,
    depth_max: float,
    depth_min: float,
) -> float:
    """Return the overall depth at distance *x* from the cantilever tip.

    The taper is linear:
        d(x) = depth_min + (depth_max - depth_min) * x / L_cant

    Parameters
    ----------
    x : float
        Distance measured from the cantilever tip towards the pier face.
        ``x = 0`` at the tip, ``x = L_cant`` at the pier face.
    L_cant : float
        Total cantilever length (tip to pier face).
    depth_max : float
        Depth at pier face (x = L_cant).
    depth_min : float
        Depth at tip (x = 0).

    Returns
    -------
    float
        Overall depth at position *x*.
    """
    if L_cant <= 0.0:
        return depth_max
    ratio = min(max(x / L_cant, 0.0), 1.0)
    return depth_min + (depth_max - depth_min) * ratio


# ---------------------------------------------------------------------------
# Internal helpers -- flexure design (IRC 112 Cl. 8.2)
# ---------------------------------------------------------------------------

def _solve_xu(
    Mu_Nmm: float,
    fck: float,
    b: float,
    d: float,
) -> float:
    """Solve for neutral-axis depth xu from the moment equilibrium.

    For the IRC 112 rectangular stress-block:
        Mu = 0.36 * fck * b * xu * (d - 0.42 * xu)

    Rearranging as a quadratic in xu:
        0.36 * fck * b * 0.42 * xu^2 - 0.36 * fck * b * d * xu + Mu = 0
        a_coeff * xu^2 - b_coeff * xu + Mu = 0

    Parameters
    ----------
    Mu_Nmm : float
        Ultimate moment in N.mm.
    fck : float
        Characteristic compressive strength in MPa.
    b : float
        Section width in mm.
    d : float
        Effective depth in mm.

    Returns
    -------
    float
        Neutral-axis depth xu in mm.  Returns the smaller root (under-
        reinforced solution).  If no real root exists the section is
        inadequate and the function returns xu_lim = 0.46 * d.
    """
    # Quadratic: 0.1512 * fck * b * xu^2  -  0.36 * fck * b * d * xu  +  Mu = 0
    a_coeff = 0.36 * fck * b * 0.42   # = 0.1512 * fck * b
    b_coeff = 0.36 * fck * b * d
    c_coeff = Mu_Nmm

    discriminant = b_coeff ** 2 - 4.0 * a_coeff * c_coeff
    if discriminant < 0.0:
        # Section is inadequate for singly-reinforced design
        return _XU_D_LIMIT_FE500 * d

    xu = (b_coeff - math.sqrt(discriminant)) / (2.0 * a_coeff)
    return max(xu, 0.0)


def _moment_capacity(
    fck: float,
    b: float,
    xu: float,
    d: float,
) -> float:
    """Return the moment capacity Mu in N.mm for a given xu.

    Mu = 0.36 * fck * b * xu * (d - 0.42 * xu)
    """
    return 0.36 * fck * b * xu * (d - 0.42 * xu)


def _ast_required(
    Mu_Nmm: float,
    fyk: float,
    xu: float,
    d: float,
) -> float:
    """Required tensile steel area for the given moment.

    Ast = Mu / (0.87 * fyk * (d - 0.42 * xu))

    Parameters
    ----------
    Mu_Nmm : float
        Ultimate moment in N.mm.
    fyk : float
        Characteristic yield strength of steel in MPa.
    xu : float
        Neutral-axis depth in mm.
    d : float
        Effective depth in mm.

    Returns
    -------
    float
        Steel area in mm2.
    """
    lever = d - 0.42 * xu
    if lever <= 0.0:
        lever = 0.1 * d  # fallback -- section inadequate
    return Mu_Nmm / (0.87 * fyk * lever)


def _provide_steel(
    Ast_reqd: float,
    bar_dia: float = _DEFAULT_BAR_DIA,
) -> tuple[float, int]:
    """Round up the required steel area to whole bars.

    Parameters
    ----------
    Ast_reqd : float
        Required steel area in mm2.
    bar_dia : float
        Bar diameter in mm.

    Returns
    -------
    tuple[float, int]
        (area_provided_mm2, number_of_bars)
    """
    A_bar = math.pi * bar_dia ** 2 / 4.0
    n_bars = max(2, math.ceil(Ast_reqd / A_bar))
    return n_bars * A_bar, n_bars


# ---------------------------------------------------------------------------
# Internal helpers -- shear design (IRC 112 Cl. 10.3)
# ---------------------------------------------------------------------------

def _VRdc(
    fck: float,
    bw: float,
    d: float,
    Ast: float,
    sigma_cp: float = 0.0,
) -> float:
    """Concrete shear resistance VRd,c per IRC 112 Cl. 10.3.2.

    Parameters
    ----------
    fck : float
        Characteristic compressive strength (MPa).
    bw : float
        Web width (mm).
    d : float
        Effective depth (mm).
    Ast : float
        Area of tensile reinforcement (mm2).
    sigma_cp : float
        Compressive stress due to axial load NEd/Ac (MPa), typically 0
        for pier caps.

    Returns
    -------
    float
        VRd,c in **N**.
    """
    k = min(1.0 + math.sqrt(200.0 / d), 2.0)
    rho_l = min(Ast / (bw * d), 0.02)

    VRdc_main = (
        0.12 * k * (80.0 * rho_l * fck) ** (1.0 / 3.0)
        + 0.15 * sigma_cp
    ) * bw * d

    VRdc_min = (
        0.031 * k ** 1.5 * fck ** 0.5
        + 0.15 * sigma_cp
    ) * bw * d

    return max(VRdc_main, VRdc_min)


def _shear_reinforcement(
    VEd_N: float,
    VRdc_N: float,
    fck: float,
    fyk: float,
    fywd: float,
    bw: float,
    d: float,
) -> tuple[float, float, float, float]:
    """Design shear reinforcement if VEd > VRdc.

    Parameters
    ----------
    VEd_N : float
        Design shear force in N.
    VRdc_N : float
        Concrete shear resistance in N.
    fck : float
        Characteristic concrete strength (MPa).
    fyk : float
        Characteristic steel yield strength (MPa).
    fywd : float
        Design yield strength for shear reinforcement (MPa).
        Per IRC 112: fywd = 0.8 * min(fyk, 500).
    bw : float
        Web width (mm).
    d : float
        Effective depth (mm).

    Returns
    -------
    tuple[float, float, float, float]
        (Asw_s_reqd, Asw_s_min, VRds_N, VRd_max_N)

        * Asw_s_reqd  -- required Asw/s (mm2/mm)
        * Asw_s_min   -- minimum Asw/s (mm2/mm)
        * VRds_N      -- shear resistance from stirrups (N)
        * VRd_max_N   -- maximum shear limited by strut crushing (N)
    """
    z = 0.9 * d  # lever arm

    # Strut angle -- start with cot(theta) = 2.5 (theta ~ 21.8 deg)
    cot_theta = 2.5
    tan_theta = 1.0 / cot_theta

    # VRd,max (crushing of compression strut)
    nu = 0.6 * (1.0 - fck / 310.0)
    fcd = 0.67 * fck / 1.5  # alpha_cc * fck / gamma_c
    VRd_max_N = bw * z * nu * fcd / (cot_theta + tan_theta)

    # If VEd exceeds VRd,max at cot(theta)=2.5, increase theta
    if VEd_N > VRd_max_N:
        # Try cot(theta) = 1.0 (theta = 45 deg) for maximum VRd,max
        cot_theta = 1.0
        tan_theta = 1.0
        VRd_max_N = bw * z * nu * fcd / (cot_theta + tan_theta)

    # Required Asw/s from VRd,s = (Asw/s) * z * fywd * cot(theta)
    if VEd_N > VRdc_N:
        Asw_s_reqd = VEd_N / (z * fywd * cot_theta)
    else:
        Asw_s_reqd = 0.0

    # Minimum shear reinforcement
    Asw_s_min = 0.072 * bw * math.sqrt(fck) / fyk

    # Governing Asw/s
    Asw_s_reqd = max(Asw_s_reqd, Asw_s_min)

    # Capacity provided (at governing Asw/s)
    VRds_N = Asw_s_reqd * z * fywd * cot_theta

    return Asw_s_reqd, Asw_s_min, VRds_N, VRd_max_N


def _provide_stirrups(
    Asw_s_reqd: float,
    stirrup_dia: float = _DEFAULT_STIRRUP_DIA,
    n_legs: int = _DEFAULT_STIRRUP_LEGS,
) -> tuple[float, float]:
    """Select stirrup spacing and return provided Asw/s.

    Parameters
    ----------
    Asw_s_reqd : float
        Required Asw/s (mm2/mm).
    stirrup_dia : float
        Stirrup bar diameter (mm).
    n_legs : int
        Number of stirrup legs.

    Returns
    -------
    tuple[float, float]
        (Asw_s_provided, spacing_mm)
    """
    Asw_one = n_legs * math.pi * stirrup_dia ** 2 / 4.0  # mm2 per set

    if Asw_s_reqd <= 0.0:
        spacing = 200.0  # nominal
    else:
        spacing = Asw_one / Asw_s_reqd

    # Round spacing down to nearest 25 mm
    spacing = max(25.0, 25.0 * math.floor(spacing / 25.0))

    # Cap spacing per code: <= min(0.75*d, 300 mm) -- approximate
    spacing = min(spacing, 300.0)

    Asw_s_provided = Asw_one / spacing
    return Asw_s_provided, spacing


# ---------------------------------------------------------------------------
# Internal helpers -- SLS crack width (IRC 112 Cl. 12.3.4)
# ---------------------------------------------------------------------------

def _crack_width(
    M_sls_Nmm: float,
    b: float,
    h: float,
    d: float,
    Ast: float,
    cover: float,
    bar_dia: float,
    fctm: float,
    Es_MPa: float,
    Ecm_GPa: float,
    *,
    alpha_e_override: float | None = None,
) -> float:
    """Calculate crack width wk per IRC 112 Cl. 12.3.4.

    Parameters
    ----------
    M_sls_Nmm : float
        SLS service moment in N.mm (quasi-permanent combination).
    b : float
        Section width (mm).
    h : float
        Overall section depth (mm).
    d : float
        Effective depth (mm).
    Ast : float
        Provided tensile steel area (mm2).
    cover : float
        Clear cover to outermost bar (mm).
    bar_dia : float
        Main bar diameter (mm).
    fctm : float
        Mean tensile strength of concrete (MPa).
    Es_MPa : float
        Steel elastic modulus (MPa).
    Ecm_GPa : float
        Concrete secant modulus (GPa).
    alpha_e_override : float or None
        If provided, use this modular ratio instead of Es/Ecm.
        Typically the long-term ratio Es / Ec_eff where
        Ec_eff = Ecm / (1 + phi_creep).

    Returns
    -------
    float
        Calculated crack width wk (mm).
    """
    if alpha_e_override is not None:
        alpha_e = alpha_e_override
    else:
        Ecm_MPa = Ecm_GPa * 1000.0
        alpha_e = Es_MPa / Ecm_MPa  # short-term modular ratio

    # ---- Steel stress under SLS ----
    # Cracked-section neutral axis depth by transformed-section method:
    #   0.5 * b * x^2 = alpha_e * Ast * (d - x)
    #   0.5 * b * x^2 + alpha_e * Ast * x - alpha_e * Ast * d = 0
    a_q = 0.5 * b
    b_q = alpha_e * Ast
    c_q = -alpha_e * Ast * d

    disc = b_q ** 2 - 4.0 * a_q * c_q
    if disc < 0.0:
        x_cr = 0.4 * d  # fallback
    else:
        x_cr = (-b_q + math.sqrt(disc)) / (2.0 * a_q)

    # Lever arm
    z_cr = d - x_cr / 3.0

    # Steel stress
    if Ast > 0.0 and z_cr > 0.0:
        sigma_s = M_sls_Nmm / (Ast * z_cr)
    else:
        sigma_s = 0.0

    # ---- Effective tension area ----
    h_c_eff = min(2.5 * (h - d), h / 2.0, (h - x_cr) / 3.0)
    h_c_eff = max(h_c_eff, 0.0)
    Ac_eff = b * h_c_eff
    rho_p_eff = Ast / Ac_eff if Ac_eff > 0.0 else 0.01

    # ---- Strain difference (epsilon_sm - epsilon_cm) ----
    kt = 0.5  # long-term loading
    if rho_p_eff > 0.0 and Es_MPa > 0.0:
        eps_diff = (
            sigma_s / Es_MPa
            - kt * fctm / rho_p_eff * (1.0 + alpha_e * rho_p_eff) / Es_MPa
        )
        eps_min = 0.6 * sigma_s / Es_MPa
        eps_diff = max(eps_diff, eps_min)
    else:
        eps_diff = 0.0

    # ---- Maximum crack spacing sr_max ----
    k1 = 0.8   # deformed bars
    k2 = 0.5   # bending
    if rho_p_eff > 0.0:
        sr_max = 3.4 * cover + 0.425 * k1 * k2 * bar_dia / rho_p_eff
    else:
        sr_max = 3.4 * cover + 0.425 * k1 * k2 * bar_dia / 0.01

    # ---- Crack width ----
    wk = sr_max * eps_diff
    return max(wk, 0.0)


# ---------------------------------------------------------------------------
# Internal helpers -- minimum reinforcement
# ---------------------------------------------------------------------------

def _min_reinforcement(
    fctm: float,
    fyk: float,
    bt: float,
    d: float,
) -> float:
    """Minimum tensile reinforcement per IRC 112 Cl. 16.5.1.1.

    As_min = max(0.26 * fctm / fyk * bt * d,  0.0013 * bt * d)

    Parameters
    ----------
    fctm : float
        Mean tensile strength of concrete (MPa).
    fyk : float
        Characteristic yield strength of steel (MPa).
    bt : float
        Mean width of the tension zone (mm).
    d : float
        Effective depth (mm).

    Returns
    -------
    float
        Minimum steel area (mm2).
    """
    return max(
        0.26 * fctm / fyk * bt * d,
        0.0013 * bt * d,
    )


# ---------------------------------------------------------------------------
# Internal helpers -- design forces at pier-cap sections
# ---------------------------------------------------------------------------

def _extract_pier_cap_forces(
    config: dict[str, Any],
    geometry: GeometryResults,
    combinations: CombinationResults,
) -> dict[str, Any]:
    """Extract design forces (BM, SF) at the pier-cap critical sections.

    The pier cap is modelled as a cantilever beam projecting from both
    faces of the pier.  The loads acting on each cantilever are:

    * Bearing reactions from the superstructure (DL + SIDL + LL) applied
      at the transverse bearing locations.
    * Self-weight of the pier cap itself as a distributed load.

    The critical cantilever is the one carrying the most load (or both
    are designed for the envelope).

    Returns a dict with keys:
        Mu_uls_pier, Mu_uls_curt, VEd_uls_pier, VEd_uls_curt,
        M_sls_pier (for crack width check),
        bearing_positions_from_pier_face (list of distances in mm).
    """
    pcap_cfg = config["pier_cap"]
    pier_cfg = config["pier"]
    mat_cfg = config["materials"]

    # -- Pier cap dimensions (convert m -> mm) --
    pcap_length_m = float(pcap_cfg["length_trans"])       # m
    pcap_width_m = float(pcap_cfg["width_long"])          # m
    depth_max_m = float(pcap_cfg["depth_max"])            # m
    depth_min_m = float(pcap_cfg["depth_min"])            # m
    density = float(mat_cfg.get("concrete_density", 25.0))  # kN/m3

    # Pier transverse dimension (to find cantilever length)
    if pier_cfg["type"] == "Circular":
        pier_dim_m = float(pier_cfg["diameter_bottom"])
    else:
        pier_dim_m = float(pier_cfg["diameter_top"])  # transverse

    L_cant_m = (pcap_length_m - pier_dim_m) / 2.0  # m  cantilever each side
    if L_cant_m <= 0.0:
        L_cant_m = 0.01  # safety fallback

    # -- Bearing positions relative to pier centre (transverse) --
    # Combine LHS and RHS bearing coordinates.  The transverse coordinate
    # is the second element [1] of each coordinate pair.  We work with
    # absolute distances from the pier centre.
    bearing_y_coords_m: list[float] = []
    for coord in geometry.bearing_coords_lhs:
        bearing_y_coords_m.append(abs(coord[1]))
    for coord in geometry.bearing_coords_rhs:
        bearing_y_coords_m.append(abs(coord[1]))

    # Remove duplicates (LHS and RHS are symmetric) and sort
    unique_y_m = sorted(set(bearing_y_coords_m))

    # Distance from pier face (not pier centre)
    pier_face_offset_m = pier_dim_m / 2.0
    bearing_dist_from_face_m = [
        y - pier_face_offset_m for y in unique_y_m if y > pier_face_offset_m
    ]

    # Curtailment section: at the outermost bearing
    if bearing_dist_from_face_m:
        curt_dist_m = max(bearing_dist_from_face_m)
    else:
        curt_dist_m = L_cant_m * 0.5  # fallback

    # -- Self-weight UDL of pier cap (trapezoidal, approximate as average) --
    avg_depth_m = (depth_max_m + depth_min_m) / 2.0
    sw_udl_kN_per_m = pcap_width_m * avg_depth_m * density  # kN/m along transverse

    # -- Bearing reactions from governing combinations --
    # ULS: use the governing ULS basic combination for maximum P
    # SLS: use the governing SLS quasi-permanent for crack width
    #      (IRC 112 Cl. 12.3.4 requires quasi-permanent combination)
    gov_uls = combinations.governing_uls_basic
    gov_sls = combinations.governing_sls_quasi

    # Per-bearing vertical reaction under ULS.
    # The total vertical load P at pier base includes DL + SIDL + LL + ...
    # For the pier cap design we need the *bearing reactions* only (loads
    # transmitted from the superstructure).  Approximate these as:
    #   R_per_bearing = (P_total - self_weight_pier - self_weight_piercap) / n_bearings
    # where P_total is from the governing combination.

    n_brg_total = (len(geometry.bearing_coords_lhs)
                   + len(geometry.bearing_coords_rhs))
    if n_brg_total == 0:
        n_brg_total = 1  # safety

    # Pier and piercap self-weight (unfactored) for subtraction
    wt_piercap_kN = geometry.piercap_section.width * geometry.piercap_section.depth * pcap_length_m * density
    # Actually piercap_section stores width in m and depth in m already
    # from geometry module:  area = width_long * depth_max, so:
    wt_piercap_kN = (geometry.piercap_section.area  # m2
                     * pcap_length_m * density)  # area is width*depth, need *length

    # Correction: piercap_section.area = width_long * depth_max (cross-section area).
    # Weight = area * length_trans * density -- but that double-counts.
    # Actually area in PierCapSection is width_long * depth_max, which is the
    # cross-section area when looking along the transverse axis.
    # Weight = cross_section_area * length_trans * density
    # But that assumes constant depth.  Use the average depth for trapezoidal:
    wt_piercap_kN = pcap_width_m * avg_depth_m * pcap_length_m * density

    pier_area_m2 = geometry.pier_section.area
    pier_height_m = geometry.pier_height
    wt_pier_kN = pier_area_m2 * pier_height_m * density

    # ULS bearing reaction per bearing
    P_uls_total = gov_uls.forces_pier_base.P
    # Subtract pier and piercap self-weight (factored approx at gamma=1.35)
    P_super_uls = P_uls_total - 1.35 * (wt_pier_kN + wt_piercap_kN)
    P_super_uls = max(P_super_uls, 0.0)
    R_per_bearing_uls = P_super_uls / n_brg_total

    # SLS bearing reaction per bearing
    P_sls_total = gov_sls.forces_pier_base.P
    P_super_sls = P_sls_total - 1.0 * (wt_pier_kN + wt_piercap_kN)
    P_super_sls = max(P_super_sls, 0.0)
    R_per_bearing_sls = P_super_sls / n_brg_total

    # -- Number of bearings per cantilever side --
    # Bearings are placed symmetrically.  Count bearings on one transverse
    # side of the pier (i.e. those beyond the pier face).
    # From the bearing coordinates, bearings with |y| > pier_face_offset
    # are on the cantilever.  Each unique |y| position typically has two
    # bearings (one LHS-of-pier, one RHS-of-pier in the longitudinal sense)
    # but they all load the *same* transverse cantilever.
    # Count bearings on ONE cantilever side only (positive y).
    # Bearings at +y and -y are on opposite cantilever sides; using abs(y)
    # would double-count.  We pick bearings with y > 0 (one side).
    n_brg_per_cant_side: dict[float, int] = {}
    for coord in geometry.bearing_coords_lhs + geometry.bearing_coords_rhs:
        y = coord[1]
        if y > 0.001:  # positive cantilever side only
            y_rounded = round(y, 3)
            if y_rounded > pier_face_offset_m - 0.001:
                n_brg_per_cant_side[y_rounded] = (
                    n_brg_per_cant_side.get(y_rounded, 0) + 1
                )

    # -- Moment and shear at pier face (ULS) --
    # M_pier = sum(R_i * e_i) + sw_moment
    # where e_i = distance of bearing from pier face, R_i = reaction at that bearing
    Mu_bearings_pier = 0.0
    VEd_bearings_pier = 0.0

    for y_abs_round, count in n_brg_per_cant_side.items():
        e_m = y_abs_round - pier_face_offset_m
        if e_m < 0.0:
            continue
        R_total_at_y = count * R_per_bearing_uls  # total reaction at this y-row
        # This acts on *one* cantilever side (we design one cantilever).
        # Each side of the pier has bearings; the count already includes both
        # longitudinal rows for this transverse position.
        # However, for the transverse cantilever, both LHS and RHS bearings
        # at the same |y| contribute to the *same* cantilever.
        Mu_bearings_pier += R_total_at_y * e_m  # kN * m = kN.m
        VEd_bearings_pier += R_total_at_y       # kN

    # Self-weight moment at pier face: UDL on cantilever
    # Factored at gamma = 1.35 for ULS
    sw_M_pier = 1.35 * sw_udl_kN_per_m * L_cant_m ** 2 / 2.0  # kN.m
    sw_V_pier = 1.35 * sw_udl_kN_per_m * L_cant_m              # kN

    Mu_uls_pier = Mu_bearings_pier + sw_M_pier
    VEd_uls_pier = VEd_bearings_pier + sw_V_pier

    # -- Moment and shear at curtailment section (ULS) --
    # Only bearings *beyond* the curtailment section contribute.
    Mu_bearings_curt = 0.0
    VEd_bearings_curt = 0.0

    for y_abs_round, count in n_brg_per_cant_side.items():
        e_from_face_m = y_abs_round - pier_face_offset_m
        if e_from_face_m < 0.0:
            continue
        e_from_curt_m = e_from_face_m - curt_dist_m
        if e_from_curt_m < 0.0:
            continue
        R_total_at_y = count * R_per_bearing_uls
        Mu_bearings_curt += R_total_at_y * e_from_curt_m
        VEd_bearings_curt += R_total_at_y

    # Self-weight beyond curtailment
    L_beyond_curt_m = L_cant_m - curt_dist_m
    if L_beyond_curt_m < 0.0:
        L_beyond_curt_m = 0.0

    sw_M_curt = 1.35 * sw_udl_kN_per_m * L_beyond_curt_m ** 2 / 2.0
    sw_V_curt = 1.35 * sw_udl_kN_per_m * L_beyond_curt_m

    Mu_uls_curt = Mu_bearings_curt + sw_M_curt
    VEd_uls_curt = VEd_bearings_curt + sw_V_curt

    # -- SLS moment at pier face (for crack width) --
    M_sls_bearings_pier = 0.0
    for y_abs_round, count in n_brg_per_cant_side.items():
        e_m = y_abs_round - pier_face_offset_m
        if e_m < 0.0:
            continue
        R_sls = count * R_per_bearing_sls
        M_sls_bearings_pier += R_sls * e_m

    sw_M_pier_sls = 1.0 * sw_udl_kN_per_m * L_cant_m ** 2 / 2.0
    M_sls_pier = M_sls_bearings_pier + sw_M_pier_sls

    return {
        "Mu_uls_pier": Mu_uls_pier,
        "Mu_uls_curt": Mu_uls_curt,
        "VEd_uls_pier": VEd_uls_pier,
        "VEd_uls_curt": VEd_uls_curt,
        "M_sls_pier": M_sls_pier,
        "L_cant_m": L_cant_m,
        "curt_dist_m": curt_dist_m,
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def design_pier_cap(
    config: dict[str, Any],
    geometry: GeometryResults,
    combinations: CombinationResults,
    creep_result: Any = None,
) -> PierCapDesignResult:
    """Design the pier cap for flexure, shear, and crack width.

    This is the primary public entry point.  It reads pier cap geometry
    from *config* and *geometry*, extracts design forces from the
    governing load combinations, and performs all IRC 112:2020 checks.

    Parameters
    ----------
    config : dict
        Validated configuration dictionary (from ``input_parser``).
    geometry : GeometryResults
        Pre-computed geometry (from ``geometry.calculate_geometry``).
    combinations : CombinationResults
        Load combinations (from ``load_combinations.calculate_combinations``).
    creep_result : CreepResult or None
        Creep calculation results.  If provided, the long-term modular
        ratio ``Es / Ec_eff`` is used for crack width calculations.

    Returns
    -------
    PierCapDesignResult
        Comprehensive design results including reinforcement areas,
        capacity ratios, crack width, and overall pass/fail status.
    """
    irc = load_irc_tables()

    # ==================================================================
    # 1.  Material properties
    # ==================================================================
    mat_cfg = config["materials"]
    fck = float(mat_cfg["pier_cap"]["fck"])
    aggregate = str(mat_cfg["pier_cap"].get("aggregate", "Quartzite"))
    exposure = str(mat_cfg["exposure"])
    fyk = float(mat_cfg["steel"]["fyk"])
    Es_input = float(mat_cfg["steel"]["Es"])  # may be in MPa (200000)

    # Normalise Es to GPa for get_steel_properties (expects GPa)
    Es_GPa = Es_input if Es_input < 1000.0 else Es_input / 1000.0

    concrete = get_concrete_properties(fck, aggregate=aggregate, exposure=exposure)
    steel = get_steel_properties(fyk=fyk, Es=Es_GPa)

    # ==================================================================
    # 2.  Pier cap geometry (all in mm)
    # ==================================================================
    pcap_cfg = config["pier_cap"]
    pier_cfg = config["pier"]

    width_mm = float(pcap_cfg["width_long"]) * 1000.0      # mm
    depth_max_mm = float(pcap_cfg["depth_max"]) * 1000.0    # mm
    depth_min_mm = float(pcap_cfg["depth_min"]) * 1000.0    # mm
    pcap_length_mm = float(pcap_cfg["length_trans"]) * 1000.0  # mm

    if pier_cfg["type"] == "Circular":
        pier_dim_mm = float(pier_cfg["diameter_bottom"]) * 1000.0
    else:
        pier_dim_mm = float(pier_cfg["diameter_top"]) * 1000.0

    L_cant_mm = (pcap_length_mm - pier_dim_mm) / 2.0

    # Cover
    if concrete.min_cover is not None:
        cover = concrete.min_cover
    else:
        cover = 45.0  # default for severe exposure

    bar_dia = _DEFAULT_BAR_DIA
    stirrup_dia = _DEFAULT_STIRRUP_DIA

    # Effective depths
    d_eff_pier = depth_max_mm - cover - stirrup_dia - bar_dia / 2.0
    d_eff_curt_depth_mm = _depth_at_x(
        x=L_cant_mm,  # pier face is at x = L_cant from tip
        L_cant=L_cant_mm,
        depth_max=depth_max_mm,
        depth_min=depth_min_mm,
    )
    # Actually we need depth at the curtailment section.
    # Curtailment is at some distance from the pier face along the cantilever.
    # Distance from tip to curtailment = L_cant - curt_dist

    # ==================================================================
    # 3.  Design forces
    # ==================================================================
    forces = _extract_pier_cap_forces(config, geometry, combinations)
    Mu_uls_pier = forces["Mu_uls_pier"]        # kN.m
    Mu_uls_curt = forces["Mu_uls_curt"]        # kN.m
    VEd_uls_pier = forces["VEd_uls_pier"]      # kN
    VEd_uls_curt = forces["VEd_uls_curt"]      # kN
    M_sls_pier = forces["M_sls_pier"]          # kN.m
    curt_dist_m = forces["curt_dist_m"]        # m
    L_cant_m_val = forces["L_cant_m"]          # m

    # Depth at curtailment section
    curt_dist_mm = curt_dist_m * 1000.0
    dist_from_tip_mm = L_cant_mm - curt_dist_mm
    depth_at_curt_mm = _depth_at_x(
        x=dist_from_tip_mm,
        L_cant=L_cant_mm,
        depth_max=depth_max_mm,
        depth_min=depth_min_mm,
    )
    d_eff_curt = depth_at_curt_mm - cover - stirrup_dia - bar_dia / 2.0

    # Ensure positive effective depths
    d_eff_pier = max(d_eff_pier, 100.0)
    d_eff_curt = max(d_eff_curt, 100.0)

    # ==================================================================
    # 4.  Flexure design (IRC 112 Cl. 8.2)
    # ==================================================================
    # Convert moments to N.mm:  kN.m -> N.mm = * 1e6
    Mu_pier_Nmm = Mu_uls_pier * 1.0e6
    Mu_curt_Nmm = Mu_uls_curt * 1.0e6

    # --- Pier face ---
    xu_pier = _solve_xu(Mu_pier_Nmm, fck, width_mm, d_eff_pier)
    xu_d_pier = xu_pier / d_eff_pier
    if xu_d_pier > _XU_D_LIMIT_FE500:
        xu_pier = _XU_D_LIMIT_FE500 * d_eff_pier
        xu_d_pier = _XU_D_LIMIT_FE500

    Mu_cap_pier_Nmm = _moment_capacity(fck, width_mm, xu_pier, d_eff_pier)
    Ast_reqd_pier = _ast_required(Mu_pier_Nmm, fyk, xu_pier, d_eff_pier)
    # Enforce minimum reinforcement per IRC 112 Cl. 16.5.1.1
    As_min_pier_early = _min_reinforcement(concrete.fctm, fyk, width_mm, d_eff_pier)
    Ast_reqd_pier = max(Ast_reqd_pier, As_min_pier_early)
    Ast_prov_pier, _ = _provide_steel(Ast_reqd_pier, bar_dia)

    # Recalculate capacity with provided steel
    xu_pier_prov = _solve_xu(
        _moment_capacity(fck, width_mm, _XU_D_LIMIT_FE500 * d_eff_pier, d_eff_pier),
        fck, width_mm, d_eff_pier,
    )
    # Simpler: compute capacity from provided Ast
    # From equilibrium: 0.87 * fyk * Ast = 0.36 * fck * b * xu
    xu_from_Ast_pier = 0.87 * fyk * Ast_prov_pier / (0.36 * fck * width_mm)
    xu_from_Ast_pier = min(xu_from_Ast_pier, _XU_D_LIMIT_FE500 * d_eff_pier)
    Mu_cap_pier_prov_Nmm = _moment_capacity(
        fck, width_mm, xu_from_Ast_pier, d_eff_pier,
    )

    flexure_util_pier = Mu_pier_Nmm / Mu_cap_pier_prov_Nmm if Mu_cap_pier_prov_Nmm > 0 else 999.0

    # --- Curtailment section ---
    xu_curt = _solve_xu(Mu_curt_Nmm, fck, width_mm, d_eff_curt)
    xu_d_curt = xu_curt / d_eff_curt
    if xu_d_curt > _XU_D_LIMIT_FE500:
        xu_curt = _XU_D_LIMIT_FE500 * d_eff_curt
        xu_d_curt = _XU_D_LIMIT_FE500

    Ast_reqd_curt = _ast_required(Mu_curt_Nmm, fyk, xu_curt, d_eff_curt)
    # Enforce minimum reinforcement per IRC 112 Cl. 16.5.1.1
    As_min_curt_early = _min_reinforcement(concrete.fctm, fyk, width_mm, d_eff_curt)
    Ast_reqd_curt = max(Ast_reqd_curt, As_min_curt_early)
    Ast_prov_curt, _ = _provide_steel(Ast_reqd_curt, bar_dia)

    xu_from_Ast_curt = 0.87 * fyk * Ast_prov_curt / (0.36 * fck * width_mm)
    xu_from_Ast_curt = min(xu_from_Ast_curt, _XU_D_LIMIT_FE500 * d_eff_curt)
    Mu_cap_curt_prov_Nmm = _moment_capacity(
        fck, width_mm, xu_from_Ast_curt, d_eff_curt,
    )

    flexure_util_curt = Mu_curt_Nmm / Mu_cap_curt_prov_Nmm if Mu_cap_curt_prov_Nmm > 0 else 999.0

    # ==================================================================
    # 5.  Shear design (IRC 112 Cl. 10.3)
    # ==================================================================
    # Use provided Ast for concrete shear capacity
    VRdc_pier_N = _VRdc(fck, width_mm, d_eff_pier, Ast_prov_pier)
    VRdc_curt_N = _VRdc(fck, width_mm, d_eff_curt, Ast_prov_curt)

    VEd_pier_N = VEd_uls_pier * 1000.0   # kN -> N
    VEd_curt_N = VEd_uls_curt * 1000.0

    # Governing shear section
    VEd_gov_N = max(VEd_pier_N, VEd_curt_N)
    VRdc_gov_N = VRdc_pier_N  # at pier face (usually more critical due to larger d)

    shear_reinf_reqd = VEd_gov_N > VRdc_gov_N

    fyd_MPa = steel.fyd  # already in MPa from get_steel_properties
    # IRC 112: fywd = 0.8 * fywk, where fywk = min(fyk, 500) for shear reinforcement
    fywd_MPa = 0.8 * min(fyk, 500.0)

    if shear_reinf_reqd:
        Asw_s_reqd, Asw_s_min, VRds_N, VRd_max_N = _shear_reinforcement(
            VEd_gov_N, VRdc_gov_N,
            fck, fyk, fywd_MPa, width_mm, d_eff_pier,
        )
    else:
        # Still provide minimum shear reinforcement
        Asw_s_min = 0.072 * width_mm * math.sqrt(fck) / fyk
        Asw_s_reqd = Asw_s_min
        z = 0.9 * d_eff_pier
        cot_theta = 2.5
        VRds_N = Asw_s_min * z * fywd_MPa * cot_theta
        nu = 0.6 * (1.0 - fck / 310.0)
        fcd_val = 0.67 * fck / 1.5
        VRd_max_N = width_mm * z * nu * fcd_val / (cot_theta + 1.0 / cot_theta)

    Asw_s_prov, stirrup_spacing = _provide_stirrups(
        Asw_s_reqd, stirrup_dia, _DEFAULT_STIRRUP_LEGS,
    )

    # Recalculate VRds with provided stirrups
    z = 0.9 * d_eff_pier
    cot_theta = 2.5
    VRds_prov_N = Asw_s_prov * z * fywd_MPa * cot_theta

    shear_capacity_N = max(VRdc_pier_N, VRds_prov_N) if shear_reinf_reqd else VRdc_pier_N
    shear_util = VEd_gov_N / shear_capacity_N if shear_capacity_N > 0 else 999.0

    # ==================================================================
    # 6.  Crack width (IRC 112 Cl. 12.3.4)
    #     Uses quasi-permanent SLS moment and long-term modular ratio
    # ==================================================================
    M_sls_Nmm = M_sls_pier * 1.0e6  # kN.m -> N.mm
    Es_MPa = Es_GPa * 1000.0

    # Long-term modular ratio: alpha_e = Es / Ec_eff
    # where Ec_eff = Ecm / (1 + phi_creep)
    alpha_e_lt = None
    if creep_result is not None and hasattr(creep_result, 'phi'):
        Ec_eff_MPa = concrete.Ecm * 1000.0 / (1.0 + creep_result.phi)
        alpha_e_lt = Es_MPa / Ec_eff_MPa

    # Crack width limit from IRC tables
    crack_limits: dict[str, float] = irc.get("crack_width_limits", {})
    wk_limit = float(crack_limits.get(exposure, 0.3))

    # Iteratively increase Ast_prov_pier until crack width is within limit.
    # In practice, designers increase reinforcement for crack control.
    _max_iter = 20
    for _iter in range(_max_iter):
        wk = _crack_width(
            M_sls_Nmm=M_sls_Nmm,
            b=width_mm,
            h=depth_max_mm,
            d=d_eff_pier,
            Ast=Ast_prov_pier,
            cover=cover,
            bar_dia=bar_dia,
            fctm=concrete.fctm,
            Es_MPa=Es_MPa,
            Ecm_GPa=concrete.Ecm,
            alpha_e_override=alpha_e_lt,
        )
        if wk <= wk_limit or M_sls_Nmm <= 0.0:
            break
        # Increase Ast by ~10% per iteration
        Ast_prov_pier *= 1.10
        Ast_prov_pier, _ = _provide_steel(Ast_prov_pier, bar_dia)

    # Recompute flexure capacity with final Ast
    xu_from_Ast_pier = 0.87 * fyk * Ast_prov_pier / (0.36 * fck * width_mm)
    xu_from_Ast_pier = min(xu_from_Ast_pier, _XU_D_LIMIT_FE500 * d_eff_pier)
    Mu_cap_pier_prov_Nmm = _moment_capacity(
        fck, width_mm, xu_from_Ast_pier, d_eff_pier,
    )
    flexure_util_pier = Mu_pier_Nmm / Mu_cap_pier_prov_Nmm if Mu_cap_pier_prov_Nmm > 0 else 999.0

    crack_ok = wk <= wk_limit

    # ==================================================================
    # 7.  Minimum reinforcement
    # ==================================================================
    As_min_pier = _min_reinforcement(concrete.fctm, fyk, width_mm, d_eff_pier)
    As_min_curt = _min_reinforcement(concrete.fctm, fyk, width_mm, d_eff_curt)
    As_min = max(As_min_pier, As_min_curt)

    min_reinf_ok = (
        Ast_prov_pier >= As_min_pier and Ast_prov_curt >= As_min_curt
    )

    # ==================================================================
    # 8.  Overall status
    # ==================================================================
    flexure_ok = flexure_util_pier <= 1.0 and flexure_util_curt <= 1.0
    shear_ok = shear_util <= 1.0
    all_ok = flexure_ok and shear_ok and crack_ok and min_reinf_ok
    status = "OK" if all_ok else "NOT OK"

    # ==================================================================
    # 9.  Assemble result
    # ==================================================================
    return PierCapDesignResult(
        # Section properties
        width=width_mm,
        depth_at_pier=depth_max_mm,
        depth_at_curtailment=depth_at_curt_mm,
        cover=cover,
        d_eff_pier=d_eff_pier,
        d_eff_curt=d_eff_curt,
        # ULS Flexure
        Mu_pier=Mu_uls_pier,
        Mu_curt=Mu_uls_curt,
        Mu_capacity_pier=Mu_cap_pier_prov_Nmm / 1.0e6,  # N.mm -> kN.m
        Mu_capacity_curt=Mu_cap_curt_prov_Nmm / 1.0e6,
        xu_d_pier=xu_from_Ast_pier / d_eff_pier,
        xu_d_curt=xu_from_Ast_curt / d_eff_curt,
        Ast_reqd_pier=Ast_reqd_pier,
        Ast_reqd_curt=Ast_reqd_curt,
        Ast_provided_pier=Ast_prov_pier,
        Ast_provided_curt=Ast_prov_curt,
        flexure_util_pier=flexure_util_pier,
        flexure_util_curt=flexure_util_curt,
        # ULS Shear
        VEd_pier=VEd_uls_pier,
        VEd_curt=VEd_uls_curt,
        VRdc_pier=VRdc_pier_N / 1000.0,  # N -> kN
        VRdc_curt=VRdc_curt_N / 1000.0,
        shear_reinf_reqd=shear_reinf_reqd,
        Asw_s_reqd=Asw_s_reqd,
        Asw_s_provided=Asw_s_prov,
        VRds=VRds_prov_N / 1000.0,
        VRd_max=VRd_max_N / 1000.0,
        shear_util=shear_util,
        # SLS Crack Width
        crack_width=wk,
        crack_width_limit=wk_limit,
        crack_ok=crack_ok,
        # Minimum reinforcement
        As_min=As_min,
        min_reinf_ok=min_reinf_ok,
        # Overall
        status=status,
    )


# ---------------------------------------------------------------------------
# Summary / reporting
# ---------------------------------------------------------------------------

def summarise_pier_cap_design(result: PierCapDesignResult) -> str:
    """Return a human-readable summary of the pier cap design.

    Useful for console output and report generation.

    Parameters
    ----------
    result : PierCapDesignResult
        Output from :func:`design_pier_cap`.

    Returns
    -------
    str
        Multi-line formatted summary string.
    """
    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("PIER CAP DESIGN SUMMARY (IRC 112:2020)")
    lines.append("=" * 72)

    lines.append("")
    lines.append("  SECTION PROPERTIES")
    lines.append(f"    Width (b)                : {result.width:>10.0f} mm")
    lines.append(f"    Depth at pier face       : {result.depth_at_pier:>10.0f} mm")
    lines.append(f"    Depth at curtailment     : {result.depth_at_curtailment:>10.0f} mm")
    lines.append(f"    Clear cover              : {result.cover:>10.0f} mm")
    lines.append(f"    d_eff at pier face       : {result.d_eff_pier:>10.1f} mm")
    lines.append(f"    d_eff at curtailment     : {result.d_eff_curt:>10.1f} mm")

    lines.append("")
    lines.append("  ULS FLEXURE (IRC 112 Cl. 8.2)")
    lines.append("  " + "-" * 68)
    lines.append(f"    {'':30s}  {'Pier Face':>14s}  {'Curtailment':>14s}")
    lines.append(f"    {'Mu applied (kN.m)':30s}  {result.Mu_pier:>14.1f}  {result.Mu_curt:>14.1f}")
    lines.append(f"    {'Mu capacity (kN.m)':30s}  {result.Mu_capacity_pier:>14.1f}  {result.Mu_capacity_curt:>14.1f}")
    lines.append(f"    {'xu/d':30s}  {result.xu_d_pier:>14.3f}  {result.xu_d_curt:>14.3f}")
    lines.append(f"    {'Ast required (mm2)':30s}  {result.Ast_reqd_pier:>14.0f}  {result.Ast_reqd_curt:>14.0f}")
    lines.append(f"    {'Ast provided (mm2)':30s}  {result.Ast_provided_pier:>14.0f}  {result.Ast_provided_curt:>14.0f}")
    lines.append(f"    {'Utilisation':30s}  {result.flexure_util_pier:>14.3f}  {result.flexure_util_curt:>14.3f}")

    xu_limit_tag = f"(limit = {_XU_D_LIMIT_FE500:.2f})"
    flexure_status = "OK" if result.flexure_util_pier <= 1.0 and result.flexure_util_curt <= 1.0 else "NOT OK"
    lines.append(f"    xu/d limit {xu_limit_tag:30s}  Status: {flexure_status}")

    lines.append("")
    lines.append("  ULS SHEAR (IRC 112 Cl. 10.3)")
    lines.append("  " + "-" * 68)
    lines.append(f"    {'VEd at pier face (kN)':30s}  {result.VEd_pier:>14.1f}")
    lines.append(f"    {'VEd at curtailment (kN)':30s}  {result.VEd_curt:>14.1f}")
    lines.append(f"    {'VRd,c at pier face (kN)':30s}  {result.VRdc_pier:>14.1f}")
    lines.append(f"    {'VRd,c at curtailment (kN)':30s}  {result.VRdc_curt:>14.1f}")
    lines.append(f"    {'Shear reinf. required':30s}  {'Yes' if result.shear_reinf_reqd else 'No':>14s}")
    lines.append(f"    {'Asw/s required (mm2/mm)':30s}  {result.Asw_s_reqd:>14.3f}")
    lines.append(f"    {'Asw/s provided (mm2/mm)':30s}  {result.Asw_s_provided:>14.3f}")
    lines.append(f"    {'VRd,s (kN)':30s}  {result.VRds:>14.1f}")
    lines.append(f"    {'VRd,max (kN)':30s}  {result.VRd_max:>14.1f}")
    lines.append(f"    {'Shear utilisation':30s}  {result.shear_util:>14.3f}")

    shear_status = "OK" if result.shear_util <= 1.0 else "NOT OK"
    lines.append(f"    {'Status':30s}  {shear_status:>14s}")

    lines.append("")
    lines.append("  SLS CRACK WIDTH (IRC 112 Cl. 12.3.4)")
    lines.append("  " + "-" * 68)
    lines.append(f"    {'Calculated wk (mm)':30s}  {result.crack_width:>14.4f}")
    lines.append(f"    {'Limit wk (mm)':30s}  {result.crack_width_limit:>14.4f}")
    lines.append(f"    {'Status':30s}  {'OK' if result.crack_ok else 'NOT OK':>14s}")

    lines.append("")
    lines.append("  MINIMUM REINFORCEMENT (IRC 112 Cl. 16.5.1.1)")
    lines.append("  " + "-" * 68)
    lines.append(f"    {'As,min (mm2)':30s}  {result.As_min:>14.0f}")
    lines.append(f"    {'As,prov pier (mm2)':30s}  {result.Ast_provided_pier:>14.0f}")
    lines.append(f"    {'As,prov curt (mm2)':30s}  {result.Ast_provided_curt:>14.0f}")
    lines.append(f"    {'Status':30s}  {'OK' if result.min_reinf_ok else 'NOT OK':>14s}")

    lines.append("")
    lines.append("=" * 72)
    lines.append(f"  OVERALL STATUS : {result.status}")
    lines.append("=" * 72)

    return "\n".join(lines)

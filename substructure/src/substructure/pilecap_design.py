"""Pile cap design (flexure, shear, punching, crack width) per IRC 112:2020.

Designs the pile cap for flexure at the face of the pier, one-way shear at
distance d from the pier face, punching shear around the pier perimeter, and
punching shear around individual piles.

Design checks
-------------
1. **ULS Flexure** -- IRC 112 Cl. 8.2 (rectangular stress-block) in both
   longitudinal and transverse directions.
2. **ULS One-Way Shear** -- IRC 112 Cl. 10.3 (concrete capacity) in both
   longitudinal and transverse directions.
3. **ULS Punching Shear (Pier)** -- IRC 112 Cl. 10.4 (punching around pier).
4. **ULS Punching Shear (Pile)** -- IRC 112 Cl. 10.4 (punching around pile).
5. **SLS Crack Width** -- IRC 112 Cl. 12.3.4 (direct calculation).
6. **Minimum Reinforcement** -- IRC 112 Cl. 16.5.1.1.

Units convention
----------------
Forces in **kN**, moments in **kN.m**, dimensions in **mm** for
section-level calculations (converted from metres as needed).

Key references
--------------
* IRC 112:2020 -- Code of Practice for Concrete Road Bridges
* IRC 78:2014 -- Standard Specifications and Code of Practice for Road
  Bridges, Section VII: Foundations and Substructure
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

from .loads import ForceVector
from .geometry import GeometryResults, PileCoordinate, PierSection
from .load_combinations import CombinationResult, CombinationResults
from .materials import (
    get_concrete_properties,
    get_steel_properties,
    ConcreteProperties,
    SteelProperties,
)
from .utils import load_irc_tables
from .pile_capacity import PileCapacityResult, PileLoad, PileCombinationResult


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_XU_D_LIMIT_FE500: float = 0.46
"""Maximum permitted xu/d ratio for Fe 500 / Fe 550 steel (IRC 112)."""

_DEFAULT_BAR_DIA: float = 25.0
"""Default main bar diameter (mm) for pilecap reinforcement."""


# ---------------------------------------------------------------------------
# Result data class
# ---------------------------------------------------------------------------

@dataclass
class PilecapDesignResult:
    """Comprehensive results of the pile cap design.

    All dimensions are in **mm**, forces in **kN**, moments in **kN.m**,
    areas in **mm2**, stresses in **MPa**, and crack widths in **mm**
    unless stated otherwise.
    """

    # -- Geometry ----------------------------------------------------------
    length_long: float
    """Length in longitudinal direction (mm)."""

    width_trans: float
    """Width in transverse direction (mm)."""

    thickness: float
    """Thickness of pile cap (mm)."""

    cover: float
    """Clear cover to main reinforcement (mm)."""

    d_eff: float
    """Effective depth (mm)."""

    fck: float
    """Characteristic concrete strength (MPa)."""

    # -- Flexure (longitudinal) --------------------------------------------
    Mu_long: float
    """Applied ULS moment in longitudinal direction (kN.m)."""

    Mu_cap_long: float
    """Moment capacity in longitudinal direction (kN.m)."""

    Ast_reqd_long: float
    """Required steel area in longitudinal direction (mm2)."""

    Ast_prov_long: float
    """Provided steel area in longitudinal direction (mm2)."""

    flexure_util_long: float
    """Flexure utilisation ratio in longitudinal direction."""

    # -- Flexure (transverse) ----------------------------------------------
    Mu_trans: float
    """Applied ULS moment in transverse direction (kN.m)."""

    Mu_cap_trans: float
    """Moment capacity in transverse direction (kN.m)."""

    Ast_reqd_trans: float
    """Required steel area in transverse direction (mm2)."""

    Ast_prov_trans: float
    """Provided steel area in transverse direction (mm2)."""

    flexure_util_trans: float
    """Flexure utilisation ratio in transverse direction."""

    # -- One-way shear -----------------------------------------------------
    VEd_long: float
    """Applied ULS shear in longitudinal direction (kN)."""

    VRdc_long: float
    """Concrete shear resistance in longitudinal direction (kN)."""

    shear_util_long: float
    """Shear utilisation ratio in longitudinal direction."""

    VEd_trans: float
    """Applied ULS shear in transverse direction (kN)."""

    VRdc_trans: float
    """Concrete shear resistance in transverse direction (kN)."""

    shear_util_trans: float
    """Shear utilisation ratio in transverse direction."""

    # -- Punching shear (pier) ---------------------------------------------
    VEd_punch_pier: float
    """Applied ULS punching load at pier (kN)."""

    VRdc_punch_pier: float
    """Punching shear resistance at pier (kN)."""

    u_pier: float
    """Control perimeter around pier (mm)."""

    punch_util_pier: float
    """Punching utilisation ratio at pier."""

    # -- Punching shear (pile) ---------------------------------------------
    VEd_punch_pile: float
    """Applied ULS punching load at pile (kN)."""

    VRdc_punch_pile: float
    """Punching shear resistance at pile (kN)."""

    u_pile: float
    """Control perimeter around pile (mm)."""

    punch_util_pile: float
    """Punching utilisation ratio at pile."""

    # -- Crack width -------------------------------------------------------
    crack_width: float
    """Calculated crack width under SLS (mm)."""

    crack_width_limit: float
    """Permissible crack width for exposure (mm)."""

    crack_ok: bool
    """True if crack width <= limit."""

    # -- Minimum reinforcement ---------------------------------------------
    As_min: float
    """Minimum reinforcement area (mm2)."""

    min_reinf_ok: bool
    """True if provided steel >= As_min."""

    # -- Overall status ----------------------------------------------------
    status: str
    """'OK' if all checks pass, otherwise 'NOT OK'."""


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
        0.1512 * fck * b * xu^2 - 0.36 * fck * b * d * xu + Mu = 0

    Returns the smaller root (under-reinforced solution).
    """
    a_coeff = 0.36 * fck * b * 0.42   # 0.1512 * fck * b
    b_coeff = 0.36 * fck * b * d
    c_coeff = Mu_Nmm

    discriminant = b_coeff ** 2 - 4.0 * a_coeff * c_coeff
    if discriminant < 0.0:
        # Section inadequate
        return _XU_D_LIMIT_FE500 * d

    xu = (b_coeff - math.sqrt(discriminant)) / (2.0 * a_coeff)
    return max(xu, 0.0)


def _moment_capacity(
    fck: float,
    b: float,
    xu: float,
    d: float,
) -> float:
    """Return moment capacity Mu in N.mm for a given xu.

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
    """
    lever = d - 0.42 * xu
    if lever <= 0.0:
        lever = 0.1 * d  # fallback
    return Mu_Nmm / (0.87 * fyk * lever)


def _provide_steel(
    Ast_reqd: float,
    bar_dia: float = _DEFAULT_BAR_DIA,
) -> tuple[float, int]:
    """Round up the required steel area to whole bars.

    Returns (area_provided_mm2, number_of_bars).
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

    Returns VRd,c in N (not kN).
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


# ---------------------------------------------------------------------------
# Internal helpers -- punching shear (IRC 112 Cl. 10.4)
# ---------------------------------------------------------------------------

def _punching_perimeter_circular(
    D_pier_mm: float,
    d_eff_mm: float,
) -> float:
    """Control perimeter for circular pier at 2d from face.

    u = pi * (D_pier + 4*d)
    """
    return math.pi * (D_pier_mm + 4.0 * d_eff_mm)


def _punching_perimeter_rectangular(
    b_pier_mm: float,
    h_pier_mm: float,
    d_eff_mm: float,
) -> float:
    """Control perimeter for rectangular pier at 2d from face.

    u = 2 * (b + h) + 2 * pi * 2 * d
    """
    return 2.0 * (b_pier_mm + h_pier_mm) + 2.0 * math.pi * 2.0 * d_eff_mm


def _punching_perimeter_pile(
    D_pile_mm: float,
    d_eff_mm: float,
) -> float:
    """Control perimeter for pile at 2d from face.

    For interior piles: u = pi * (D_pile + 4*d)
    """
    return math.pi * (D_pile_mm + 4.0 * d_eff_mm)


def _vRdc_punching(
    fck: float,
    rho_l: float,
    d_eff_mm: float,
) -> float:
    """Punching shear stress capacity v_Rd,c in MPa.

    v_Rd,c = max(
        0.12 * k * (80 * rho_l * fck)^(1/3),
        0.031 * k^(1.5) * sqrt(fck)
    )
    """
    k = min(1.0 + math.sqrt(200.0 / d_eff_mm), 2.0)

    v1 = 0.12 * k * (80.0 * rho_l * fck) ** (1.0 / 3.0)
    v2 = 0.031 * k ** 1.5 * math.sqrt(fck)

    return max(v1, v2)


# ---------------------------------------------------------------------------
# Internal helpers -- crack width (IRC 112 Cl. 12.3.4)
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

    Returns crack width in mm.
    """
    if alpha_e_override is not None:
        alpha_e = alpha_e_override
    else:
        Ecm_MPa = Ecm_GPa * 1000.0
        alpha_e = Es_MPa / Ecm_MPa

    # Cracked-section neutral axis
    a_q = 0.5 * b
    b_q = alpha_e * Ast
    c_q = -alpha_e * Ast * d

    disc = b_q ** 2 - 4.0 * a_q * c_q
    if disc < 0.0:
        x_cr = 0.4 * d
    else:
        x_cr = (-b_q + math.sqrt(disc)) / (2.0 * a_q)

    z_cr = d - x_cr / 3.0

    # Steel stress
    if Ast > 0.0 and z_cr > 0.0:
        sigma_s = M_sls_Nmm / (Ast * z_cr)
    else:
        sigma_s = 0.0

    # Effective tension area
    h_c_eff = min(2.5 * (h - d), h / 2.0, (h - x_cr) / 3.0)
    h_c_eff = max(h_c_eff, 0.0)
    Ac_eff = b * h_c_eff
    rho_p_eff = Ast / Ac_eff if Ac_eff > 0.0 else 0.01

    # Strain difference
    kt = 0.5
    if rho_p_eff > 0.0 and Es_MPa > 0.0:
        eps_diff = (
            sigma_s / Es_MPa
            - kt * fctm / rho_p_eff * (1.0 + alpha_e * rho_p_eff) / Es_MPa
        )
        eps_min = 0.6 * sigma_s / Es_MPa
        eps_diff = max(eps_diff, eps_min)
    else:
        eps_diff = 0.0

    # Maximum crack spacing
    k1 = 0.8
    k2 = 0.5
    if rho_p_eff > 0.0:
        sr_max = 3.4 * cover + 0.425 * k1 * k2 * bar_dia / rho_p_eff
    else:
        sr_max = 3.4 * cover + 0.425 * k1 * k2 * bar_dia / 0.01

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

    As_min = max(0.26 * fctm / fyk * bt * d, 0.0013 * bt * d)
    """
    return max(
        0.26 * fctm / fyk * bt * d,
        0.0013 * bt * d,
    )


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def design_pilecap(
    config: dict[str, Any],
    geometry: GeometryResults,
    combinations: CombinationResults,
    pile_capacity_result: PileCapacityResult,
    creep_result: Any = None,
) -> PilecapDesignResult:
    """Design the pile cap for flexure, shear, punching, and crack width.

    This is the primary public entry point. It reads pile cap geometry
    from *config* and *geometry*, extracts design forces from the
    governing load combinations and pile capacity results, and performs
    all IRC 112:2020 checks.

    Parameters
    ----------
    config : dict
        Validated configuration dictionary (from input_parser).
    geometry : GeometryResults
        Pre-computed geometry (from geometry.calculate_geometry).
    combinations : CombinationResults
        Load combinations (from load_combinations.calculate_combinations).
    pile_capacity_result : PileCapacityResult
        Pile capacity analysis results (from pile_capacity.calculate_pile_capacity).
    creep_result : CreepResult or None
        Creep calculation results.  If provided, the long-term modular
        ratio is used for crack width calculations.

    Returns
    -------
    PilecapDesignResult
        Comprehensive design results including reinforcement areas,
        capacity ratios, crack width, and overall pass/fail status.
    """
    irc = load_irc_tables()

    # ==================================================================
    # 1. Material properties
    # ==================================================================
    mat_cfg = config["materials"]
    fck = float(mat_cfg["pile_pilecap"]["fck"])
    aggregate = str(mat_cfg.get("pile_pilecap", {}).get("aggregate", "Quartzite"))
    exposure = str(mat_cfg["exposure"])
    fyk = float(mat_cfg["steel"]["fyk"])
    Es_input = float(mat_cfg["steel"]["Es"])

    # Normalise Es to GPa
    Es_GPa = Es_input if Es_input < 1000.0 else Es_input / 1000.0

    concrete = get_concrete_properties(fck, aggregate=aggregate, exposure=exposure)
    steel = get_steel_properties(fyk=fyk, Es=Es_GPa)

    # ==================================================================
    # 2. Pile cap geometry (all in mm)
    # ==================================================================
    pilecap_length_long_mm = geometry.pilecap_length_long * 1000.0
    pilecap_width_trans_mm = geometry.pilecap_width_trans * 1000.0
    pilecap_thickness_mm = geometry.pilecap_thickness * 1000.0

    # Cover
    cover_mm = float(config["foundation"].get("cover", 75.0))  # mm

    bar_dia = _DEFAULT_BAR_DIA

    # Effective depth
    d_eff_mm = pilecap_thickness_mm - cover_mm - bar_dia / 2.0
    d_eff_mm = max(d_eff_mm, 100.0)  # safety

    # ==================================================================
    # 3. Find governing ULS combination
    # ==================================================================
    # Use the combo that gives maximum pile reaction
    gov_combo_name = pile_capacity_result.max_compression_combo
    gov_pile_loads: list[PileLoad] = []

    for pile_combo in pile_capacity_result.all_pile_results:
        if pile_combo.combo_name == gov_combo_name:
            gov_pile_loads = pile_combo.pile_loads
            break

    # Get total applied forces at pilecap top
    gov_uls_combo: CombinationResult | None = None
    for combo in combinations.all_combinations:
        if combo.name == gov_combo_name:
            gov_uls_combo = combo
            break

    if gov_uls_combo is None:
        # Fallback to governing ULS basic
        gov_uls_combo = combinations.governing_uls_basic
        # Re-find pile loads
        for pile_combo in pile_capacity_result.all_pile_results:
            if pile_combo.combo_name == gov_uls_combo.name:
                gov_pile_loads = pile_combo.pile_loads
                break

    total_P_uls = gov_uls_combo.forces_pilecap_top.P  # kN

    # Get SLS forces for crack width (quasi-permanent per IRC 112 Cl. 12.3.4)
    gov_sls_combo = combinations.governing_sls_quasi
    total_P_sls = gov_sls_combo.forces_pilecap_top.P

    # ==================================================================
    # 4. Pier dimensions (from geometry)
    # ==================================================================
    pier_section: PierSection = geometry.pier_section
    pier_type = pier_section.pier_type

    if pier_type == "Circular":
        pier_dim_long_mm = pier_section.width_long * 1000.0
        pier_dim_trans_mm = pier_section.width_trans * 1000.0
    else:
        pier_dim_long_mm = pier_section.width_long * 1000.0
        pier_dim_trans_mm = pier_section.width_trans * 1000.0

    # ==================================================================
    # 5. Flexure design (both directions)
    # ==================================================================
    # Critical section at face of pier
    # Longitudinal direction: piles beyond pier face in x-direction
    # Transverse direction: piles beyond pier face in y-direction

    pier_face_offset_long_m = pier_section.width_long / 2.0  # m
    pier_face_offset_trans_m = pier_section.width_trans / 2.0  # m

    Mu_long = 0.0  # kN.m
    Mu_trans = 0.0  # kN.m

    for pile_load in gov_pile_loads:
        x_pile = pile_load.pile_coord.x  # m
        y_pile = pile_load.pile_coord.y  # m
        P_pile = pile_load.P  # kN

        # Longitudinal moment: piles beyond pier face in x-direction
        if abs(x_pile) > pier_face_offset_long_m:
            lever_x = abs(x_pile) - pier_face_offset_long_m
            Mu_long += P_pile * lever_x

        # Transverse moment: piles beyond pier face in y-direction
        if abs(y_pile) > pier_face_offset_trans_m:
            lever_y = abs(y_pile) - pier_face_offset_trans_m
            Mu_trans += P_pile * lever_y

    # Convert to N.mm
    Mu_long_Nmm = Mu_long * 1.0e6
    Mu_trans_Nmm = Mu_trans * 1.0e6

    # --- Longitudinal direction ---
    b_long = pilecap_width_trans_mm
    d_long = d_eff_mm

    xu_long = _solve_xu(Mu_long_Nmm, fck, b_long, d_long)
    xu_d_long = xu_long / d_long
    if xu_d_long > _XU_D_LIMIT_FE500:
        xu_long = _XU_D_LIMIT_FE500 * d_long

    Ast_reqd_long = _ast_required(Mu_long_Nmm, fyk, xu_long, d_long)
    # Enforce minimum reinforcement per IRC 112 Cl. 16.5.1.1
    As_min_long_early = _min_reinforcement(
        get_concrete_properties(fck).fctm, fyk, b_long, d_long,
    )
    Ast_reqd_long = max(Ast_reqd_long, As_min_long_early)
    Ast_prov_long, _ = _provide_steel(Ast_reqd_long, bar_dia)

    # Capacity with provided steel
    xu_from_Ast_long = 0.87 * fyk * Ast_prov_long / (0.36 * fck * b_long)
    xu_from_Ast_long = min(xu_from_Ast_long, _XU_D_LIMIT_FE500 * d_long)
    Mu_cap_long_Nmm = _moment_capacity(fck, b_long, xu_from_Ast_long, d_long)

    flexure_util_long = Mu_long_Nmm / Mu_cap_long_Nmm if Mu_cap_long_Nmm > 0 else 999.0

    # --- Transverse direction ---
    b_trans = pilecap_length_long_mm
    d_trans = d_eff_mm

    xu_trans = _solve_xu(Mu_trans_Nmm, fck, b_trans, d_trans)
    xu_d_trans = xu_trans / d_trans
    if xu_d_trans > _XU_D_LIMIT_FE500:
        xu_trans = _XU_D_LIMIT_FE500 * d_trans

    Ast_reqd_trans = _ast_required(Mu_trans_Nmm, fyk, xu_trans, d_trans)
    As_min_trans_early = _min_reinforcement(
        get_concrete_properties(fck).fctm, fyk, b_trans, d_trans,
    )
    Ast_reqd_trans = max(Ast_reqd_trans, As_min_trans_early)
    Ast_prov_trans, _ = _provide_steel(Ast_reqd_trans, bar_dia)

    # Capacity with provided steel
    xu_from_Ast_trans = 0.87 * fyk * Ast_prov_trans / (0.36 * fck * b_trans)
    xu_from_Ast_trans = min(xu_from_Ast_trans, _XU_D_LIMIT_FE500 * d_trans)
    Mu_cap_trans_Nmm = _moment_capacity(fck, b_trans, xu_from_Ast_trans, d_trans)

    flexure_util_trans = Mu_trans_Nmm / Mu_cap_trans_Nmm if Mu_cap_trans_Nmm > 0 else 999.0

    # ==================================================================
    # 6. One-way shear (both directions)
    # ==================================================================
    # Standard critical section at d from pier face (piles beyond d)
    d_m = d_eff_mm / 1000.0  # convert to m

    VEd_long = 0.0  # kN
    VEd_trans = 0.0  # kN

    for pile_load in gov_pile_loads:
        x_pile = pile_load.pile_coord.x
        y_pile = pile_load.pile_coord.y
        P_pile = pile_load.P

        # Longitudinal: count piles on the MORE loaded side only
        if abs(x_pile) > (pier_face_offset_long_m + d_m):
            VEd_long += P_pile

        # Transverse: count piles on the MORE loaded side only
        if abs(y_pile) > (pier_face_offset_trans_m + d_m):
            VEd_trans += P_pile

    # Concrete shear capacity
    VRdc_long_N = _VRdc(fck, b_long, d_long, Ast_prov_long)
    VRdc_trans_N = _VRdc(fck, b_trans, d_trans, Ast_prov_trans)

    VRdc_long = VRdc_long_N / 1000.0  # N -> kN
    VRdc_trans = VRdc_trans_N / 1000.0

    shear_util_long = VEd_long / VRdc_long if VRdc_long > 0 else 999.0
    shear_util_trans = VEd_trans / VRdc_trans if VRdc_trans > 0 else 999.0

    # ==================================================================
    # 7-8. Punching shear per IRC 112 Cl. 10.4
    # ==================================================================
    # Four checks following the Excel methodology:
    #   A) Pier punching at pier face → vs VRd,max
    #   B) Pier punching at 2d control perimeter → vs VRd,c (enhanced 2d/av)
    #   C) Pile punching at pile face → vs VRd,max
    #   D) Pile punching at perimeter "c" at av → vs VRd,c (enhanced 2d/av)
    # Governing = max utilization across all checks.

    # --- Common punching parameters ---
    k_punch = min(1.0 + math.sqrt(200.0 / d_eff_mm), 2.0)
    rho_lx_p = Ast_prov_long / (b_long * d_long)
    rho_ly_p = Ast_prov_trans / (b_trans * d_trans)
    rho_l_p = min(math.sqrt(rho_lx_p * rho_ly_p), 0.02)

    # Base VRd,c in MPa (without enhancement)
    _v1 = 0.12 * k_punch * (80.0 * rho_l_p * fck) ** (1.0 / 3.0)
    _v2 = 0.031 * k_punch ** 1.5 * math.sqrt(fck)
    vRdc_base_MPa = max(_v1, _v2)

    # VRd,max in MPa (strut crushing limit)
    _nu = 0.6 * (1.0 - fck / 310.0)
    _fcd = 0.67 * fck / 1.5
    vRd_max_MPa = 0.5 * _nu * _fcd

    D_pile_mm = geometry.pile_props.diameter * 1000.0
    pile_R_mm = D_pile_mm / 2.0
    density_kNm3 = float(config["materials"].get("concrete_density", 25.0))

    # --- av: clear distance from pile face to pier face (diagonal) ---
    max_pile = pile_capacity_result.max_compression_pile
    px_m = max_pile.pile_coord.x
    py_m = max_pile.pile_coord.y
    dist_pile_pier_m = math.sqrt(px_m ** 2 + py_m ** 2)

    if pier_type == "Circular":
        pier_R_m = pier_dim_long_mm / 2000.0
        av_m = dist_pile_pier_m - pier_R_m - pile_R_mm / 1000.0
    else:
        pier_half_L = pier_dim_long_mm / 2000.0
        pier_half_T = pier_dim_trans_mm / 2000.0
        angle = math.atan2(abs(py_m), max(abs(px_m), 1e-9))
        d_face = min(
            pier_half_L / max(abs(math.cos(angle)), 1e-9),
            pier_half_T / max(abs(math.sin(angle)), 1e-9),
        )
        av_m = dist_pile_pier_m - d_face - pile_R_mm / 1000.0

    av_m = max(av_m, 0.001)
    av_mm = av_m * 1000.0

    # Enhancement factor 2d/av (IRC 112 Cl. 10.4.5 / Cl. 10.3.3.3)
    enhancement_2d_av = 2.0 * d_eff_mm / av_mm
    vRdc_enhanced_MPa = vRdc_base_MPa * enhancement_2d_av

    # --- Check A: Pier punching at face ---
    if pier_type == "Circular":
        u0_pier_mm = math.pi * pier_dim_long_mm
    else:
        u0_pier_mm = 2.0 * (pier_dim_long_mm + pier_dim_trans_mm)

    VEd_pier_total = sum(pl.P for pl in gov_pile_loads)  # kN
    M_gov = math.sqrt(
        gov_uls_combo.forces_pilecap_top.ML ** 2
        + gov_uls_combo.forces_pilecap_top.MT ** 2
    )
    e_pier_m = M_gov / max(VEd_pier_total, 1.0)
    D_eff_pier_m = (pier_dim_long_mm + 4.0 * d_eff_mm) / 1000.0
    beta_pier_face = min(1.0 + 0.6 * math.pi * e_pier_m / max(pier_dim_long_mm / 1000.0, 0.1), 2.5)

    ved_A = beta_pier_face * VEd_pier_total * 1e3 / (u0_pier_mm * d_eff_mm)
    util_A = ved_A / vRd_max_MPa if vRd_max_MPa > 0 else 999.0

    # --- Check B: Pier punching at 2d perimeter ---
    if pier_type == "Circular":
        u1_pier_mm = math.pi * (pier_dim_long_mm + 4.0 * d_eff_mm)
    else:
        u1_pier_mm = (2.0 * (pier_dim_long_mm + pier_dim_trans_mm)
                      + 2.0 * math.pi * 2.0 * d_eff_mm)

    beta_pier_2d = min(1.0 + 0.6 * math.pi * e_pier_m / D_eff_pier_m, 2.5)
    ved_B = beta_pier_2d * VEd_pier_total * 1e3 / (u1_pier_mm * d_eff_mm)
    util_B = ved_B / vRdc_enhanced_MPa if vRdc_enhanced_MPa > 0 else 999.0

    # --- Check C: Pile punching at pile face ---
    VEd_pile_raw = max_pile.P  # kN
    pile_plan_area_m2 = math.pi * (geometry.pile_props.diameter / 2.0) ** 2
    sw_deduct_face = pile_plan_area_m2 * (pilecap_thickness_mm / 1000.0) * density_kNm3
    VEd_pile_face = VEd_pile_raw - sw_deduct_face

    u0_pile_mm = math.pi * D_pile_mm
    beta_pile_face = 1.15  # typical for corner piles

    ved_C = beta_pile_face * VEd_pile_face * 1e3 / (u0_pile_mm * d_eff_mm)
    util_C = ved_C / vRd_max_MPa if vRd_max_MPa > 0 else 999.0

    # --- Check D: Pile punching at perimeter "c" at av from pile face ---
    # Truncated control perimeter at av from pile face.
    # For corner piles (truncated by 2 pilecap edges): quarter-circle + 2 straights
    # For edge piles (truncated by 1 edge): half-circle + 1 straight
    # For interior piles: full circle
    R_control = pile_R_mm + av_mm  # radius of control circle (mm)
    c_edge_x = pilecap_length_long_mm / 2.0 - abs(px_m) * 1000.0  # mm
    c_edge_y = pilecap_width_trans_mm / 2.0 - abs(py_m) * 1000.0  # mm

    n_edges_trunc = int(c_edge_x < R_control) + int(c_edge_y < R_control)

    if n_edges_trunc >= 2:
        c_edge_mm = min(c_edge_x, c_edge_y)
        u_pile_c_mm = (2.0 * c_edge_mm
                       + math.pi / 4.0 * (D_pile_mm + 2.0 * av_mm))
    elif n_edges_trunc == 1:
        c_edge_mm = min(c_edge_x, c_edge_y)
        u_pile_c_mm = (c_edge_mm
                       + math.pi / 2.0 * (D_pile_mm + 2.0 * av_mm))
    else:
        u_pile_c_mm = math.pi * (D_pile_mm + 2.0 * av_mm)

    # Deduct pilecap self-weight over control perimeter plan area
    plan_area_c_m2 = u_pile_c_mm / 1000.0  # Excel convention
    sw_deduct_c = plan_area_c_m2 * (pilecap_thickness_mm / 1000.0) * density_kNm3
    VEd_pile_c = VEd_pile_raw - sw_deduct_c

    # Beta for pile eccentricity per IRC 112 Eq. 10.29:
    # β = 1 + k * e * u1/W1.  For a pile in a group, the eccentricity is
    # primarily from the moment variation across the pile group.
    # Compute from the actual forces: e = M_resultant / V_total at pilecap
    if VEd_pier_total > 0:
        e_pile_m = M_gov / VEd_pier_total
        D_control_m = (D_pile_mm + 2.0 * av_mm) / 1000.0
        beta_pile_c = 1.0 + 0.6 * e_pile_m * dist_pile_pier_m / (D_control_m ** 2)
        beta_pile_c = min(max(beta_pile_c, 1.0), 1.5)
    else:
        beta_pile_c = 1.10
    ved_D = beta_pile_c * VEd_pile_c * 1e3 / (u_pile_c_mm * d_eff_mm)

    # --- Ensure reinforcement is sufficient for punching ---
    # Back-calculate the minimum rho_l needed so vRdc_enhanced >= ved_D.
    # vRdc_base * (2d/av) >= ved_D  →  vRdc_base >= ved_D / enhancement_2d_av
    vRdc_target = ved_D / enhancement_2d_av if enhancement_2d_av > 0 else 999.0
    # From 0.12*k*(80*rho*fck)^(1/3) >= vRdc_target
    rhs_cubed = (vRdc_target / (0.12 * k_punch)) ** 3  # = 80*rho*fck
    rho_punch_reqd = rhs_cubed / (80.0 * fck) if fck > 0 else 0.0
    rho_punch_reqd = max(rho_punch_reqd, 0.0)

    # If current rho is insufficient, increase reinforcement
    Ast_punch_reqd_long = rho_punch_reqd * b_long * d_long
    Ast_punch_reqd_trans = rho_punch_reqd * b_trans * d_trans

    if Ast_prov_long < Ast_punch_reqd_long:
        Ast_prov_long, _ = _provide_steel(Ast_punch_reqd_long, bar_dia)
        # Recompute xu for capacity
        xu_from_Ast_long = 0.87 * fyk * Ast_prov_long / (0.36 * fck * b_long)
        xu_from_Ast_long = min(xu_from_Ast_long, _XU_D_LIMIT_FE500 * d_long)
        Mu_cap_long_Nmm = _moment_capacity(fck, b_long, xu_from_Ast_long, d_long)
        flexure_util_long = Mu_long_Nmm / Mu_cap_long_Nmm if Mu_cap_long_Nmm > 0 else 999.0

    if Ast_prov_trans < Ast_punch_reqd_trans:
        Ast_prov_trans, _ = _provide_steel(Ast_punch_reqd_trans, bar_dia)
        xu_from_Ast_trans = 0.87 * fyk * Ast_prov_trans / (0.36 * fck * b_trans)
        xu_from_Ast_trans = min(xu_from_Ast_trans, _XU_D_LIMIT_FE500 * d_trans)
        Mu_cap_trans_Nmm = _moment_capacity(fck, b_trans, xu_from_Ast_trans, d_trans)
        flexure_util_trans = Mu_trans_Nmm / Mu_cap_trans_Nmm if Mu_cap_trans_Nmm > 0 else 999.0

    # Recompute punching with updated reinforcement
    rho_lx_p = Ast_prov_long / (b_long * d_long)
    rho_ly_p = Ast_prov_trans / (b_trans * d_trans)
    rho_l_p = min(math.sqrt(rho_lx_p * rho_ly_p), 0.02)

    _v1 = 0.12 * k_punch * (80.0 * rho_l_p * fck) ** (1.0 / 3.0)
    _v2 = 0.031 * k_punch ** 1.5 * math.sqrt(fck)
    vRdc_base_MPa = max(_v1, _v2)
    vRdc_enhanced_MPa = vRdc_base_MPa * enhancement_2d_av

    # Recompute all utilizations with updated VRd,c
    util_A = ved_A / vRd_max_MPa if vRd_max_MPa > 0 else 999.0
    util_B = ved_B / vRdc_enhanced_MPa if vRdc_enhanced_MPa > 0 else 999.0
    util_C = ved_C / vRd_max_MPa if vRd_max_MPa > 0 else 999.0
    util_D = ved_D / vRdc_enhanced_MPa if vRdc_enhanced_MPa > 0 else 999.0

    # --- Governing utilizations ---
    punch_util_pier = max(util_A, util_B)
    punch_util_pile = max(util_C, util_D)

    # Store representative values for reporting
    VEd_punch_pier = VEd_pier_total
    VRdc_punch_pier = vRdc_enhanced_MPa * u1_pier_mm * d_eff_mm / 1e3  # kN
    u_pier_mm = u1_pier_mm

    VEd_punch_pile = VEd_pile_c
    VRdc_punch_pile = vRdc_enhanced_MPa * u_pile_c_mm * d_eff_mm / 1e3  # kN
    u_pile_mm = u_pile_c_mm

    # ==================================================================
    # 9. Crack width (SLS)
    # ==================================================================
    # Use longitudinal direction at pier face for crack width check
    # SLS moment: approximate from SLS total load
    ratio_sls_to_uls = total_P_sls / total_P_uls if total_P_uls > 0 else 1.0
    M_sls_long_Nmm = Mu_long_Nmm * ratio_sls_to_uls

    Es_MPa = Es_GPa * 1000.0

    # Long-term modular ratio: alpha_e = Es / Ec_eff
    # where Ec_eff = Ecm / (1 + phi_creep)
    alpha_e_lt = None
    if creep_result is not None and hasattr(creep_result, 'phi'):
        Ec_eff_MPa = concrete.Ecm * 1000.0 / (1.0 + creep_result.phi)
        alpha_e_lt = Es_MPa / Ec_eff_MPa

    wk = _crack_width(
        M_sls_Nmm=M_sls_long_Nmm,
        b=b_long,
        h=pilecap_thickness_mm,
        d=d_eff_mm,
        Ast=Ast_prov_long,
        cover=cover_mm,
        bar_dia=bar_dia,
        fctm=concrete.fctm,
        Es_MPa=Es_MPa,
        Ecm_GPa=concrete.Ecm,
        alpha_e_override=alpha_e_lt,
    )

    # Crack width limit
    crack_limits: dict[str, float] = irc.get("crack_width_limits", {})
    wk_limit = float(crack_limits.get(exposure, 0.3))
    crack_ok = wk <= wk_limit

    # ==================================================================
    # 10. Minimum reinforcement
    # ==================================================================
    As_min_long = _min_reinforcement(concrete.fctm, fyk, b_long, d_long)
    As_min_trans = _min_reinforcement(concrete.fctm, fyk, b_trans, d_trans)
    As_min = max(As_min_long, As_min_trans)

    min_reinf_ok = (
        Ast_prov_long >= As_min_long and Ast_prov_trans >= As_min_trans
    )

    # ==================================================================
    # 11. Overall status
    # ==================================================================
    flexure_ok = flexure_util_long <= 1.0 and flexure_util_trans <= 1.0
    shear_ok = shear_util_long <= 1.0 and shear_util_trans <= 1.0
    punch_ok = punch_util_pier <= 1.0 and punch_util_pile <= 1.0
    all_ok = flexure_ok and shear_ok and punch_ok and crack_ok and min_reinf_ok
    status = "OK" if all_ok else "NOT OK"

    # ==================================================================
    # 12. Assemble result
    # ==================================================================
    return PilecapDesignResult(
        # Geometry
        length_long=pilecap_length_long_mm,
        width_trans=pilecap_width_trans_mm,
        thickness=pilecap_thickness_mm,
        cover=cover_mm,
        d_eff=d_eff_mm,
        fck=fck,
        # Flexure longitudinal
        Mu_long=Mu_long,
        Mu_cap_long=Mu_cap_long_Nmm / 1.0e6,
        Ast_reqd_long=Ast_reqd_long,
        Ast_prov_long=Ast_prov_long,
        flexure_util_long=flexure_util_long,
        # Flexure transverse
        Mu_trans=Mu_trans,
        Mu_cap_trans=Mu_cap_trans_Nmm / 1.0e6,
        Ast_reqd_trans=Ast_reqd_trans,
        Ast_prov_trans=Ast_prov_trans,
        flexure_util_trans=flexure_util_trans,
        # One-way shear
        VEd_long=VEd_long,
        VRdc_long=VRdc_long,
        shear_util_long=shear_util_long,
        VEd_trans=VEd_trans,
        VRdc_trans=VRdc_trans,
        shear_util_trans=shear_util_trans,
        # Punching pier
        VEd_punch_pier=VEd_punch_pier,
        VRdc_punch_pier=VRdc_punch_pier,
        u_pier=u_pier_mm,
        punch_util_pier=punch_util_pier,
        # Punching pile
        VEd_punch_pile=VEd_punch_pile,
        VRdc_punch_pile=VRdc_punch_pile,
        u_pile=u_pile_mm,
        punch_util_pile=punch_util_pile,
        # Crack width
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

def summarise_pilecap_design(result: PilecapDesignResult) -> str:
    """Return a human-readable summary of the pile cap design.

    Useful for console output and report generation.

    Parameters
    ----------
    result : PilecapDesignResult
        Output from design_pilecap.

    Returns
    -------
    str
        Multi-line formatted summary string.
    """
    lines: list[str] = []
    lines.append("=" * 80)
    lines.append("PILE CAP DESIGN SUMMARY (IRC 112:2020)")
    lines.append("=" * 80)

    lines.append("")
    lines.append("  GEOMETRY")
    lines.append(f"    Length (longitudinal)    : {result.length_long:>10.0f} mm")
    lines.append(f"    Width (transverse)       : {result.width_trans:>10.0f} mm")
    lines.append(f"    Thickness                : {result.thickness:>10.0f} mm")
    lines.append(f"    Clear cover              : {result.cover:>10.0f} mm")
    lines.append(f"    Effective depth          : {result.d_eff:>10.1f} mm")
    lines.append(f"    fck                      : {result.fck:>10.1f} MPa")

    lines.append("")
    lines.append("  ULS FLEXURE (IRC 112 Cl. 8.2)")
    lines.append("  " + "-" * 76)
    lines.append(f"    {'':35s}  {'Longitudinal':>18s}  {'Transverse':>18s}")
    lines.append(f"    {'Mu applied (kN.m)':35s}  {result.Mu_long:>18.1f}  {result.Mu_trans:>18.1f}")
    lines.append(f"    {'Mu capacity (kN.m)':35s}  {result.Mu_cap_long:>18.1f}  {result.Mu_cap_trans:>18.1f}")
    lines.append(f"    {'Ast required (mm2)':35s}  {result.Ast_reqd_long:>18.0f}  {result.Ast_reqd_trans:>18.0f}")
    lines.append(f"    {'Ast provided (mm2)':35s}  {result.Ast_prov_long:>18.0f}  {result.Ast_prov_trans:>18.0f}")
    lines.append(f"    {'Utilisation':35s}  {result.flexure_util_long:>18.3f}  {result.flexure_util_trans:>18.3f}")

    flexure_status = "OK" if result.flexure_util_long <= 1.0 and result.flexure_util_trans <= 1.0 else "NOT OK"
    lines.append(f"    {'Status':35s}  {flexure_status:>39s}")

    lines.append("")
    lines.append("  ULS ONE-WAY SHEAR (IRC 112 Cl. 10.3)")
    lines.append("  " + "-" * 76)
    lines.append(f"    {'':35s}  {'Longitudinal':>18s}  {'Transverse':>18s}")
    lines.append(f"    {'VEd (kN)':35s}  {result.VEd_long:>18.1f}  {result.VEd_trans:>18.1f}")
    lines.append(f"    {'VRd,c (kN)':35s}  {result.VRdc_long:>18.1f}  {result.VRdc_trans:>18.1f}")
    lines.append(f"    {'Utilisation':35s}  {result.shear_util_long:>18.3f}  {result.shear_util_trans:>18.3f}")

    shear_status = "OK" if result.shear_util_long <= 1.0 and result.shear_util_trans <= 1.0 else "NOT OK"
    lines.append(f"    {'Status':35s}  {shear_status:>39s}")

    lines.append("")
    lines.append("  ULS PUNCHING SHEAR (IRC 112 Cl. 10.4)")
    lines.append("  " + "-" * 76)
    lines.append(f"    {'':35s}  {'Around Pier':>18s}  {'Around Pile':>18s}")
    lines.append(f"    {'VEd (kN)':35s}  {result.VEd_punch_pier:>18.1f}  {result.VEd_punch_pile:>18.1f}")
    lines.append(f"    {'VRd,c (kN)':35s}  {result.VRdc_punch_pier:>18.1f}  {result.VRdc_punch_pile:>18.1f}")
    lines.append(f"    {'Control perimeter (mm)':35s}  {result.u_pier:>18.1f}  {result.u_pile:>18.1f}")
    lines.append(f"    {'Utilisation':35s}  {result.punch_util_pier:>18.3f}  {result.punch_util_pile:>18.3f}")

    punch_status = "OK" if result.punch_util_pier <= 1.0 and result.punch_util_pile <= 1.0 else "NOT OK"
    lines.append(f"    {'Status':35s}  {punch_status:>39s}")

    lines.append("")
    lines.append("  SLS CRACK WIDTH (IRC 112 Cl. 12.3.4)")
    lines.append("  " + "-" * 76)
    lines.append(f"    {'Calculated wk (mm)':35s}  {result.crack_width:>18.4f}")
    lines.append(f"    {'Limit wk (mm)':35s}  {result.crack_width_limit:>18.4f}")
    lines.append(f"    {'Status':35s}  {'OK' if result.crack_ok else 'NOT OK':>18s}")

    lines.append("")
    lines.append("  MINIMUM REINFORCEMENT (IRC 112 Cl. 16.5.1.1)")
    lines.append("  " + "-" * 76)
    lines.append(f"    {'As,min (mm2)':35s}  {result.As_min:>18.0f}")
    lines.append(f"    {'As,prov long (mm2)':35s}  {result.Ast_prov_long:>18.0f}")
    lines.append(f"    {'As,prov trans (mm2)':35s}  {result.Ast_prov_trans:>18.0f}")
    lines.append(f"    {'Status':35s}  {'OK' if result.min_reinf_ok else 'NOT OK':>18s}")

    lines.append("")
    lines.append("=" * 80)
    lines.append(f"  OVERALL STATUS : {result.status}")
    lines.append("=" * 80)

    return "\n".join(lines)

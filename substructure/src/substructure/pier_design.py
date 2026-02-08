"""Pier column design under biaxial bending per IRC 112:2020.

Designs the pier shaft for the governing ULS and SLS load combinations,
checking slenderness, second-order effects, biaxial interaction, stress
limits, crack widths, and ductile detailing requirements.

Key references
--------------
* IRC 112:2020, Cl. 8.3   -- Second-order effects and slenderness
* IRC 112:2020, Cl. 15.2  -- Design of columns under axial load + bending
* IRC 112:2020, Cl. 17.2  -- Ductile detailing for seismic zones
* IRC 112:2020, Cl. 12    -- Serviceability limit state checks
* IRC 112:2020, Table 14.1 -- Minimum cover for durability

Sign conventions
----------------
* P   : axial compression, positive downward (kN)
* ML  : moment about the longitudinal (traffic) axis (kN.m)
* MT  : moment about the transverse axis (kN.m)
* All section dimensions in mm unless otherwise noted.
* Forces in kN, moments in kN.m, stresses in MPa.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from .loads import ForceVector
from .geometry import GeometryResults, PierSection
from .load_combinations import CombinationResult, CombinationResults
from .materials import (
    get_concrete_properties,
    get_steel_properties,
    ConcreteProperties,
    SteelProperties,
)
from .interaction import (
    generate_interaction_circular,
    generate_interaction_rectangular,
    InteractionDiagram,
    check_utilisation,
    check_biaxial,
)
from .creep import calculate_creep, calculate_notional_size
from .utils import load_irc_tables


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class PierDesignResult:
    """Comprehensive results from pier column design.

    Contains section geometry, reinforcement details, slenderness check,
    ULS capacity and utilisation, SLS stress and crack width checks,
    ductile detailing assessment, and the interaction diagram data for
    downstream plotting.
    """

    # -- Section geometry ---------------------------------------------------
    pier_type: str          # "Circular" or "Rectangular"
    diameter: float         # mm (diameter for circular, width_long for rectangular)
    width_trans: float      # mm (same as diameter for circular)
    height: float           # m  (pier height)
    fck: float              # MPa
    cover: float            # mm

    # -- Reinforcement layout -----------------------------------------------
    n_bars: int
    bar_dia: float          # mm
    Ast_provided: float     # mm2
    rho_l: float            # longitudinal reinforcement ratio (As / Ac)

    # -- Slenderness --------------------------------------------------------
    lambda_: float          # slenderness ratio
    lambda_lim: float       # limiting slenderness per IRC 112 Cl. 8.3.2
    second_order: bool      # True if second-order effects must be included

    # -- ULS capacity (from interaction diagram) ----------------------------
    P_capacity: float       # kN  -- pure compression capacity
    M_capacity_xx: float    # kN.m -- moment capacity about xx at balanced point
    M_capacity_yy: float    # kN.m -- moment capacity about yy at balanced point

    # -- Governing ULS combination ------------------------------------------
    governing_combo_name: str
    P_applied: float        # kN
    ML_applied: float       # kN.m (first-order moment about longitudinal axis)
    MT_applied: float       # kN.m (first-order moment about transverse axis)
    ML_with_2nd_order: float  # kN.m (including second-order if applicable)
    MT_with_2nd_order: float  # kN.m (including second-order if applicable)

    # -- Utilisation --------------------------------------------------------
    util_uniaxial_xx: float  # ML / MuL at given P
    util_uniaxial_yy: float  # MT / MuT at given P
    util_biaxial: float      # Bresler combined utilisation

    # -- SLS checks ---------------------------------------------------------
    sigma_c_rare: float     # MPa  -- concrete stress under SLS rare
    sigma_c_qp: float       # MPa  -- concrete stress under SLS quasi-permanent
    sigma_s_rare: float     # MPa  -- steel stress under SLS rare
    crack_width: float      # mm   -- estimated crack width under SLS QP
    sls_ok: bool            # True if all SLS checks pass

    # -- Ductile detailing (IRC 112 Cl. 17.2.1) ----------------------------
    Lp: float               # mm  -- plastic hinge length
    spiral_dia: float       # mm  -- confining reinforcement diameter
    spiral_spacing: float   # mm  -- spacing of hoops / spirals
    omega_wd: float         # provided confining reinforcement ratio
    omega_wd_required: float  # required confining reinforcement ratio
    ductile_ok: bool        # True if ductile detailing is adequate

    # -- Overall status -----------------------------------------------------
    status: str             # "OK" or "NOT OK"

    # -- Interaction diagrams for plotting ----------------------------------
    interaction_xx: InteractionDiagram
    interaction_yy: InteractionDiagram


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _get_cover(config: dict[str, Any]) -> float:
    """Determine clear cover in mm from config or IRC 112 Table 14.1.

    Uses the exposure condition from config["materials"]["exposure"] to
    look up the minimum cover, then adds a tolerance of 10 mm per IRC 112
    Cl. 14.4.2 (deviation allowance).

    If a manual ``pier.cover`` key is provided, that value is used directly.
    """
    pier_cfg = config.get("pier", {})
    if "cover" in pier_cfg and pier_cfg["cover"] is not None:
        return float(pier_cfg["cover"])

    exposure = config["materials"].get("exposure", "Moderate")
    irc = load_irc_tables()
    cover_table = irc["cover_requirements"]
    base_cover = float(cover_table.get(exposure, 45))
    # IRC 112 Cl. 14.4.2: add deviation tolerance of 10 mm
    return base_cover + 10.0


def _compute_Ac_mm2(pier_section: PierSection) -> float:
    """Gross concrete area in mm2 from PierSection (which stores area in m2)."""
    return pier_section.area * 1e6


def _compute_Ic_mm4(pier_section: PierSection, axis: str = "xx") -> float:
    """Gross moment of inertia in mm4 from PierSection (stored in m4)."""
    if axis == "xx":
        return pier_section.inertia_xx * 1e12
    return pier_section.inertia_yy * 1e12


def _radius_of_gyration(I_mm4: float, A_mm2: float) -> float:
    """Radius of gyration i = sqrt(I / A) in mm."""
    if A_mm2 <= 0:
        return 1.0
    return math.sqrt(I_mm4 / A_mm2)


def _bar_area(dia_mm: float) -> float:
    """Cross-sectional area of one bar in mm2."""
    return math.pi * dia_mm ** 2 / 4.0


def _select_reinforcement(
    config: dict[str, Any],
    Ac_mm2: float,
    fcd: float,
    fyd: float,
    max_NED_kN: float,
) -> tuple[int, float, float]:
    """Select longitudinal reinforcement for the pier.

    Returns (n_bars, bar_dia_mm, Ast_mm2).

    Strategy:
    1. If config["pier"] specifies n_bars and bar_dia, use them directly.
    2. If config["pier"] specifies reinforcement_pct, size from percentage.
    3. Otherwise start from a default of 20 bars of 32 mm diameter and
       verify against IRC 112 minimum requirements.

    Minimum reinforcement (IRC 112 Cl. 15.2.3.2):
        As_min = max(0.002 * Ac, 0.1 * NED / fyd)
    Maximum reinforcement (IRC 112 Cl. 15.2.3.2):
        As_max = 0.04 * Ac
    """
    pier_cfg = config.get("pier", {})

    # Minimum and maximum per IRC 112
    NEd_N = abs(max_NED_kN) * 1000.0  # convert kN to N
    As_min = max(0.002 * Ac_mm2, 0.1 * NEd_N / fyd)
    As_max = 0.04 * Ac_mm2

    # Option 1: explicit bar specification
    if "n_bars" in pier_cfg and "bar_dia" in pier_cfg:
        n = int(pier_cfg["n_bars"])
        dia = float(pier_cfg["bar_dia"])
        Ast = n * _bar_area(dia)
        # Clip to min/max
        if Ast < As_min:
            # Increase number of bars to meet minimum
            n = math.ceil(As_min / _bar_area(dia))
            Ast = n * _bar_area(dia)
        return n, dia, Ast

    # Option 2: reinforcement percentage
    if "reinforcement_pct" in pier_cfg:
        pct = float(pier_cfg["reinforcement_pct"]) / 100.0
        Ast_target = pct * Ac_mm2
        Ast_target = max(Ast_target, As_min)
        Ast_target = min(Ast_target, As_max)
        # Pick bar size and count
        for dia in [25.0, 28.0, 32.0, 36.0, 40.0]:
            n = math.ceil(Ast_target / _bar_area(dia))
            if n >= 6:  # reasonable minimum number of bars
                Ast = n * _bar_area(dia)
                return n, dia, Ast
        # Fallback
        dia = 32.0
        n = math.ceil(Ast_target / _bar_area(dia))
        n = max(n, 6)
        return n, dia, n * _bar_area(dia)

    # Option 3: default -- start with 20 bars of 32 mm
    n = 20
    dia = 32.0
    Ast = n * _bar_area(dia)
    if Ast < As_min:
        # Increase count
        n = math.ceil(As_min / _bar_area(dia))
        Ast = n * _bar_area(dia)
    if Ast > As_max:
        # Reduce count (should not normally happen with 20x32)
        n = max(6, int(As_max / _bar_area(dia)))
        Ast = n * _bar_area(dia)
    return n, dia, Ast


def _effective_length(config: dict[str, Any], pier_height_m: float) -> float:
    """Compute effective length l0 in mm.

    The end condition is read from config["pier"]["end_condition"]:
      - "fixed-free"   (cantilever): l0 = 2.0 * L
      - "fixed-guided" (sway prevented): l0 = 1.0 * L
      - "fixed-fixed"               : l0 = 0.7 * L

    Default assumption: fixed at base, free at top (cantilever), so l0 = 2.0 * L.
    """
    pier_cfg = config.get("pier", {})
    end_condition = pier_cfg.get("end_condition", "fixed-free")
    L_mm = pier_height_m * 1000.0

    if end_condition == "fixed-fixed":
        return 0.7 * L_mm
    if end_condition == "fixed-guided":
        return 1.0 * L_mm
    # Default: cantilever (fixed base, free top)
    return 2.0 * L_mm


def _slenderness_check(
    l0_mm: float,
    I_mm4: float,
    Ac_mm2: float,
    fck: float,
    fcd: float,
    fyd: float,
    Ast_mm2: float,
    NEd_kN: float,
    phi_eff: float,
    M01_kNm: float,
    M02_kNm: float,
) -> tuple[float, float, bool]:
    """Check slenderness per IRC 112 Cl. 8.3.2.

    Parameters
    ----------
    l0_mm : effective length in mm
    I_mm4 : gross moment of inertia about the axis being checked, mm4
    Ac_mm2 : gross concrete area, mm2
    fck : characteristic concrete strength, MPa
    fcd : design concrete strength, MPa
    fyd : design steel yield strength, MPa
    Ast_mm2 : total longitudinal reinforcement area, mm2
    NEd_kN : design axial force (compression positive), kN
    phi_eff : effective creep ratio
    M01_kNm : smaller end moment (absolute value), kN.m
    M02_kNm : larger end moment (absolute value), kN.m

    Returns
    -------
    (lambda_, lambda_lim, second_order_needed)
    """
    # Radius of gyration
    i = _radius_of_gyration(I_mm4, Ac_mm2)  # mm

    # Slenderness ratio
    lambda_ = l0_mm / i if i > 0 else 999.0

    # Normalised axial force
    NEd_N = abs(NEd_kN) * 1000.0
    n = NEd_N / (Ac_mm2 * fcd) if (Ac_mm2 * fcd) > 0 else 0.01
    n = max(n, 0.01)  # avoid division by zero

    # Mechanical reinforcement ratio
    omega = Ast_mm2 * fyd / (Ac_mm2 * fcd) if (Ac_mm2 * fcd) > 0 else 0.0

    # Factor A: accounts for creep
    A = 1.0 / (1.0 + 0.2 * phi_eff)

    # Factor B: accounts for reinforcement
    B = math.sqrt(1.0 + 2.0 * omega)

    # Factor C: accounts for end-moment ratio
    if abs(M02_kNm) > 1e-6:
        rm = M01_kNm / M02_kNm
    else:
        rm = 1.0
    C = 1.7 - rm

    # Limiting slenderness
    lambda_lim = 20.0 * A * B * C / math.sqrt(n)

    second_order = lambda_ > lambda_lim

    return lambda_, lambda_lim, second_order


def _compute_creep_coefficient(
    config: dict[str, Any],
    concrete: ConcreteProperties,
    pier_section: PierSection,
) -> float:
    """Compute creep coefficient for the pier section.

    Uses the creep module with pier section notional size and parameters
    from config["creep"].
    """
    creep_cfg = config.get("creep", {})
    t0 = float(creep_cfg.get("age_of_loading", 28))
    RH = float(creep_cfg.get("relative_humidity", 70.0))
    t = float(creep_cfg.get("design_life", 25550))

    # Notional size: h0 = 2 * Ac / u in mm
    Ac_mm2 = pier_section.area * 1e6
    perimeter_mm = pier_section.perimeter * 1e3
    h0 = calculate_notional_size(Ac_mm2, perimeter_mm)

    fcm = concrete.fcm
    Ecm = concrete.Ecm  # in GPa
    Es = 200.0           # GPa

    creep_result = calculate_creep(
        fcm=fcm, Ecm=Ecm, Es=Es,
        h0=h0, RH=RH, t0=t0, t=t,
    )
    return creep_result.phi


def _effective_creep_ratio(
    phi: float,
    M0_qp: float,
    M0_Ed: float,
) -> float:
    """Compute effective creep ratio per IRC 112.

    phi_eff = phi * (M0_qp / M0_Ed)

    where M0_qp is the first-order moment under quasi-permanent load
    and M0_Ed is the first-order moment under the design (ULS) load.
    """
    if abs(M0_Ed) < 1e-6:
        return phi  # conservative: assume full creep
    return phi * abs(M0_qp) / abs(M0_Ed)


def _second_order_moment(
    NEd_kN: float,
    l0_mm: float,
    fyd: float,
    Es_GPa: float,
    d_mm: float,
    Kr: float,
    Kphi: float,
) -> float:
    """Compute second-order moment M2 per IRC 112 Cl. 8.3.2.2.

    M2 = NEd * e2
    e2 = (1/r) * l0^2 / c
    1/r = Kr * Kphi * 1/r0
    1/r0 = fyd / (Es * 0.45 * d)
    c = 10 for constant cross-section (pi^2 approximation)

    Parameters
    ----------
    NEd_kN : design axial force, kN (compression positive)
    l0_mm : effective length, mm
    fyd : design yield stress, MPa
    Es_GPa : steel modulus, GPa
    d_mm : effective depth, mm
    Kr : correction factor for axial force
    Kphi : correction factor for creep

    Returns
    -------
    M2 in kN.m
    """
    Es_MPa = Es_GPa * 1000.0
    c = 10.0  # constant cross-section

    # Base curvature
    one_over_r0 = fyd / (Es_MPa * 0.45 * d_mm) if d_mm > 0 else 0.0

    # Design curvature
    one_over_r = Kr * Kphi * one_over_r0

    # Second-order eccentricity
    e2 = one_over_r * l0_mm ** 2 / c  # mm

    # Second-order moment: NEd (kN) * e2 (mm) / 1000 -> kN.m
    M2_kNm = abs(NEd_kN) * e2 / 1000.0

    return M2_kNm


def _compute_Kr(
    NEd_kN: float,
    Ac_mm2: float,
    fcd: float,
    fyd: float,
    Ast_mm2: float,
) -> float:
    """Compute correction factor Kr for nominal curvature method.

    Kr = min((nu - n) / (nu - n_bal), 1.0)
    nu = 1 + omega
    n_bal = 0.4
    n = NEd / (Ac * fcd)
    omega = As * fyd / (Ac * fcd)
    """
    NEd_N = abs(NEd_kN) * 1000.0
    if Ac_mm2 * fcd <= 0:
        return 1.0

    n = NEd_N / (Ac_mm2 * fcd)
    omega = Ast_mm2 * fyd / (Ac_mm2 * fcd)
    nu = 1.0 + omega
    n_bal = 0.4

    if abs(nu - n_bal) < 1e-9:
        return 1.0

    Kr = (nu - n) / (nu - n_bal)
    return min(max(Kr, 0.0), 1.0)


def _compute_Kphi(
    phi_eff: float,
    fck: float,
    lambda_: float,
) -> float:
    """Compute creep correction factor Kphi per IRC 112 Cl. 8.3.2.2.

    Kphi = 1 + beta * phi_eff
    beta = 0.35 + fck/200 - lambda/150
    """
    beta = 0.35 + fck / 200.0 - lambda_ / 150.0
    beta = max(beta, 0.0)
    return 1.0 + beta * phi_eff


def _sls_stresses(
    P_kN: float,
    M_kNm: float,
    Ac_mm2: float,
    I_mm4: float,
    y_max_mm: float,
    Ast_mm2: float,
    modular_ratio: float,
) -> tuple[float, float]:
    """Estimate concrete and steel stresses under SLS using elastic analysis.

    Assumes uncracked transformed section for a simplified estimate.

    Returns (sigma_c, sigma_s) both in MPa.
    sigma_c is the maximum compressive fibre stress.
    sigma_s is the stress in the extreme reinforcement.
    """
    N = abs(P_kN) * 1000.0  # N
    M = abs(M_kNm) * 1e6     # N.mm

    # Transformed section properties (approximate: add (m-1)*As to Ac)
    m = modular_ratio
    A_tr = Ac_mm2 + (m - 1.0) * Ast_mm2
    # For inertia, approximate bars on the perimeter at ~0.8 * y_max
    r_bar = 0.8 * y_max_mm
    I_tr = I_mm4 + (m - 1.0) * Ast_mm2 * r_bar ** 2

    # Concrete stress at extreme fibre
    sigma_c = N / A_tr + M * y_max_mm / I_tr  # MPa (N/mm2)

    # Steel stress at outermost bar
    if Ast_mm2 > 0:
        sigma_s = m * (N / A_tr + M * r_bar / I_tr)
    else:
        sigma_s = 0.0

    return abs(sigma_c), abs(sigma_s)


def _estimate_crack_width(
    sigma_s: float,
    bar_dia: float,
    cover: float,
    fctm: float,
    rho_eff: float,
    Es_GPa: float,
) -> float:
    """Estimate crack width per IRC 112 Cl. 12.3.4.

    Simplified formula:
        w_k = s_r_max * (epsilon_sm - epsilon_cm)

    s_r_max = 3.4 * c + 0.425 * k1 * k2 * phi / rho_eff
        k1 = 0.8 (deformed bars)
        k2 = 0.5 (bending)

    (epsilon_sm - epsilon_cm) = max(
        [sigma_s - kt * fctm / rho_eff * (1 + modular_ratio * rho_eff)] / Es,
        0.6 * sigma_s / Es
    )
    kt = 0.4 for long-term loading
    """
    if rho_eff <= 0 or sigma_s <= 0:
        return 0.0

    Es_MPa = Es_GPa * 1000.0
    k1 = 0.8
    k2 = 0.5
    kt = 0.4

    s_r_max = 3.4 * cover + 0.425 * k1 * k2 * bar_dia / rho_eff

    # Strain difference
    alpha_e = Es_GPa / 35.0  # approximate modular ratio for strain calc
    term1 = (sigma_s - kt * fctm / rho_eff * (1.0 + alpha_e * rho_eff)) / Es_MPa
    term2 = 0.6 * sigma_s / Es_MPa
    eps_diff = max(term1, term2)

    return s_r_max * eps_diff


def _ductile_detailing(
    pier_type: str,
    diameter_mm: float,
    width_trans_mm: float,
    pier_height_m: float,
    bar_dia: float,
    n_bars: int,
    Ac_mm2: float,
    fck: float,
    fcd: float,
    fyd: float,
    Ast_mm2: float,
    NEd_kN: float,
    config: dict[str, Any],
) -> tuple[float, float, float, float, float, bool]:
    """Compute ductile detailing requirements per IRC 112 Cl. 17.2.1.

    Returns (Lp, spiral_dia, spiral_spacing, omega_wd, omega_wd_required, ok).

    Plastic hinge length: Lp = max(D, L/6, 450 mm)
    Maximum spacing within plastic hinge: s_max = min(D/5, 6*phi_L, 150 mm)
    Minimum confining reinforcement ratio: omega_wd_min = 0.12
    Required: omega_w_req = 0.37 * Ac/Acc * nk + 0.13 * fyd/fcd * (rho_L - 0.01)
    """
    L_mm = pier_height_m * 1000.0
    D = diameter_mm if pier_type == "Circular" else max(diameter_mm, width_trans_mm)

    # Plastic hinge length
    Lp = max(D, L_mm / 6.0, 450.0)

    # Maximum spacing of confining reinforcement in plastic hinge
    s_max = min(D / 5.0, 6.0 * bar_dia, 150.0)

    # Default spiral/hoop diameter: typically 12 mm or 16 mm
    pier_cfg = config.get("pier", {})
    spiral_dia = float(pier_cfg.get("spiral_dia", 12.0))

    # Use s_max as the provided spacing (conservative design)
    spiral_spacing = s_max

    # Confined concrete area (core within hoops)
    cover = _get_cover(config)
    if pier_type == "Circular":
        # Confined core diameter = D - 2 * (cover + spiral_dia)
        Dc = D - 2.0 * (cover + spiral_dia)
        Acc = math.pi * Dc ** 2 / 4.0  # mm2
    else:
        # Rectangular: confined core dimensions
        bc = diameter_mm - 2.0 * (cover + spiral_dia)
        dc = width_trans_mm - 2.0 * (cover + spiral_dia)
        bc = max(bc, 1.0)
        dc = max(dc, 1.0)
        Acc = bc * dc

    Acc = max(Acc, 1.0)

    # Normalised axial force ratio
    NEd_N = abs(NEd_kN) * 1000.0
    nk = NEd_N / (Ac_mm2 * fck) if (Ac_mm2 * fck) > 0 else 0.0

    # Longitudinal reinforcement ratio
    rho_L = Ast_mm2 / Ac_mm2 if Ac_mm2 > 0 else 0.0

    # Required confining reinforcement ratio
    omega_wd_min = 0.12
    omega_w_req = (
        0.37 * (Ac_mm2 / Acc) * nk
        + 0.13 * (fyd / fcd) * max(rho_L - 0.01, 0.0)
    )
    omega_wd_required = max(omega_wd_min, omega_w_req)

    # Auto-size confining reinforcement to satisfy omega_wd_required.
    # Try standard spiral diameters in increasing order until the required
    # spacing is practical (>= 75 mm).
    _MIN_SPIRAL_SPACING = 75.0  # mm -- practical minimum
    _CANDIDATE_DIAS = [spiral_dia, 12.0, 16.0, 20.0, 25.0, 32.0]
    # Deduplicate while preserving order
    _seen: set[float] = set()
    candidate_dias: list[float] = []
    for _d in _CANDIDATE_DIAS:
        if _d not in _seen:
            _seen.add(_d)
            candidate_dias.append(_d)

    for try_dia in candidate_dias:
        Ash_try = _bar_area(try_dia)

        if pier_type == "Circular":
            Dc_try = D - 2.0 * (cover + try_dia)
            Dc_try = max(Dc_try, 1.0)
            Acc_try = math.pi * Dc_try ** 2 / 4.0
        else:
            bc_try = diameter_mm - 2.0 * (cover + try_dia)
            dc_try = width_trans_mm - 2.0 * (cover + try_dia)
            bc_try = max(bc_try, 1.0)
            dc_try = max(dc_try, 1.0)
            Acc_try = bc_try * dc_try
        Acc_try = max(Acc_try, 1.0)

        # Back-calculate the spacing needed for omega_wd = omega_wd_required
        denom = omega_wd_required * Acc_try * fcd
        if pier_type == "Circular":
            s_needed = (Ash_try * math.pi * Dc_try * fyd / denom) if denom > 0 else 0.0
        else:
            n_legs = 2
            s_needed = (n_legs * Ash_try * (bc_try + dc_try) * fyd / denom) if denom > 0 else 0.0

        if s_needed >= _MIN_SPIRAL_SPACING:
            # This diameter works -- use min(s_needed, s_max) as the spacing
            spiral_dia = try_dia
            spiral_spacing = min(s_needed, s_max)
            # Round spacing down to nearest 5 mm for practicality
            spiral_spacing = max(math.floor(spiral_spacing / 5.0) * 5.0, _MIN_SPIRAL_SPACING)
            Acc = Acc_try
            if pier_type == "Circular":
                Dc = Dc_try
            break
    # else: keep original spiral_dia / s_max (fallback -- largest dia still insufficient)

    # Compute provided omega_wd with final dia and spacing
    Ash = _bar_area(spiral_dia)
    if pier_type == "Circular":
        if Acc * spiral_spacing * fcd > 0:
            omega_wd = (Ash * math.pi * max(Dc, 0.0) * fyd) / (Acc * spiral_spacing * fcd)
        else:
            omega_wd = 0.0
    else:
        n_legs = 2
        if Acc * spiral_spacing * fcd > 0:
            omega_wd = (
                n_legs * Ash * (max(bc, 0.0) + max(dc, 0.0)) * fyd
                / (Acc * spiral_spacing * fcd)
            )
        else:
            omega_wd = 0.0

    ok = omega_wd >= omega_wd_required

    return Lp, spiral_dia, spiral_spacing, omega_wd, omega_wd_required, ok


def _find_governing_uls(
    combinations: CombinationResults,
) -> CombinationResult:
    """Find the governing ULS combination for pier design.

    Selects the combination with the maximum resultant moment
    sqrt(ML^2 + MT^2) from both ULS basic and ULS seismic categories.
    """
    uls_combos = [
        c for c in combinations.all_combinations
        if c.category in ("uls_basic", "uls_seismic")
    ]
    if not uls_combos:
        # Fallback to any combination
        uls_combos = combinations.all_combinations

    return max(
        uls_combos,
        key=lambda c: math.sqrt(
            c.forces_pier_base.ML ** 2 + c.forces_pier_base.MT ** 2
        ),
    )


def _find_max_axial_uls(
    combinations: CombinationResults,
) -> CombinationResult:
    """Find the ULS combination with maximum axial load P."""
    uls_combos = [
        c for c in combinations.all_combinations
        if c.category in ("uls_basic", "uls_seismic")
    ]
    if not uls_combos:
        uls_combos = combinations.all_combinations
    return max(uls_combos, key=lambda c: c.forces_pier_base.P)


def _find_sls_rare_governing(
    combinations: CombinationResults,
) -> CombinationResult:
    """Find governing SLS Rare combination (max resultant moment)."""
    return combinations.governing_sls_rare


def _find_sls_qp_governing(
    combinations: CombinationResults,
) -> CombinationResult:
    """Find governing SLS Quasi-Permanent combination."""
    return combinations.governing_sls_quasi


def _resultant_moment(fv: ForceVector) -> float:
    """Compute sqrt(ML^2 + MT^2) from a force vector."""
    return math.sqrt(fv.ML ** 2 + fv.MT ** 2)


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def design_pier(
    config: dict[str, Any],
    geometry: GeometryResults,
    combinations: CombinationResults,
) -> PierDesignResult:
    """Design the pier column under biaxial bending per IRC 112:2020.

    This is the primary public API of the pier design module.  It reads
    section geometry and material properties from the configuration,
    determines reinforcement, checks slenderness and second-order effects,
    generates P-M interaction diagrams, verifies ULS biaxial utilisation,
    checks SLS stress and crack width limits, and assesses ductile
    detailing requirements.

    Parameters
    ----------
    config : dict
        Validated configuration dictionary from :func:`input_parser.parse_input`.
    geometry : GeometryResults
        Pre-computed geometry from :func:`geometry.calculate_geometry`.
    combinations : CombinationResults
        Load combination results from :func:`load_combinations.calculate_combinations`.

    Returns
    -------
    PierDesignResult
        Comprehensive design results including status, utilisations,
        interaction diagrams, and detailing requirements.
    """
    # ==================================================================
    # 1. Section properties and materials
    # ==================================================================
    pier_cfg = config["pier"]
    mat_cfg = config["materials"]
    pier_section = geometry.pier_section
    pier_type = pier_section.pier_type

    fck = float(mat_cfg["pier"]["fck"])
    aggregate = mat_cfg["pier"].get("aggregate", "Quartzite")
    exposure = mat_cfg.get("exposure", "Moderate")

    concrete = get_concrete_properties(fck, aggregate=aggregate, exposure=exposure)
    fcd = concrete.fcd   # MPa
    fctm = concrete.fctm

    fyk = float(mat_cfg["steel"]["fyk"])
    Es_input = mat_cfg["steel"].get("Es", 200000)
    # Es might be in MPa or GPa; normalise to GPa
    Es_GPa = float(Es_input)
    if Es_GPa > 1000:
        Es_GPa = Es_GPa / 1000.0  # was in MPa, convert to GPa
    steel = get_steel_properties(fyk=fyk, Es=Es_GPa)
    fyd = steel.fyd

    cover = _get_cover(config)

    # Section dimensions in mm
    diameter_mm = pier_section.width_long * 1000.0   # m -> mm
    width_trans_mm = pier_section.width_trans * 1000.0

    Ac_mm2 = _compute_Ac_mm2(pier_section)
    Ixx_mm4 = _compute_Ic_mm4(pier_section, "xx")  # about longitudinal axis
    Iyy_mm4 = _compute_Ic_mm4(pier_section, "yy")  # about transverse axis

    pier_height_m = geometry.pier_height

    # ==================================================================
    # 2. Reinforcement selection
    # ==================================================================
    # Find maximum axial load among ULS combos for minimum reinforcement check
    max_axial_combo = _find_max_axial_uls(combinations)
    max_NEd = max_axial_combo.forces_pier_base.P  # kN

    n_bars, bar_dia, Ast_mm2 = _select_reinforcement(
        config, Ac_mm2, fcd, fyd, max_NEd,
    )
    rho_l = Ast_mm2 / Ac_mm2 if Ac_mm2 > 0 else 0.0

    # ==================================================================
    # 3. Creep coefficient
    # ==================================================================
    phi_creep = _compute_creep_coefficient(config, concrete, pier_section)

    # Effective creep ratio: needs quasi-permanent and ULS moments
    gov_uls = _find_governing_uls(combinations)
    gov_qp = _find_sls_qp_governing(combinations)

    M0_Ed = _resultant_moment(gov_uls.forces_pier_base)
    M0_qp = _resultant_moment(gov_qp.forces_pier_base)
    phi_eff = _effective_creep_ratio(phi_creep, M0_qp, M0_Ed)

    # ==================================================================
    # 4. Slenderness check (both axes)
    # ==================================================================
    l0_mm = _effective_length(config, pier_height_m)

    # End moments for slenderness formula: use governing ULS moments
    # M01 = smaller end moment, M02 = larger end moment
    # For cantilever: M at base = M02, M at top = 0 -> M01 = 0
    end_condition = pier_cfg.get("end_condition", "fixed-free")
    if end_condition in ("fixed-free", "fixed-guided"):
        # Cantilever or guided top: moment at base, zero at top
        M01_xx = 0.0
        M02_xx = abs(gov_uls.forces_pier_base.ML)
        M01_yy = 0.0
        M02_yy = abs(gov_uls.forces_pier_base.MT)
    else:
        # Fixed-fixed: assume symmetric moments, rm = -1
        M01_xx = abs(gov_uls.forces_pier_base.ML)
        M02_xx = abs(gov_uls.forces_pier_base.ML)
        M01_yy = abs(gov_uls.forces_pier_base.MT)
        M02_yy = abs(gov_uls.forces_pier_base.MT)

    NEd = gov_uls.forces_pier_base.P  # kN

    # Check about xx (transverse axis -- moment ML bends about xx)
    lambda_xx, lambda_lim_xx, second_order_xx = _slenderness_check(
        l0_mm, Ixx_mm4, Ac_mm2, fck, fcd, fyd, Ast_mm2,
        NEd, phi_eff, M01_xx, M02_xx,
    )
    # Check about yy (longitudinal axis -- moment MT bends about yy)
    lambda_yy, lambda_lim_yy, second_order_yy = _slenderness_check(
        l0_mm, Iyy_mm4, Ac_mm2, fck, fcd, fyd, Ast_mm2,
        NEd, phi_eff, M01_yy, M02_yy,
    )

    # Governing slenderness (report the more critical axis)
    if lambda_xx >= lambda_yy:
        lambda_report = lambda_xx
        lambda_lim_report = lambda_lim_xx
    else:
        lambda_report = lambda_yy
        lambda_lim_report = lambda_lim_yy

    second_order_needed = second_order_xx or second_order_yy

    # ==================================================================
    # 5. Second-order effects (if needed)
    # ==================================================================
    ML_applied = gov_uls.forces_pier_base.ML  # kN.m
    MT_applied = gov_uls.forces_pier_base.MT  # kN.m
    P_applied = gov_uls.forces_pier_base.P    # kN

    M2_xx = 0.0
    M2_yy = 0.0

    if second_order_needed:
        # Effective depth for curvature calculation
        if pier_type == "Circular":
            # d = D/2 + r_bar_circle for circular, approximate as 0.9 * D/2
            d_xx = 0.45 * diameter_mm + (diameter_mm / 2.0 - cover - bar_dia / 2.0)
            d_xx = max(d_xx, 0.5 * diameter_mm)
            d_yy = d_xx  # symmetric for circular
        else:
            d_xx = width_trans_mm - cover - bar_dia / 2.0
            d_yy = diameter_mm - cover - bar_dia / 2.0

        Kr = _compute_Kr(NEd, Ac_mm2, fcd, fyd, Ast_mm2)

        if second_order_xx:
            Kphi_xx = _compute_Kphi(phi_eff, fck, lambda_xx)
            M2_xx = _second_order_moment(
                NEd, l0_mm, fyd, Es_GPa, d_xx, Kr, Kphi_xx,
            )

        if second_order_yy:
            Kphi_yy = _compute_Kphi(phi_eff, fck, lambda_yy)
            M2_yy = _second_order_moment(
                NEd, l0_mm, fyd, Es_GPa, d_yy, Kr, Kphi_yy,
            )

    ML_design = abs(ML_applied) + M2_xx  # total design moment about xx
    MT_design = abs(MT_applied) + M2_yy  # total design moment about yy

    # Ensure minimum eccentricity: e_min = max(h/30, 20 mm) per IRC 112
    e_min_xx = max(width_trans_mm / 30.0, 20.0)  # mm
    e_min_yy = max(diameter_mm / 30.0, 20.0)
    M_min_xx = abs(P_applied) * e_min_xx / 1000.0  # kN.m
    M_min_yy = abs(P_applied) * e_min_yy / 1000.0
    ML_design = max(ML_design, M_min_xx)
    MT_design = max(MT_design, M_min_yy)

    # ==================================================================
    # 6. Generate interaction diagrams
    # ==================================================================
    if pier_type == "Circular":
        interaction_xx = generate_interaction_circular(
            diameter=diameter_mm,
            n_bars=n_bars,
            bar_diameter=bar_dia,
            cover=cover,
            fck=fck,
            fyk=fyk,
        )
        # For circular, xx and yy are the same
        interaction_yy = interaction_xx
    else:
        # Rectangular: xx axis (bending about longitudinal axis,
        # width = width_long, depth = width_trans)
        # For rectangular pier, split bars equally between tension/compression faces
        n_bars_face = max(n_bars // 2, 2)
        interaction_xx = generate_interaction_rectangular(
            width=diameter_mm,       # width_long in mm
            depth=width_trans_mm,    # width_trans in mm
            n_bars_tension=n_bars_face,
            n_bars_compression=n_bars_face,
            bar_diameter=bar_dia,
            cover=cover,
            fck=fck,
            fyk=fyk,
        )
        # yy axis (bending about transverse axis,
        # width = width_trans, depth = width_long)
        interaction_yy = generate_interaction_rectangular(
            width=width_trans_mm,
            depth=diameter_mm,
            n_bars_tension=n_bars_face,
            n_bars_compression=n_bars_face,
            bar_diameter=bar_dia,
            cover=cover,
            fck=fck,
            fyk=fyk,
        )

    # Extract capacities
    P_capacity = interaction_xx.P_max  # kN (pure compression)
    M_capacity_xx = interaction_xx.M_max  # kN.m at max moment (near balanced)
    M_capacity_yy = interaction_yy.M_max

    # ==================================================================
    # 7. ULS utilisation check — evaluate ALL ULS combinations
    # ==================================================================
    # Evaluate every ULS combination, not just the max-SRSS one, because
    # a combo with lower P (relieving DL) may govern due to reduced
    # section capacity at low axial loads (proportional-scaling effect).
    uls_combos = [
        c for c in combinations.all_combinations
        if c.category in ("uls_basic", "uls_seismic")
    ]
    if not uls_combos:
        uls_combos = [gov_uls]

    # Pre-compute constants for per-combo second-order effects
    if pier_type == "Circular":
        _d_so_xx = 0.45 * diameter_mm + (diameter_mm / 2.0 - cover - bar_dia / 2.0)
        _d_so_xx = max(_d_so_xx, 0.5 * diameter_mm)
        _d_so_yy = _d_so_xx
    else:
        _d_so_xx = width_trans_mm - cover - bar_dia / 2.0
        _d_so_yy = diameter_mm - cover - bar_dia / 2.0
    _e_min_xx = max(width_trans_mm / 30.0, 20.0)
    _e_min_yy = max(diameter_mm / 30.0, 20.0)

    util_biaxial = -1.0
    util_xx = 0.0
    util_yy = 0.0
    gov_combo_name = gov_uls.name

    for _combo in uls_combos:
        _fv = _combo.forces_pier_base
        _P = _fv.P
        _ML_abs = abs(_fv.ML)
        _MT_abs = abs(_fv.MT)

        # Slenderness & second-order determination for this combo
        if end_condition in ("fixed-free", "fixed-guided"):
            _M01_xx, _M02_xx = 0.0, _ML_abs
            _M01_yy, _M02_yy = 0.0, _MT_abs
        else:
            _M01_xx = _M02_xx = _ML_abs
            _M01_yy = _M02_yy = _MT_abs

        _, _, _so_xx = _slenderness_check(
            l0_mm, Ixx_mm4, Ac_mm2, fck, fcd, fyd, Ast_mm2,
            _P, phi_eff, _M01_xx, _M02_xx,
        )
        _, _, _so_yy = _slenderness_check(
            l0_mm, Iyy_mm4, Ac_mm2, fck, fcd, fyd, Ast_mm2,
            _P, phi_eff, _M01_yy, _M02_yy,
        )

        _M2_xx = 0.0
        _M2_yy = 0.0
        if _so_xx or _so_yy:
            _Kr = _compute_Kr(_P, Ac_mm2, fcd, fyd, Ast_mm2)
            if _so_xx:
                _Kphi = _compute_Kphi(phi_eff, fck, lambda_xx)
                _M2_xx = _second_order_moment(
                    _P, l0_mm, fyd, Es_GPa, _d_so_xx, _Kr, _Kphi,
                )
            if _so_yy:
                _Kphi = _compute_Kphi(phi_eff, fck, lambda_yy)
                _M2_yy = _second_order_moment(
                    _P, l0_mm, fyd, Es_GPa, _d_so_yy, _Kr, _Kphi,
                )

        _ML_des = max(_ML_abs + _M2_xx, abs(_P) * _e_min_xx / 1000.0)
        _MT_des = max(_MT_abs + _M2_yy, abs(_P) * _e_min_yy / 1000.0)

        # Utilization
        if pier_type == "Circular":
            _M_res = math.sqrt(_ML_des ** 2 + _MT_des ** 2)
            _util = check_utilisation(interaction_xx, _P, _M_res)
        else:
            _util = check_biaxial(
                interaction_xx, interaction_yy,
                _P, _ML_des, _MT_des,
            )

        if _util > util_biaxial:
            util_biaxial = _util
            P_applied = _P
            ML_applied = _fv.ML
            MT_applied = _fv.MT
            ML_design = _ML_des
            MT_design = _MT_des
            util_xx = check_utilisation(interaction_xx, _P, _ML_des)
            util_yy = check_utilisation(interaction_yy, _P, _MT_des)
            gov_combo_name = _combo.name

    # ==================================================================
    # 8. SLS checks
    # ==================================================================
    # SLS Rare -- stress limits
    gov_sls_rare = _find_sls_rare_governing(combinations)
    P_sls_rare = gov_sls_rare.forces_pier_base.P
    M_sls_rare = _resultant_moment(gov_sls_rare.forces_pier_base)

    # Use y_max as half the section depth for the axis of max moment
    y_max_xx = width_trans_mm / 2.0
    y_max_yy = diameter_mm / 2.0
    y_max = max(y_max_xx, y_max_yy)

    m_ratio = concrete.modular_ratio_short

    sigma_c_rare, sigma_s_rare = _sls_stresses(
        P_sls_rare, M_sls_rare,
        Ac_mm2, max(Ixx_mm4, Iyy_mm4), y_max,
        Ast_mm2, m_ratio,
    )

    # SLS Quasi-Permanent -- stress and crack width
    P_sls_qp = gov_qp.forces_pier_base.P
    M_sls_qp = _resultant_moment(gov_qp.forces_pier_base)

    sigma_c_qp, sigma_s_qp = _sls_stresses(
        P_sls_qp, M_sls_qp,
        Ac_mm2, max(Ixx_mm4, Iyy_mm4), y_max,
        Ast_mm2, m_ratio,
    )

    # Stress limits per IRC 112
    sigma_c_rare_limit = 0.48 * fck   # Rare combination
    sigma_c_qp_limit = 0.36 * fck     # Quasi-permanent
    sigma_s_limit = 0.8 * fyk         # Steel stress under rare

    # Effective reinforcement ratio for crack width
    # A_ct_eff = effective tension area ~ 0.5 * Ac for columns
    A_ct_eff = 0.5 * Ac_mm2
    rho_eff = Ast_mm2 / A_ct_eff if A_ct_eff > 0 else 0.01

    crack_width = _estimate_crack_width(
        sigma_s_qp, bar_dia, cover, fctm, rho_eff, Es_GPa,
    )

    # Crack width limit from IRC tables
    irc_tables = load_irc_tables()
    crack_limits = irc_tables.get("crack_width_limits", {})
    crack_limit = float(crack_limits.get(exposure, 0.3))

    sls_ok = (
        sigma_c_rare <= sigma_c_rare_limit
        and sigma_c_qp <= sigma_c_qp_limit
        and sigma_s_rare <= sigma_s_limit
        and crack_width <= crack_limit
    )

    # ==================================================================
    # 9. Ductile detailing
    # ==================================================================
    Lp, spiral_dia, spiral_spacing, omega_wd, omega_wd_req, ductile_ok = (
        _ductile_detailing(
            pier_type, diameter_mm, width_trans_mm,
            pier_height_m, bar_dia, n_bars,
            Ac_mm2, fck, fcd, fyd, Ast_mm2,
            P_applied, config,
        )
    )

    # ==================================================================
    # 10. Overall status
    # ==================================================================
    uls_ok = util_biaxial <= 1.0
    status = "OK" if (uls_ok and sls_ok and ductile_ok) else "NOT OK"

    # ==================================================================
    # Assemble result
    # ==================================================================
    return PierDesignResult(
        # Section
        pier_type=pier_type,
        diameter=diameter_mm,
        width_trans=width_trans_mm,
        height=pier_height_m,
        fck=fck,
        cover=cover,
        # Reinforcement
        n_bars=n_bars,
        bar_dia=bar_dia,
        Ast_provided=Ast_mm2,
        rho_l=rho_l,
        # Slenderness
        lambda_=lambda_report,
        lambda_lim=lambda_lim_report,
        second_order=second_order_needed,
        # ULS capacity
        P_capacity=P_capacity,
        M_capacity_xx=M_capacity_xx,
        M_capacity_yy=M_capacity_yy,
        # Governing combination
        governing_combo_name=gov_combo_name,
        P_applied=P_applied,
        ML_applied=ML_applied,
        MT_applied=MT_applied,
        ML_with_2nd_order=ML_design,
        MT_with_2nd_order=MT_design,
        # Utilisation
        util_uniaxial_xx=util_xx,
        util_uniaxial_yy=util_yy,
        util_biaxial=util_biaxial,
        # SLS
        sigma_c_rare=sigma_c_rare,
        sigma_c_qp=sigma_c_qp,
        sigma_s_rare=sigma_s_rare,
        crack_width=crack_width,
        sls_ok=sls_ok,
        # Ductile detailing
        Lp=Lp,
        spiral_dia=spiral_dia,
        spiral_spacing=spiral_spacing,
        omega_wd=omega_wd,
        omega_wd_required=omega_wd_req,
        ductile_ok=ductile_ok,
        # Overall
        status=status,
        # Interaction diagrams
        interaction_xx=interaction_xx,
        interaction_yy=interaction_yy,
    )

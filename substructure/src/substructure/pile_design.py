"""
Pile structural design: bending, shear, and crack width per IRC 112:2020.

This module designs the pile structurally using:
- Interaction diagram for flexure (P-M)
- IRC 112 shear provisions for circular sections
- IRC 112 crack width calculations for SLS
"""

from dataclasses import dataclass
from typing import Any
import math

from .loads import ForceVector
from .geometry import GeometryResults, PileProperties
from .load_combinations import CombinationResult, CombinationResults
from .materials import (
    get_concrete_properties,
    get_steel_properties,
    ConcreteProperties,
    SteelProperties,
)
from .interaction import (
    generate_interaction_circular,
    InteractionDiagram,
    check_utilisation,
)
from .utils import load_irc_tables
from .pile_capacity import PileCapacityResult, PileLoad, PileCombinationResult


@dataclass
class PileDesignResult:
    """Results from pile structural design."""

    diameter: float  # mm
    cover: float  # mm
    d_eff: float  # mm
    fck: float  # MPa

    n_bars: int
    bar_dia: float  # mm
    Ast_provided: float  # mm2
    rho_l: float

    governing_combo_name: str
    P_applied: float  # kN
    H_applied: float  # kN (resultant horizontal)
    M_applied: float  # kN.m

    interaction: InteractionDiagram
    util_PM: float

    VEd: float  # kN
    VRdc: float  # kN
    shear_util: float

    crack_width: float  # mm
    crack_width_limit: float  # mm
    crack_ok: bool

    curtailment_depth: float  # m (where moment reduces enough to curtail bars)

    status: str


def design_pile(
    config: dict[str, Any],
    geometry: GeometryResults,
    combinations: CombinationResults,
    pile_capacity_result: PileCapacityResult,
) -> PileDesignResult:
    """
    Design pile structurally for bending, shear, and crack width.

    Args:
        config: Configuration dictionary
        geometry: Geometry results
        combinations: Load combinations
        pile_capacity_result: Pile capacity and load distribution results

    Returns:
        PileDesignResult with design checks
    """
    # Step 1: Material properties
    fck = config["materials"]["pile_pilecap"]["fck"]  # MPa
    fyk = config["materials"]["steel"]["fyk"]  # MPa
    Es = config["materials"]["steel"]["Es"]  # MPa

    aggregate = config["materials"].get("pile_pilecap", {}).get("aggregate", "Quartzite")
    exposure = config["materials"].get("exposure", "Moderate")
    concrete = get_concrete_properties(fck, aggregate=aggregate, exposure=exposure)
    Es_GPa = Es / 1000.0 if Es > 1000 else Es
    steel = get_steel_properties(fyk, Es_GPa)

    # Step 2: Pile dimensions
    D_m = geometry.pile_props.diameter  # m
    D_mm = D_m * 1000  # mm
    cover = config["foundation"]["cover"]  # mm (already in mm!)

    # Get reinforcement from config or use default
    fnd_cfg = config.get("foundation", {})
    n_bars = fnd_cfg.get("pile_n_bars") or 12
    bar_dia = fnd_cfg.get("pile_bar_dia") or 25.0  # mm

    Ast_provided = n_bars * math.pi * (bar_dia / 2) ** 2  # mm2
    Ac = math.pi * (D_mm / 2) ** 2  # mm2

    # Effective depth (conservative: to centroid of tension steel)
    # For circular section, approx d_eff = D - cover - bar_dia/2
    d_eff = D_mm - cover - bar_dia / 2  # mm

    rho_l = min(Ast_provided / (D_mm * d_eff), 0.02)

    # Step 3: Check minimum reinforcement
    fyd = fyk / 1.15  # MPa
    # We'll check this after finding governing load

    # Step 4: Find governing pile load (ULS)
    # Look through all ULS combinations, find pile with max H
    max_H = 0.0
    governing_combo_name = ""
    governing_P = 0.0
    governing_H = 0.0
    governing_HL = 0.0
    governing_HT = 0.0

    for combo_result in pile_capacity_result.all_pile_results:
        # Check if this is ULS combo
        combo_name = combo_result.combo_name
        is_uls = any(
            x in combo_name.upper()
            for x in ["BASIC", "SEISMIC", "ACCIDENTAL", "ULS"]
        )

        if not is_uls:
            continue

        # Find pile with max H in this combo
        for pile_load in combo_result.pile_loads:
            HL = pile_load.HL
            HT = pile_load.HT
            H = math.sqrt(HL**2 + HT**2)

            if H > max_H:
                max_H = H
                governing_combo_name = combo_name
                governing_P = pile_load.P
                governing_H = H
                governing_HL = HL
                governing_HT = HT

    # If no ULS combo found, use the first combo
    if governing_combo_name == "":
        combo_result = pile_capacity_result.all_pile_results[0]
        governing_combo_name = combo_result.combo_name
        pile_load = combo_result.pile_loads[0]
        governing_P = pile_load.P
        governing_HL = pile_load.HL
        governing_HT = pile_load.HT
        governing_H = math.sqrt(governing_HL**2 + governing_HT**2)

    # Check minimum reinforcement with governing load
    As_min_1 = 0.002 * Ac  # mm2
    As_min_2 = 0.1 * abs(governing_P) * 1000 / fyd  # kN -> N, then / MPa = mm2
    As_min = max(As_min_1, As_min_2)

    if Ast_provided < As_min:
        # Increase reinforcement (adjust number of bars)
        n_bars = math.ceil(As_min / (math.pi * (bar_dia / 2) ** 2))
        Ast_provided = n_bars * math.pi * (bar_dia / 2) ** 2
        rho_l = min(Ast_provided / (D_mm * d_eff), 0.02)

    # Step 5: Compute applied moment (IS:2911 coefficient method)
    # For a fixed-head pile (embedded in pile cap), the maximum moment at
    # the pile head is M = Am * H * T where Am is the IS:2911 moment
    # coefficient.  For fixed-head piles with zero free length (e/T = 0),
    # Am ≈ 0.93.  This replaces the overly-conservative M = H * Lf.
    T = geometry.pile_props.stiffness_factor_T  # m
    Am = config.get("foundation", {}).get("pile_moment_coeff", 0.93)
    M_applied = Am * governing_H * T  # kN.m

    # Step 6: Generate interaction diagram and check utilisation
    interaction = generate_interaction_circular(
        diameter=D_mm,
        n_bars=n_bars,
        bar_diameter=bar_dia,
        cover=cover,
        fck=fck,
        fyk=fyk,
    )

    P_kN = governing_P  # kN (compression positive)
    M_kNm = M_applied  # kN.m

    util_PM = check_utilisation(interaction, P_kN, M_kNm)

    # Step 7: Shear check (IRC 112 Cl. 10.3.2)
    VEd = governing_H  # kN (shear at fixity point equals horizontal load)

    # For circular section, use equivalent rectangular
    bw = D_mm  # mm
    d = d_eff  # mm

    k = min(1.0 + math.sqrt(200.0 / d), 2.0)

    # Axial stress (compression positive in MPa)
    sigma_cp = min(
        abs(P_kN) * 1000 / Ac,  # N/mm2 = MPa
        0.2 * concrete.fcd,
    )

    # VRd,c formulas
    term1 = 0.12 * k * (80 * rho_l * fck) ** (1.0 / 3.0) + 0.15 * sigma_cp
    term2 = 0.031 * k ** (3.0 / 2.0) * math.sqrt(fck) + 0.15 * sigma_cp

    VRdc_1 = term1 * bw * d / 1000  # kN
    VRdc_2 = term2 * bw * d / 1000  # kN

    VRdc = max(VRdc_1, VRdc_2)  # kN

    shear_util = abs(VEd) / VRdc if VRdc > 0 else 999.0

    # Step 8: Crack width check (SLS)
    # Use quasi-permanent SLS combination per IRC 112 Cl. 12.3.4
    max_H_sls = 0.0
    sls_P = 0.0
    sls_H = 0.0

    for combo_result in pile_capacity_result.all_pile_results:
        # Filter for quasi-permanent SLS combinations
        if combo_result.category != "sls_quasi_permanent":
            continue

        for pile_load in combo_result.pile_loads:
            HL = pile_load.HL
            HT = pile_load.HT
            H = math.sqrt(HL**2 + HT**2)

            if H > max_H_sls:
                max_H_sls = H
                sls_P = pile_load.P
                sls_H = H

    # If no quasi-permanent combo found, try any SLS combo, then fallback
    if max_H_sls == 0.0:
        for combo_result in pile_capacity_result.all_pile_results:
            if "sls" in combo_result.category:
                for pile_load in combo_result.pile_loads:
                    H = math.sqrt(pile_load.HL**2 + pile_load.HT**2)
                    if H > max_H_sls:
                        max_H_sls = H
                        sls_P = pile_load.P
                        sls_H = H

    # If still no SLS combo found, use 60% of ULS (quasi-permanent fraction)
    if max_H_sls == 0.0:
        sls_P = 0.6 * governing_P
        sls_H = 0.6 * governing_H

    M_sls = Am * sls_H * T  # kN.m (same IS:2911 formula as ULS)

    # Calculate crack width
    crack_width = calculate_crack_width(
        M_kNm=M_sls,
        P_kN=sls_P,
        D_mm=D_mm,
        cover=cover,
        n_bars=n_bars,
        bar_dia=bar_dia,
        Ast=Ast_provided,
        fck=fck,
        Es=Es,
    )

    crack_width_limit = 0.3  # mm for moderate exposure (IRC 112 Table 12.1)
    crack_ok = crack_width <= crack_width_limit

    # Step 9: Curtailment depth
    # Approximate as depth where moment = 50% of max
    # For fixed-head pile, moment decreases linearly, so curtailment at ~Lf/2
    Lf = geometry.pile_props.fixity_depth_Lf  # m
    curtailment_depth = Lf / 2.0  # m

    # Step 10: Overall status
    all_ok = util_PM <= 1.0 and shear_util <= 1.0 and crack_ok
    status = "OK" if all_ok else "NOT OK"

    return PileDesignResult(
        diameter=D_mm,
        cover=cover,
        d_eff=d_eff,
        fck=fck,
        n_bars=n_bars,
        bar_dia=bar_dia,
        Ast_provided=Ast_provided,
        rho_l=rho_l,
        governing_combo_name=governing_combo_name,
        P_applied=P_kN,
        H_applied=governing_H,
        M_applied=M_kNm,
        interaction=interaction,
        util_PM=util_PM,
        VEd=VEd,
        VRdc=VRdc,
        shear_util=shear_util,
        crack_width=crack_width,
        crack_width_limit=crack_width_limit,
        crack_ok=crack_ok,
        curtailment_depth=curtailment_depth,
        status=status,
    )


def calculate_crack_width(
    M_kNm: float,
    P_kN: float,
    D_mm: float,
    cover: float,
    n_bars: int,
    bar_dia: float,
    Ast: float,
    fck: float,
    Es: float,
) -> float:
    """
    Calculate crack width per IRC 112 Cl. 12.3.4.

    Args:
        M_kNm: Applied moment (kN.m)
        P_kN: Applied axial load (kN, compression positive)
        D_mm: Pile diameter (mm)
        cover: Concrete cover (mm)
        n_bars: Number of reinforcement bars
        bar_dia: Bar diameter (mm)
        Ast: Total steel area (mm2)
        fck: Concrete strength (MPa)
        Es: Steel modulus (MPa)

    Returns:
        Crack width (mm)
    """
    # Cross-section properties
    Ac = math.pi * (D_mm / 2) ** 2  # mm2
    I = math.pi * (D_mm / 4) ** 4 / 4  # mm4
    y_max = D_mm / 2  # mm

    # Steel stress (simplified: assume tension at extreme fiber)
    # M/I * y + P/A (compression positive, so tension is negative P)
    M_Nmm = M_kNm * 1e6  # N.mm
    P_N = P_kN * 1000  # N

    # Extreme fiber stress in concrete
    sigma_c = M_Nmm * y_max / I - P_N / Ac  # N/mm2 = MPa (tension positive here)

    # Transform to steel stress (approximate)
    # For cracked section, neutral axis shifts
    # Simplified: sigma_s ≈ M / (Ast * lever_arm)
    # lever_arm ≈ 0.9 * D/2 for circular section
    lever_arm = 0.9 * D_mm / 2  # mm

    if abs(M_Nmm) > 1e-6:
        sigma_s = abs(M_Nmm) / (Ast * lever_arm)  # MPa
    else:
        sigma_s = 0.0

    # Limit to reasonable values
    sigma_s = min(sigma_s, 0.8 * 550)  # MPa (80% of fyk)

    # IRC 112 crack width formulas (Cl. 12.3.4)
    c = cover  # mm
    phi_bar = bar_dia  # mm

    k1 = 0.8  # deformed bars
    k2 = 0.5  # bending
    kt = 0.4  # long-term loading

    # Effective reinforcement ratio
    # For circular section, use effective tension area
    # Approximate: A_eff = 2.5 * cover * perimeter
    perimeter = math.pi * D_mm  # mm
    A_eff = 2.5 * c * perimeter  # mm2
    rho_eff = min(Ast / A_eff, 0.02) if A_eff > 0 else 0.01

    # Maximum crack spacing
    sr_max = 3.4 * c + 0.425 * k1 * k2 * phi_bar / rho_eff  # mm

    # Mean strain
    fctm = 0.3 * fck ** (2.0 / 3.0)  # MPa (IRC 112 Table 6.5)
    alpha_e = Es / (5000 * math.sqrt(fck))  # modular ratio (approx)

    term1 = sigma_s - kt * fctm / rho_eff * (1 + alpha_e * rho_eff)
    term1 = term1 / Es if sigma_s > 0 else 0.0

    term2 = 0.6 * sigma_s / Es if sigma_s > 0 else 0.0

    eps_sm_minus_eps_cm = max(term1, term2)

    # Crack width
    wk = sr_max * eps_sm_minus_eps_cm  # mm

    return max(wk, 0.0)


def summarise_pile_design(result: PileDesignResult) -> str:
    """
    Generate formatted summary of pile design results.

    Args:
        result: PileDesignResult object

    Returns:
        Formatted string summary
    """
    lines = []
    lines.append("=" * 70)
    lines.append("PILE STRUCTURAL DESIGN SUMMARY (IRC 112:2020)")
    lines.append("=" * 70)
    lines.append("")

    # Section properties
    lines.append("SECTION PROPERTIES:")
    lines.append(f"  Diameter:           {result.diameter:.0f} mm")
    lines.append(f"  Cover:              {result.cover:.0f} mm")
    lines.append(f"  Effective depth:    {result.d_eff:.0f} mm")
    lines.append(f"  Concrete grade:     M{result.fck:.0f}")
    lines.append("")

    # Reinforcement
    lines.append("REINFORCEMENT:")
    lines.append(f"  Configuration:      {result.n_bars} bars of {result.bar_dia:.0f}mm dia")
    lines.append(f"  Ast provided:       {result.Ast_provided:.0f} mm²")
    lines.append(f"  Reinforcement ratio: {result.rho_l:.4f}")
    lines.append("")

    # Applied loads
    lines.append("GOVERNING LOADS (ULS):")
    lines.append(f"  Combination:        {result.governing_combo_name}")
    lines.append(f"  Axial load P:       {result.P_applied:+.1f} kN")
    lines.append(f"  Horizontal H:       {result.H_applied:.1f} kN")
    lines.append(f"  Moment M:           {result.M_applied:.1f} kN.m")
    lines.append("")

    # Flexure check
    lines.append("FLEXURE CHECK (P-M Interaction):")
    lines.append(f"  P applied:          {result.P_applied:+.1f} kN")
    lines.append(f"  M applied:          {result.M_applied:.1f} kN.m")
    lines.append(f"  Utilisation:        {result.util_PM:.3f}")
    if result.util_PM <= 1.0:
        lines.append(f"  Status:             OK ✓")
    else:
        lines.append(f"  Status:             FAIL ✗ (exceeds capacity)")
    lines.append("")

    # Shear check
    lines.append("SHEAR CHECK (IRC 112 Cl. 10.3.2):")
    lines.append(f"  VEd (applied):      {result.VEd:.1f} kN")
    lines.append(f"  VRd,c (capacity):   {result.VRdc:.1f} kN")
    lines.append(f"  Utilisation:        {result.shear_util:.3f}")
    if result.shear_util <= 1.0:
        lines.append(f"  Status:             OK ✓")
    else:
        lines.append(f"  Status:             FAIL ✗ (shear reinforcement required)")
    lines.append("")

    # Crack width check
    lines.append("CRACK WIDTH CHECK (SLS, IRC 112 Cl. 12.3.4):")
    lines.append(f"  Calculated width:   {result.crack_width:.3f} mm")
    lines.append(f"  Limit:              {result.crack_width_limit:.3f} mm")
    if result.crack_ok:
        lines.append(f"  Status:             OK ✓")
    else:
        lines.append(f"  Status:             FAIL ✗ (excessive cracking)")
    lines.append("")

    # Curtailment
    lines.append("DETAILING:")
    lines.append(f"  Curtailment depth:  {result.curtailment_depth:.2f} m")
    lines.append("")

    # Overall status
    lines.append("OVERALL STATUS:")
    lines.append(f"  {result.status}")
    lines.append("")
    lines.append("=" * 70)

    return "\n".join(lines)

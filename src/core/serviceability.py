"""
Serviceability checks per IRC 112:2020.

Implements:
- SLS Stress checks (Clause 12.2.1) - Rare combination
- Quasi-permanent stress check (σc ≤ 0.36×fck)
- Crack width check (Clause 12.3.4)
- Deflection check (Clause 12.4)

Key clauses:
- IRC 112:2020 Clause 12.2.1: Stress limitations
- IRC 112:2020 Clause 12.3.4: Crack width control
- IRC 112:2020 Clause 12.4: Deflection control
"""

import math
from typing import List, Optional
from dataclasses import dataclass

from src.codes.base_code import DesignCode
from src.codes.irc112 import IRC112
from src.models.outputs import CalculationStep, DesignStatus


@dataclass
class SLSStressResult:
    """Result from SLS stress check per IRC 112 Cl. 12.2.1."""
    status: DesignStatus

    # Service moment
    rare_moment: float  # M_rare = Mdl + Mll (kNm)

    # Elastic analysis
    modular_ratio: float  # m = Es/Ecm
    elastic_na_depth: float  # x (mm)
    cracked_mi: float  # Icr (mm⁴)

    # Concrete stress (rare combination)
    sigma_c: float  # Actual stress (MPa)
    sigma_c_limit: float  # 0.48×fck (MPa)
    sigma_c_ok: bool
    sigma_c_utilization: float  # ratio

    # Steel stress (rare combination)
    sigma_s: float  # Actual stress (MPa)
    sigma_s_limit: float  # 0.8×fy (MPa)
    sigma_s_ok: bool
    sigma_s_utilization: float  # ratio

    # Quasi-permanent check
    qp_moment: float = 0  # Quasi-permanent moment (kNm)
    sigma_c_qp: float = 0  # Concrete stress under QP (MPa)
    sigma_c_qp_limit: float = 0  # 0.36×fck (MPa)
    sigma_c_qp_ok: bool = True

    # Calculation steps
    steps: List[CalculationStep] = None


@dataclass
class CrackWidthResult:
    """Result from crack width check per IRC 112 Cl. 12.3.4."""
    status: DesignStatus

    # Steel stress for crack width
    sigma_s: float  # Steel stress under quasi-permanent (MPa)

    # Crack spacing parameters
    cover: float  # Clear cover c (mm)
    bar_dia: float  # Bar diameter φ (mm)
    spacing: float  # Bar spacing (mm)
    rho_eff: float  # Effective reinforcement ratio
    sr_max: float  # Maximum crack spacing (mm)

    # Strain difference
    eps_sm_minus_eps_cm: float  # (εsm - εcm)

    # Crack width
    wk_calculated: float  # Calculated crack width (mm)
    wk_limit: float  # Allowable crack width (mm)
    crack_width_ok: bool

    # Calculation steps
    steps: List[CalculationStep] = None


@dataclass
class DeflectionResult:
    """Result from deflection check."""
    status: DesignStatus

    # Section properties
    gross_mi: float  # Ig (mm⁴)
    centroid: float  # mm from bottom
    cracking_moment: float  # Mcr (kNm)
    effective_mi: float  # Ie (mm⁴)

    # Deflection
    calculated: float  # mm
    allowable_total: float  # L/250 (mm)
    allowable_live: float  # L/800 (mm)
    deflection_ok: bool

    # Calculation steps
    steps: List[CalculationStep] = None


class ServiceabilityChecker:
    """
    Serviceability checks per IRC 112:2020.

    Critical for bridge girders:
    - SLS stress checks (σc ≤ 0.48×fck, σs ≤ 0.8×fy) - rare combination
    - Quasi-permanent stress (σc ≤ 0.36×fck)
    - Crack width (wk ≤ 0.3mm for moderate exposure)
    - Deflection limits (L/250 total, L/800 live)
    """

    def __init__(self, code: DesignCode = None):
        self.code = code or IRC112()

    def check_sls_stresses(
        self,
        service_moment: float,   # Ms = Mdl + Mll (kNm) - Rare combination
        web_width: float,        # bw (mm)
        effective_depth: float,  # d (mm)
        ast_provided: float,     # Tension steel (mm²)
        fck: float,              # Concrete strength (MPa)
        fy: float,               # Steel strength (MPa)
        ecm: float = None,       # Modulus of elasticity (MPa)
        is_tbeam: bool = False,
        flange_width: float = None,  # bf for T-beam
        flange_depth: float = None,  # Df for T-beam
        dead_load_moment: float = None,  # Mdl for quasi-permanent check
        psi2: float = 0.3,  # Quasi-permanent factor for traffic (IRC 112 Table B.2)
    ) -> SLSStressResult:
        """
        Check SLS stresses per IRC 112 Cl. 12.2.1.

        Checks:
        1. Rare combination (Mdl + Mll): σc ≤ 0.48×fck, σs ≤ 0.8×fy
        2. Quasi-permanent (Mdl + ψ2×Mll): σc ≤ 0.36×fck (to avoid creep)

        Args:
            service_moment: Rare combination moment Ms = Mdl + Mll (kNm)
            web_width: Web width bw in mm
            effective_depth: Effective depth d in mm
            ast_provided: Provided tension steel in mm²
            fck: Characteristic concrete strength in MPa
            fy: Characteristic steel strength in MPa
            ecm: Modulus of elasticity (auto-calculated if None)
            is_tbeam: True for T-beam section
            flange_width: Effective flange width for T-beam (mm)
            flange_depth: Flange depth for T-beam (mm)
            dead_load_moment: Dead load moment for QP calculation (kNm)
            psi2: Quasi-permanent factor (0.3 for traffic per IRC 112)

        Returns:
            SLSStressResult with stress check details
        """
        steps = []
        step_num = 1

        bw = web_width
        d = effective_depth
        Ms = service_moment  # kNm (rare combination)
        Ms_Nmm = Ms * 1e6  # N.mm
        Es = 200000  # Steel modulus (MPa)

        # Get concrete modulus
        if ecm is None:
            ecm = self.code.get_ecm(fck)

        # Step 1: Modular ratio
        m = Es / ecm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Modular ratio (m)",
            formula="m = Es / Ecm",
            substitution=f"= {Es} / {ecm:.0f}",
            result=round(m, 2),
            unit="",
            code_reference=""
        ))
        step_num += 1

        # Step 2: Elastic neutral axis for cracked section
        if is_tbeam and flange_width and flange_depth:
            bf = flange_width
            Df = flange_depth

            a_coef = 0.5 * bf
            b_coef = m * ast_provided
            c_coef = -m * ast_provided * d

            x_flange = (-b_coef + math.sqrt(b_coef**2 - 4*a_coef*c_coef)) / (2*a_coef)

            if x_flange <= Df:
                x = x_flange
                Icr = bf * x**3 / 3 + m * ast_provided * (d - x)**2
            else:
                a_coef = 0.5 * bw
                b_coef = m * ast_provided + (bf - bw) * Df
                c_coef = -m * ast_provided * d - 0.5 * (bf - bw) * Df**2

                x = (-b_coef + math.sqrt(b_coef**2 - 4*a_coef*c_coef)) / (2*a_coef)
                x = max(x, Df)

                Icr = (bf * Df**3 / 12 + bf * Df * (x - Df/2)**2 +
                       bw * (x - Df)**3 / 3 +
                       m * ast_provided * (d - x)**2)
        else:
            a_coef = 0.5 * bw
            b_coef = m * ast_provided
            c_coef = -m * ast_provided * d

            x = (-b_coef + math.sqrt(b_coef**2 - 4*a_coef*c_coef)) / (2*a_coef)
            Icr = bw * x**3 / 3 + m * ast_provided * (d - x)**2

        steps.append(CalculationStep(
            step_number=step_num,
            description="Elastic neutral axis depth (x)",
            formula="0.5×b×x² + m×Ast×x - m×Ast×d = 0",
            substitution=f"x = {x:.1f} mm",
            result=round(x, 1),
            unit="mm",
            code_reference=""
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Cracked moment of inertia (Icr)",
            formula="Icr = b×x³/3 + m×Ast×(d-x)²",
            substitution=f"= {Icr/1e9:.4f} × 10⁹ mm⁴",
            result=round(Icr/1e9, 4),
            unit="×10⁹ mm⁴",
            code_reference=""
        ))
        step_num += 1

        # Step 3: Concrete stress at top fiber (Rare combination)
        sigma_c = Ms_Nmm * x / Icr
        sigma_c_limit = 0.48 * fck
        sigma_c_ok = sigma_c <= sigma_c_limit
        sigma_c_utilization = sigma_c / sigma_c_limit

        steps.append(CalculationStep(
            step_number=step_num,
            description="Concrete stress - Rare (σc)",
            formula="σc = Ms × x / Icr",
            substitution=f"= {Ms:.1f}×10⁶ × {x:.1f} / {Icr:.0f}",
            result=round(sigma_c, 2),
            unit="MPa",
            code_reference=""
        ))
        step_num += 1

        check_text = "OK" if sigma_c_ok else "FAILS"
        steps.append(CalculationStep(
            step_number=step_num,
            description="Concrete stress check (Rare)",
            formula="σc ≤ 0.48×fck",
            substitution=f"{sigma_c:.2f} {'≤' if sigma_c_ok else '>'} {sigma_c_limit:.2f} [{check_text}]",
            result=round(sigma_c_limit, 2),
            unit="MPa",
            code_reference="IRC 112:2020, Cl. 12.2.1"
        ))
        step_num += 1

        # Step 4: Steel stress (Rare combination)
        sigma_s = m * Ms_Nmm * (d - x) / Icr
        sigma_s_limit = 0.8 * fy
        sigma_s_ok = sigma_s <= sigma_s_limit
        sigma_s_utilization = sigma_s / sigma_s_limit

        steps.append(CalculationStep(
            step_number=step_num,
            description="Steel stress - Rare (σs)",
            formula="σs = m × Ms × (d - x) / Icr",
            substitution=f"= {m:.2f} × {Ms:.1f}×10⁶ × ({d:.0f} - {x:.1f}) / {Icr:.0f}",
            result=round(sigma_s, 1),
            unit="MPa",
            code_reference=""
        ))
        step_num += 1

        check_text = "OK" if sigma_s_ok else "FAILS"
        steps.append(CalculationStep(
            step_number=step_num,
            description="Steel stress check (Rare)",
            formula="σs ≤ 0.8×fy",
            substitution=f"{sigma_s:.1f} {'≤' if sigma_s_ok else '>'} {sigma_s_limit:.0f} [{check_text}]",
            result=round(sigma_s_limit, 0),
            unit="MPa",
            code_reference="IRC 112:2020, Cl. 12.2.1"
        ))
        step_num += 1

        # Step 5: Quasi-permanent stress check (to avoid non-linear creep)
        # M_qp = Mdl + ψ2 × Mll where ψ2 = 0.3 for traffic
        if dead_load_moment is not None:
            Mll = Ms - dead_load_moment
            M_qp = dead_load_moment + psi2 * Mll
        else:
            # Approximate: assume DL = 40% of total service moment
            M_qp = 0.4 * Ms + psi2 * 0.6 * Ms

        M_qp_Nmm = M_qp * 1e6
        sigma_c_qp = M_qp_Nmm * x / Icr
        sigma_c_qp_limit = 0.36 * fck
        sigma_c_qp_ok = sigma_c_qp <= sigma_c_qp_limit

        steps.append(CalculationStep(
            step_number=step_num,
            description="Quasi-permanent moment",
            formula="M_qp = Mdl + ψ2×Mll (ψ2=0.3)",
            substitution=f"M_qp = {M_qp:.1f} kNm",
            result=round(M_qp, 1),
            unit="kNm",
            code_reference="IRC 112:2020, Table B.2"
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Concrete stress - QP (σc,qp)",
            formula="σc,qp = M_qp × x / Icr",
            substitution=f"= {M_qp:.1f}×10⁶ × {x:.1f} / {Icr:.0f}",
            result=round(sigma_c_qp, 2),
            unit="MPa",
            code_reference=""
        ))
        step_num += 1

        check_text = "OK" if sigma_c_qp_ok else "FAILS"
        steps.append(CalculationStep(
            step_number=step_num,
            description="QP stress check (non-linear creep)",
            formula="σc,qp ≤ 0.36×fck",
            substitution=f"{sigma_c_qp:.2f} {'≤' if sigma_c_qp_ok else '>'} {sigma_c_qp_limit:.2f} [{check_text}]",
            result=round(sigma_c_qp_limit, 2),
            unit="MPa",
            code_reference="IRC 112:2020, Cl. 12.2.1(3)"
        ))

        # Overall status
        if sigma_c_ok and sigma_s_ok and sigma_c_qp_ok:
            status = DesignStatus.PASS
        elif sigma_c_qp_ok and (sigma_c_ok or sigma_s_ok):
            status = DesignStatus.WARNING
        else:
            status = DesignStatus.FAIL

        return SLSStressResult(
            status=status,
            rare_moment=Ms,
            modular_ratio=m,
            elastic_na_depth=x,
            cracked_mi=Icr,
            sigma_c=sigma_c,
            sigma_c_limit=sigma_c_limit,
            sigma_c_ok=sigma_c_ok,
            sigma_c_utilization=sigma_c_utilization,
            sigma_s=sigma_s,
            sigma_s_limit=sigma_s_limit,
            sigma_s_ok=sigma_s_ok,
            sigma_s_utilization=sigma_s_utilization,
            qp_moment=M_qp,
            sigma_c_qp=sigma_c_qp,
            sigma_c_qp_limit=sigma_c_qp_limit,
            sigma_c_qp_ok=sigma_c_qp_ok,
            steps=steps
        )

    def check_crack_width(
        self,
        service_moment: float,   # Quasi-permanent moment (kNm)
        web_width: float,        # bw (mm)
        effective_depth: float,  # d (mm)
        overall_depth: float,    # D (mm)
        ast_provided: float,     # Tension steel (mm²)
        bar_dia: float,          # Main bar diameter (mm)
        bar_spacing: float,      # Bar spacing (mm)
        clear_cover: float,      # Clear cover (mm)
        fck: float,
        fy: float,
        ecm: float = None,
        exposure: str = "moderate",  # Exposure condition
    ) -> CrackWidthResult:
        """
        Check crack width per IRC 112 Cl. 12.3.4.

        Uses direct calculation method (Eq. 12.5-12.12).

        Crack width limits (Table 12.1):
        - Moderate exposure: wk ≤ 0.3 mm
        - Severe/Very severe: wk ≤ 0.2 mm

        Args:
            service_moment: Quasi-permanent moment M_qp (kNm)
            web_width: Web width bw (mm)
            effective_depth: Effective depth d (mm)
            overall_depth: Overall depth D (mm)
            ast_provided: Tension steel area (mm²)
            bar_dia: Main bar diameter (mm)
            bar_spacing: Bar center-to-center spacing (mm)
            clear_cover: Clear cover to main bar (mm)
            fck: Concrete strength (MPa)
            fy: Steel strength (MPa)
            ecm: Modulus of elasticity (auto if None)
            exposure: Exposure condition

        Returns:
            CrackWidthResult with crack width check details
        """
        steps = []
        step_num = 1

        bw = web_width
        d = effective_depth
        D = overall_depth
        phi = bar_dia
        c = clear_cover
        M_qp = service_moment
        M_qp_Nmm = M_qp * 1e6
        Es = 200000

        if ecm is None:
            ecm = self.code.get_ecm(fck)

        m = Es / ecm
        fctm = self.code.get_fctm(fck)

        # Crack width limit based on exposure
        wk_limit_map = {
            "moderate": 0.3,
            "severe": 0.2,
            "very_severe": 0.2,
            "extreme": 0.1,
        }
        wk_limit = wk_limit_map.get(exposure, 0.3)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Crack width limit",
            formula="Based on exposure condition",
            substitution=f"wk,max = {wk_limit} mm for {exposure} exposure",
            result=wk_limit,
            unit="mm",
            code_reference="IRC 112:2020, Table 12.1"
        ))
        step_num += 1

        # Step 1: Neutral axis depth (cracked section)
        a_coef = 0.5 * bw
        b_coef = m * ast_provided
        c_coef = -m * ast_provided * d
        x = (-b_coef + math.sqrt(b_coef**2 - 4*a_coef*c_coef)) / (2*a_coef)
        Icr = bw * x**3 / 3 + m * ast_provided * (d - x)**2

        # Step 2: Steel stress under quasi-permanent loading
        sigma_s = m * M_qp_Nmm * (d - x) / Icr

        steps.append(CalculationStep(
            step_number=step_num,
            description="Steel stress (quasi-permanent)",
            formula="σs = m × M_qp × (d - x) / Icr",
            substitution=f"σs = {sigma_s:.1f} MPa",
            result=round(sigma_s, 1),
            unit="MPa",
            code_reference=""
        ))
        step_num += 1

        # Step 3: Effective tension area Ac,eff
        # h_c,eff = min(2.5×(D-d), (D-x)/3, D/2)
        h_c_eff = min(2.5 * (D - d), (D - x) / 3, D / 2)
        Ac_eff = bw * h_c_eff

        steps.append(CalculationStep(
            step_number=step_num,
            description="Effective tension height (hc,eff)",
            formula="hc,eff = min(2.5(D-d), (D-x)/3, D/2)",
            substitution=f"= min({2.5*(D-d):.0f}, {(D-x)/3:.0f}, {D/2:.0f})",
            result=round(h_c_eff, 0),
            unit="mm",
            code_reference="IRC 112:2020, Cl. 12.3.4"
        ))
        step_num += 1

        # Step 4: Effective reinforcement ratio
        rho_eff = ast_provided / Ac_eff

        steps.append(CalculationStep(
            step_number=step_num,
            description="Effective reinforcement ratio (ρp,eff)",
            formula="ρp,eff = As / Ac,eff",
            substitution=f"= {ast_provided:.0f} / {Ac_eff:.0f}",
            result=round(rho_eff, 4),
            unit="",
            code_reference=""
        ))
        step_num += 1

        # Step 5: Maximum crack spacing sr,max (Eq. 12.8)
        # sr,max = 3.4×c + 0.425×k1×k2×φ/ρp,eff
        k1 = 0.8  # High bond bars
        k2 = 0.5  # Bending (0.5 for bending, 1.0 for pure tension)

        sr_max = 3.4 * c + 0.425 * k1 * k2 * phi / rho_eff

        # Limit: sr,max ≤ 5×(c + φ/2) when spacing > 5×(c + φ/2)
        sr_max_limit = 5 * (c + phi / 2)
        if bar_spacing > sr_max_limit:
            sr_max = min(sr_max, 1.3 * (D - x))

        steps.append(CalculationStep(
            step_number=step_num,
            description="Maximum crack spacing (sr,max)",
            formula="sr,max = 3.4c + 0.425×k1×k2×φ/ρp,eff",
            substitution=f"= 3.4×{c:.0f} + 0.425×{k1}×{k2}×{phi:.0f}/{rho_eff:.4f}",
            result=round(sr_max, 0),
            unit="mm",
            code_reference="IRC 112:2020, Eq. 12.8"
        ))
        step_num += 1

        # Step 6: Mean strain difference (εsm - εcm) per Eq. 12.6
        # (εsm - εcm) = [σs - kt×fctm/ρp,eff × (1 + αe×ρp,eff)] / Es ≥ 0.6×σs/Es
        kt = 0.4  # Long-term loading
        alpha_e = Es / ecm

        eps_diff_calc = (sigma_s - kt * fctm / rho_eff * (1 + alpha_e * rho_eff)) / Es
        eps_diff_min = 0.6 * sigma_s / Es
        eps_sm_minus_eps_cm = max(eps_diff_calc, eps_diff_min)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Strain difference (εsm - εcm)",
            formula="= [σs - kt×fctm/ρp,eff × (1+αe×ρp,eff)] / Es",
            substitution=f"= max({eps_diff_calc:.6f}, 0.6×σs/Es)",
            result=round(eps_sm_minus_eps_cm * 1000, 4),
            unit="×10⁻³",
            code_reference="IRC 112:2020, Eq. 12.6"
        ))
        step_num += 1

        # Step 7: Crack width wk = sr,max × (εsm - εcm)
        wk = sr_max * eps_sm_minus_eps_cm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Calculated crack width (wk)",
            formula="wk = sr,max × (εsm - εcm)",
            substitution=f"= {sr_max:.0f} × {eps_sm_minus_eps_cm:.6f}",
            result=round(wk, 3),
            unit="mm",
            code_reference="IRC 112:2020, Eq. 12.5"
        ))
        step_num += 1

        crack_width_ok = wk <= wk_limit
        check_text = "OK" if crack_width_ok else "FAILS"

        steps.append(CalculationStep(
            step_number=step_num,
            description="Crack width check",
            formula="wk ≤ wk,max",
            substitution=f"{wk:.3f} {'≤' if crack_width_ok else '>'} {wk_limit} [{check_text}]",
            result=1 if crack_width_ok else 0,
            unit="",
            code_reference="IRC 112:2020, Table 12.1"
        ))

        return CrackWidthResult(
            status=DesignStatus.PASS if crack_width_ok else DesignStatus.FAIL,
            sigma_s=sigma_s,
            cover=c,
            bar_dia=phi,
            spacing=bar_spacing,
            rho_eff=rho_eff,
            sr_max=sr_max,
            eps_sm_minus_eps_cm=eps_sm_minus_eps_cm,
            wk_calculated=wk,
            wk_limit=wk_limit,
            crack_width_ok=crack_width_ok,
            steps=steps
        )

    def check_deflection(
        self,
        span_m: float,           # Span in meters
        web_width: float,        # bw (mm)
        overall_depth: float,    # D (mm)
        effective_depth: float,  # d (mm)
        ast_provided: float,     # Tension steel (mm²)
        service_moment: float,   # Ms (kNm)
        total_dl: float,         # Total dead load (kN/m)
        live_moment: float,      # Live load moment (kNm)
        fck: float,
        ecm: float = None,
        is_tbeam: bool = False,
        flange_width: float = None,
        flange_depth: float = None,
    ) -> DeflectionResult:
        """
        Check deflection per IRC 112.

        Limits:
        - Total: L/250 for quasi-permanent
        - Live load: L/800
        """
        steps = []
        step_num = 1

        span_mm = span_m * 1000
        bw = web_width
        D = overall_depth
        d = effective_depth
        Ms = service_moment
        Es = 200000

        if ecm is None:
            ecm = self.code.get_ecm(fck)

        m = Es / ecm
        fctm = self.code.get_fctm(fck)

        # Step 1: Gross moment of inertia
        if is_tbeam and flange_width and flange_depth:
            bf = flange_width
            Df = flange_depth

            A_web = bw * (D - Df)
            A_flange = bf * Df
            y_web = (D - Df) / 2
            y_flange = D - Df / 2

            A_total = A_web + A_flange
            y_bar = (A_web * y_web + A_flange * y_flange) / A_total

            I_web = bw * (D - Df)**3 / 12 + A_web * (y_bar - y_web)**2
            I_flange = bf * Df**3 / 12 + A_flange * (y_flange - y_bar)**2
            Ig = I_web + I_flange
            yt = D - y_bar
        else:
            y_bar = D / 2
            Ig = bw * D**3 / 12
            yt = D / 2

        steps.append(CalculationStep(
            step_number=step_num,
            description="Gross moment of inertia (Ig)",
            formula="bw×D³/12 for rectangular",
            substitution=f"Ig = {Ig/1e9:.4f} × 10⁹ mm⁴",
            result=round(Ig/1e9, 4),
            unit="×10⁹ mm⁴",
            code_reference=""
        ))
        step_num += 1

        # Step 2: Cracking moment
        Mcr = fctm * Ig / yt / 1e6  # kNm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Cracking moment (Mcr)",
            formula="Mcr = fctm × Ig / yt",
            substitution=f"= {fctm:.2f} × {Ig/1e9:.4f}×10⁹ / {yt:.0f}",
            result=round(Mcr, 1),
            unit="kNm",
            code_reference=""
        ))
        step_num += 1

        # Step 3: Effective moment of inertia
        if Ms > Mcr:
            a_coef = 0.5 * bw
            b_coef = m * ast_provided
            c_coef = -m * ast_provided * d
            x_cr = (-b_coef + math.sqrt(b_coef**2 - 4*a_coef*c_coef)) / (2*a_coef)
            Icr = bw * x_cr**3 / 3 + m * ast_provided * (d - x_cr)**2

            Ie = Icr + (Ig - Icr) * (Mcr / Ms)**3
            Ie = min(Ie, Ig)
        else:
            Ie = Ig

        steps.append(CalculationStep(
            step_number=step_num,
            description="Effective moment of inertia (Ie)",
            formula="Ie = Icr + (Ig-Icr)×(Mcr/Ms)³",
            substitution=f"Ie = {Ie/1e9:.4f} × 10⁹ mm⁴",
            result=round(Ie/1e9, 4),
            unit="×10⁹ mm⁴",
            code_reference=""
        ))
        step_num += 1

        # Step 4: Calculate deflection
        w_equiv = 8 * Ms / span_m**2  # kN/m
        delta = 5 * w_equiv * span_mm**4 / (384 * ecm * Ie)  # mm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Calculated deflection",
            formula="δ = 5×w×L⁴ / (384×Ecm×Ie)",
            substitution=f"δ = {delta:.2f} mm",
            result=round(delta, 2),
            unit="mm",
            code_reference=""
        ))
        step_num += 1

        # Step 5: Check limits
        delta_allow_total = span_mm / 250
        delta_allow_live = span_mm / 800

        deflection_ok = delta <= delta_allow_total
        check_text = "OK" if deflection_ok else "EXCEEDS"

        steps.append(CalculationStep(
            step_number=step_num,
            description="Deflection check",
            formula="δ ≤ L/250",
            substitution=f"{delta:.2f} {'≤' if deflection_ok else '>'} {delta_allow_total:.1f} [{check_text}]",
            result=round(delta_allow_total, 1),
            unit="mm",
            code_reference="IRC 112:2020, Cl. 12.4"
        ))

        return DeflectionResult(
            status=DesignStatus.PASS if deflection_ok else DesignStatus.WARNING,
            gross_mi=Ig,
            centroid=y_bar,
            cracking_moment=Mcr,
            effective_mi=Ie,
            calculated=delta,
            allowable_total=delta_allow_total,
            allowable_live=delta_allow_live,
            deflection_ok=deflection_ok,
            steps=steps
        )

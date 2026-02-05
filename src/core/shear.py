"""
Shear design calculations per IRC 112:2020.

Implements limit state design method for shear reinforcement using:
- Variable strut inclination method
- Size effect factor (k)
- VRd,c and VRd,max calculations

Key clauses:
- IRC 112:2020 Clause 10.3: ULS - Shear
"""

import math
from typing import List
from dataclasses import dataclass

from src.codes.base_code import DesignCode
from src.codes.irc112 import IRC112
from src.models.outputs import CalculationStep, DesignStatus
from src.utils.constants import BAR_AREAS, STIRRUP_BAR_SIZES


@dataclass
class ShearResult:
    """Internal result from shear design calculations."""
    status: DesignStatus

    # Design forces
    design_shear: float  # Vu at support (kN)
    shear_at_d: float    # Vu at d from support (kN)

    # IRC 112 specific
    size_effect_k: float   # k = 1 + sqrt(200/d)
    rho_l: float          # Reinforcement ratio
    VRdc: float           # Concrete shear resistance (kN)
    VRdc_min: float       # Minimum concrete shear (kN)
    VRd_max: float        # Maximum shear capacity (kN)

    # Reinforcement
    shear_reinf_required: bool
    stirrup_dia: int      # mm
    stirrup_legs: int
    Asv: float           # Stirrup area (mm²)
    spacing_required: float   # mm
    spacing_provided: float   # mm
    spacing_max: float        # mm

    # Checks
    shear_adequacy: bool  # Vu <= VRd,max
    min_reinf_ok: bool

    # Calculation steps
    steps: List[CalculationStep]


class ShearDesigner:
    """
    Shear reinforcement design per IRC 112:2020.

    Uses the variable strut inclination method with size effect.
    """

    def __init__(self, code: DesignCode = None):
        self.code = code or IRC112()

    def design(
        self,
        shear_force: float,      # Vu at support (kN)
        span_m: float,           # Span in meters
        web_width: float,        # Web width bw (mm)
        effective_depth: float,  # Effective depth d (mm)
        fck: float,              # Concrete strength (MPa)
        fy: float,               # Steel yield strength (MPa)
        ast_provided: float,     # Provided tension steel area (mm²)
        factored_udl: float,     # Factored UDL for dead load (kN/m) - for shear at d
        stirrup_dia: int = 12,   # Preferred stirrup diameter (mm)
    ) -> ShearResult:
        """
        Design shear reinforcement per IRC 112:2020.

        Args:
            shear_force: Design shear force at support (kN)
            span_m: Span in meters
            web_width: Web width bw in mm
            effective_depth: Effective depth d in mm
            fck: Characteristic concrete strength in MPa
            fy: Characteristic steel yield strength in MPa
            ast_provided: Provided tension steel area in mm²
            factored_udl: Factored UDL in kN/m (for calculating shear variation)
            stirrup_dia: Preferred stirrup diameter in mm

        Returns:
            ShearResult with complete design details
        """
        steps = []
        step_num = 1

        bw = web_width
        d = effective_depth
        Vu = shear_force  # kN

        # Step 1: Shear at d from support
        d_m = d / 1000  # Convert to meters
        Vu_d = Vu - factored_udl * d_m

        steps.append(CalculationStep(
            step_number=step_num,
            description="Shear at d from support face",
            formula="Vu,d = Vu - wu,DL × d",
            substitution=f"= {Vu:.2f} - {factored_udl:.2f} × {d_m:.3f}",
            result=round(Vu_d, 2),
            unit="kN",
            code_reference="IRC 112:2020"
        ))
        step_num += 1

        # Step 2: Size effect factor k
        # k = 1 + sqrt(200/d) ≤ 2.0
        k = min(1 + math.sqrt(200 / d), 2.0)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Size effect factor (k)",
            formula="k = 1 + √(200/d) ≤ 2.0",
            substitution=f"= 1 + √(200/{d:.0f}) = {1 + math.sqrt(200/d):.3f}",
            result=round(k, 3),
            unit="",
            code_reference="IRC 112:2020, Cl. 10.3.2"
        ))
        step_num += 1

        # Step 3: Reinforcement ratio
        # ρl = Ast / (bw × d) ≤ 0.02
        rho_l = min(ast_provided / (bw * d), 0.02)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Reinforcement ratio (ρl)",
            formula="ρl = Ast / (bw × d) ≤ 0.02",
            substitution=f"= {ast_provided:.0f} / ({bw:.0f} × {d:.0f})",
            result=round(rho_l, 5),
            unit="",
            code_reference=""
        ))
        step_num += 1

        # Step 4: Concrete shear resistance VRd,c
        # VRd,c = 0.12 × k × (80 × ρl × fck)^(1/3) × bw × d
        VRdc_calc = 0.12 * k * (80 * rho_l * fck) ** (1/3) * bw * d / 1000  # kN

        # Minimum value
        # VRd,c,min = 0.031 × k^(3/2) × fck^0.5 × bw × d
        VRdc_min = 0.031 * k ** (3/2) * fck ** 0.5 * bw * d / 1000  # kN

        VRdc = max(VRdc_calc, VRdc_min)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Concrete shear resistance (VRd,c)",
            formula="VRd,c = 0.12 × k × (80×ρl×fck)^(1/3) × bw × d",
            substitution=f"= 0.12 × {k:.3f} × (80×{rho_l:.5f}×{fck})^(1/3) × {bw:.0f} × {d:.0f}",
            result=round(VRdc, 2),
            unit="kN",
            code_reference="IRC 112:2020, Cl. 10.3.2"
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Minimum concrete shear (VRd,c,min)",
            formula="VRd,c,min = 0.031 × k^(3/2) × √fck × bw × d",
            substitution=f"= 0.031 × {k:.3f}^1.5 × √{fck} × {bw:.0f} × {d:.0f}",
            result=round(VRdc_min, 2),
            unit="kN",
            code_reference=""
        ))
        step_num += 1

        # Step 5: Maximum shear resistance VRd,max
        # VRd,max = αcw × bw × z × ν1 × fcd / (cotθ + tanθ)
        # For θ = 45°: cotθ + tanθ = 2
        alpha_cw = 1.0  # For non-prestressed members
        nu1 = 0.6 * (1 - fck / 250)  # Strength reduction factor
        fcd = 0.67 * fck / 1.5  # Design compressive strength
        z = 0.9 * d  # Lever arm

        VRd_max = alpha_cw * bw * z * nu1 * fcd / 2 / 1000  # kN (for θ = 45°)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Maximum shear resistance (VRd,max)",
            formula="VRd,max = αcw × bw × 0.9d × ν1 × fcd / 2",
            substitution=f"= {alpha_cw} × {bw:.0f} × {z:.0f} × {nu1:.3f} × {fcd:.2f} / 2",
            result=round(VRd_max, 2),
            unit="kN",
            code_reference="IRC 112:2020, Cl. 10.3.3"
        ))
        step_num += 1

        # Step 6: Check shear adequacy
        shear_adequacy = Vu_d <= VRd_max

        if not shear_adequacy:
            steps.append(CalculationStep(
                step_number=step_num,
                description="Shear adequacy check",
                formula="Vu,d ≤ VRd,max",
                substitution=f"{Vu_d:.2f} > {VRd_max:.2f} ✗ (SECTION INADEQUATE)",
                result=0,
                unit="",
                code_reference="Increase web width"
            ))

            return ShearResult(
                status=DesignStatus.FAIL,
                design_shear=Vu,
                shear_at_d=Vu_d,
                size_effect_k=k,
                rho_l=rho_l,
                VRdc=VRdc,
                VRdc_min=VRdc_min,
                VRd_max=VRd_max,
                shear_reinf_required=True,
                stirrup_dia=0,
                stirrup_legs=0,
                Asv=0,
                spacing_required=0,
                spacing_provided=0,
                spacing_max=0,
                shear_adequacy=False,
                min_reinf_ok=False,
                steps=steps
            )

        steps.append(CalculationStep(
            step_number=step_num,
            description="Shear adequacy check",
            formula="Vu,d ≤ VRd,max",
            substitution=f"{Vu_d:.2f} ≤ {VRd_max:.2f} ✓",
            result=VRd_max,
            unit="kN",
            code_reference=""
        ))
        step_num += 1

        # Step 7: Check if shear reinforcement required
        shear_reinf_required = Vu_d > VRdc

        if Vu_d <= VRdc:
            steps.append(CalculationStep(
                step_number=step_num,
                description="Shear reinforcement check",
                formula="Vu,d ≤ VRd,c",
                substitution=f"{Vu_d:.2f} ≤ {VRdc:.2f} ✓",
                result=VRdc,
                unit="kN",
                code_reference="Minimum shear reinforcement sufficient"
            ))
            step_num += 1

        else:
            steps.append(CalculationStep(
                step_number=step_num,
                description="Shear reinforcement required",
                formula="Vu,d > VRd,c",
                substitution=f"{Vu_d:.2f} > {VRdc:.2f}",
                result=Vu_d - VRdc,
                unit="kN",
                code_reference="Design shear reinforcement"
            ))
            step_num += 1

        # Step 8: Design stirrups
        n_legs = 2
        Asv = n_legs * BAR_AREAS.get(stirrup_dia, BAR_AREAS[12])

        # Minimum shear reinforcement ratio
        # ρw,min = 0.072 × √fck / fy
        rho_w_min = 0.072 * math.sqrt(fck) / fy

        steps.append(CalculationStep(
            step_number=step_num,
            description="Minimum shear reinforcement ratio",
            formula="ρw,min = 0.072 × √fck / fy",
            substitution=f"= 0.072 × √{fck} / {fy}",
            result=round(rho_w_min, 5),
            unit="",
            code_reference="IRC 112:2020, Cl. 16.3.2"
        ))
        step_num += 1

        if shear_reinf_required:
            # Design stirrup spacing for required shear
            # VRd,s = (Asv/s) × z × fywd × cotθ
            # For θ = 45°: cotθ = 1
            fywd = fy / 1.15  # Design yield strength
            cot_theta = 1.0  # For 45° strut

            # Required Asv/s
            Asv_s_req = Vu_d * 1000 / (z * fywd * cot_theta)

            sv_required = Asv / Asv_s_req

            steps.append(CalculationStep(
                step_number=step_num,
                description="Required stirrup spacing",
                formula="Sv = Asv × z × fywd × cotθ / Vu,d",
                substitution=f"= {Asv:.0f} × {z:.0f} × {fywd:.0f} × 1.0 / {Vu_d*1000:.0f}",
                result=round(sv_required, 0),
                unit="mm",
                code_reference=""
            ))
            step_num += 1
        else:
            # Minimum reinforcement spacing
            sv_required = Asv / (rho_w_min * bw)

        # Maximum spacing
        sv_max = min(0.75 * d, 300)

        # Practical spacing
        sv_provided = min(sv_required, sv_max)
        sv_provided = max(sv_provided, 75)  # Minimum practical
        sv_provided = math.floor(sv_provided / 25) * 25  # Round to 25mm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Maximum stirrup spacing",
            formula="Sv,max = min(0.75×d, 300)",
            substitution=f"= min(0.75×{d:.0f}, 300)",
            result=round(sv_max, 0),
            unit="mm",
            code_reference="IRC 112:2020"
        ))
        step_num += 1

        # Check minimum reinforcement
        Asv_s_provided = Asv / sv_provided
        Asv_s_min = rho_w_min * bw
        min_reinf_ok = Asv_s_provided >= Asv_s_min

        steps.append(CalculationStep(
            step_number=step_num,
            description="Stirrup arrangement",
            formula=f"Provide {n_legs}L-{stirrup_dia}φ @ {sv_provided:.0f}mm c/c",
            substitution=f"Asv = {Asv:.0f} mm², Asv/s = {Asv_s_provided:.3f}",
            result=sv_provided,
            unit="mm",
            code_reference=""
        ))

        return ShearResult(
            status=DesignStatus.PASS if (shear_adequacy and min_reinf_ok) else DesignStatus.WARNING,
            design_shear=Vu,
            shear_at_d=Vu_d,
            size_effect_k=k,
            rho_l=rho_l,
            VRdc=VRdc,
            VRdc_min=VRdc_min,
            VRd_max=VRd_max,
            shear_reinf_required=shear_reinf_required,
            stirrup_dia=stirrup_dia,
            stirrup_legs=n_legs,
            Asv=Asv,
            spacing_required=sv_required,
            spacing_provided=sv_provided,
            spacing_max=sv_max,
            shear_adequacy=shear_adequacy,
            min_reinf_ok=min_reinf_ok,
            steps=steps
        )

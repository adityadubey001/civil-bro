"""
Flexural design calculations per IRC 112:2020.

Implements limit state design method for:
- T-beam sections (with deck slab as flange)
- Rectangular beam sections
- Singly and doubly reinforced designs

Key clauses:
- IRC 112:2020 Clause 10.3: ULS - Flexure
"""

import math
from typing import Tuple, List, Optional
from dataclasses import dataclass

from src.codes.base_code import DesignCode
from src.codes.irc112 import IRC112
from src.models.outputs import CalculationStep, DesignStatus
from src.utils.constants import BAR_AREAS, STANDARD_BAR_SIZES


@dataclass
class FlexureResult:
    """Internal result from flexural design calculations."""
    # Status
    status: DesignStatus

    # Design moment
    design_moment: float  # kNm
    moment_capacity_limit: float  # Mu_lim in kNm

    # T-beam specific
    is_tbeam: bool
    flange_width: Optional[float]  # bf in mm
    flange_depth: Optional[float]  # Df in mm
    moment_if_na_at_flange: Optional[float]  # Mf in kNm
    na_in_flange: Optional[bool]

    # Tension reinforcement
    required_ast: float  # mm²
    provided_ast: float  # mm²
    tension_bar_dia: int  # mm
    tension_bar_count: int
    tension_bar_layers: int
    tension_arrangement: str

    # Compression reinforcement (doubly reinforced)
    is_doubly: bool
    required_asc: Optional[float]  # mm²
    provided_asc: Optional[float]  # mm²
    compression_bar_dia: Optional[int]
    compression_bar_count: Optional[int]
    compression_arrangement: Optional[str]

    # Section analysis
    xu: float  # Neutral axis depth (mm)
    xu_max: float  # Maximum neutral axis depth (mm)
    xu_d_ratio: float
    lever_arm: float  # z (mm)
    moment_capacity: float  # kNm

    # Checks
    is_under_reinforced: bool
    min_ast_ok: bool
    max_ast_ok: bool
    pt_provided: float  # %

    # Calculation steps
    steps: List[CalculationStep]


class FlexureDesigner:
    """
    Flexural reinforcement design per IRC 112:2020.

    Supports both T-beam (with deck flange) and rectangular sections.
    Handles singly and doubly reinforced design.
    """

    def __init__(self, code: DesignCode = None):
        self.code = code or IRC112()

    def design(
        self,
        moment: float,           # Design moment Mu (kNm)
        web_width: float,        # Web width bw (mm)
        depth: float,            # Overall depth D (mm)
        effective_depth: float,  # Effective depth d (mm)
        fck: float,              # Concrete strength (MPa)
        fy: float,               # Steel yield strength (MPa)
        is_tbeam: bool = False,
        flange_width: float = None,   # bf for T-beam (mm)
        flange_depth: float = None,   # Df for T-beam (mm)
        effective_cover: float = 50,  # d' for compression steel (mm)
        allow_doubly: bool = True,
        preferred_bars: List[int] = None,
    ) -> FlexureResult:
        """
        Design flexural reinforcement per IRC 112:2020.

        Args:
            moment: Design bending moment in kNm
            web_width: Web width (bw) in mm
            depth: Overall depth (D) in mm
            effective_depth: Effective depth to tension steel (d) in mm
            fck: Characteristic concrete strength in MPa
            fy: Characteristic steel yield strength in MPa
            is_tbeam: True for T-beam, False for rectangular
            flange_width: Effective flange width (bf) for T-beam in mm
            flange_depth: Flange thickness (Df) for T-beam in mm
            effective_cover: Cover to compression steel centroid in mm
            allow_doubly: If True, use doubly reinforced design when needed
            preferred_bars: List of preferred bar diameters in mm

        Returns:
            FlexureResult with complete design details
        """
        if preferred_bars is None:
            preferred_bars = [25, 28, 32, 36]  # Larger bars for bridge girders

        steps = []
        step_num = 1

        bw = web_width
        d = effective_depth
        d_prime = effective_cover
        D = depth

        # Design constants for IRC 112
        # xu,max/d ratio for Fe500 = 0.46
        xu_max_d = 0.46 if fy >= 500 else 0.48
        xu_max = xu_max_d * d

        Mu = moment  # kNm
        Mu_Nmm = Mu * 1e6  # N.mm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Maximum neutral axis depth ratio",
            formula="xu,max/d = 0.46 (for Fe500)",
            substitution=f"xu,max = {xu_max_d} × {d:.0f} = {xu_max:.0f} mm",
            result=round(xu_max, 0),
            unit="mm",
            code_reference="IRC 112:2020"
        ))
        step_num += 1

        if is_tbeam and flange_width and flange_depth:
            return self._design_tbeam(
                Mu_Nmm, bw, d, D, d_prime, fck, fy, xu_max, xu_max_d,
                flange_width, flange_depth, allow_doubly, preferred_bars, steps, step_num
            )
        else:
            return self._design_rectangular(
                Mu_Nmm, bw, d, D, d_prime, fck, fy, xu_max, xu_max_d,
                allow_doubly, preferred_bars, steps, step_num
            )

    def _design_tbeam(
        self,
        Mu: float,  # N.mm
        bw: float,  # Web width (mm)
        d: float,   # Effective depth (mm)
        D: float,   # Overall depth (mm)
        d_prime: float,  # mm
        fck: float,
        fy: float,
        xu_max: float,
        xu_max_d: float,
        bf: float,  # Flange width (mm)
        Df: float,  # Flange depth (mm)
        allow_doubly: bool,
        preferred_bars: List[int],
        steps: List[CalculationStep],
        step_num: int,
    ) -> FlexureResult:
        """Design T-beam section per IRC 112."""

        Mu_kNm = Mu / 1e6

        # Step: Check if NA is in flange or web
        # Moment capacity if NA is at bottom of flange
        Mf = 0.36 * fck * bf * Df * (d - 0.42 * Df) / 1e6  # kNm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Moment capacity if NA at flange bottom (Mf)",
            formula="Mf = 0.36 × fck × bf × Df × (d - 0.42×Df)",
            substitution=f"= 0.36 × {fck} × {bf:.0f} × {Df:.0f} × ({d:.0f} - 0.42×{Df:.0f})",
            result=round(Mf, 2),
            unit="kNm",
            code_reference="IRC 112:2020"
        ))
        step_num += 1

        if Mu_kNm <= Mf:
            # NA is within flange - design as rectangular beam with width bf
            steps.append(CalculationStep(
                step_number=step_num,
                description="NA location check",
                formula="Mu ≤ Mf → NA is within flange",
                substitution=f"{Mu_kNm:.2f} kNm ≤ {Mf:.2f} kNm ✓",
                result=Mu_kNm,
                unit="kNm",
                code_reference="Design as rectangular beam with width = bf"
            ))
            step_num += 1

            # Calculate limiting moment with flange width
            Mu_lim = 0.36 * fck * bf * xu_max * (d - 0.42 * xu_max) / 1e6  # kNm

            steps.append(CalculationStep(
                step_number=step_num,
                description="Limiting moment (Mu,lim)",
                formula="Mu,lim = 0.36 × fck × bf × xu,max × (d - 0.42×xu,max)",
                substitution=f"= 0.36 × {fck} × {bf:.0f} × {xu_max:.0f} × ({d:.0f} - 0.42×{xu_max:.0f})",
                result=round(Mu_lim, 2),
                unit="kNm",
                code_reference="IRC 112:2020"
            ))
            step_num += 1

            if Mu_kNm <= Mu_lim:
                # Singly reinforced T-beam (NA in flange)
                return self._design_singly_tbeam_na_in_flange(
                    Mu, bw, d, D, fck, fy, xu_max, xu_max_d, Mu_lim, Mf,
                    bf, Df, preferred_bars, steps, step_num
                )
            else:
                if allow_doubly:
                    return self._design_doubly_reinforced(
                        Mu, bf, d, D, d_prime, fck, fy, xu_max, xu_max_d, Mu_lim,
                        preferred_bars, steps, step_num, is_tbeam=True, bw=bw, Df=Df, Mf=Mf
                    )
                else:
                    return self._return_fail_result(
                        Mu_kNm, Mu_lim, xu_max, steps, step_num,
                        is_tbeam=True, bf=bf, Df=Df, Mf=Mf
                    )
        else:
            # NA is in web - T-beam analysis required
            steps.append(CalculationStep(
                step_number=step_num,
                description="NA location check",
                formula="Mu > Mf → NA is in web",
                substitution=f"{Mu_kNm:.2f} kNm > {Mf:.2f} kNm",
                result=Mu_kNm,
                unit="kNm",
                code_reference="T-beam analysis with NA in web"
            ))
            step_num += 1

            return self._design_tbeam_na_in_web(
                Mu, bw, d, D, d_prime, fck, fy, xu_max, xu_max_d, Mf,
                bf, Df, allow_doubly, preferred_bars, steps, step_num
            )

    def _design_singly_tbeam_na_in_flange(
        self,
        Mu: float,  # N.mm
        bw: float,
        d: float,
        D: float,
        fck: float,
        fy: float,
        xu_max: float,
        xu_max_d: float,
        Mu_lim: float,  # kNm
        Mf: float,  # kNm
        bf: float,
        Df: float,
        preferred_bars: List[int],
        steps: List[CalculationStep],
        step_num: int,
    ) -> FlexureResult:
        """Design singly reinforced T-beam with NA in flange."""

        Mu_kNm = Mu / 1e6

        steps.append(CalculationStep(
            step_number=step_num,
            description="Singly reinforced design applicable",
            formula="Mu ≤ Mu,lim",
            substitution=f"{Mu_kNm:.2f} kNm ≤ {Mu_lim:.2f} kNm ✓",
            result=Mu_kNm,
            unit="kNm",
            code_reference=""
        ))
        step_num += 1

        # Calculate Ast using quadratic formula
        # Mu = 0.87 × fy × Ast × (d - 0.42×xu)
        # where xu = 0.87×fy×Ast / (0.36×fck×bf)
        # Substituting: 0.42 × 0.87×fy×Ast / (0.36×fck×bf) = xu
        # Mu = 0.87×fy×Ast×d - 0.42×(0.87×fy)²×Ast² / (0.36×fck×bf)

        A = 0.42 * (0.87 * fy) ** 2 / (0.36 * fck * bf)
        B = -0.87 * fy * d
        C = Mu

        discriminant = B**2 - 4*A*C
        Ast_required = (-B - math.sqrt(discriminant)) / (2*A)

        # Calculate xu
        xu = 0.87 * fy * Ast_required / (0.36 * fck * bf)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Neutral axis depth (xu)",
            formula="xu = 0.87×fy×Ast / (0.36×fck×bf)",
            substitution=f"xu = {xu:.0f} mm (< Df = {Df:.0f} mm ✓)",
            result=round(xu, 0),
            unit="mm",
            code_reference=""
        ))
        step_num += 1

        z = d - 0.42 * xu

        steps.append(CalculationStep(
            step_number=step_num,
            description="Required tension steel area (Ast)",
            formula="Ast = Mu / (0.87×fy×z)",
            substitution=f"= {Mu/1e6:.2f}×10⁶ / (0.87×{fy}×{z:.0f})",
            result=round(Ast_required, 0),
            unit="mm²",
            code_reference=""
        ))
        step_num += 1

        # Select bars
        bar_dia, bar_count, bar_layers, Ast_provided = self._select_bars_with_layers(
            Ast_required, bw, preferred_bars
        )
        arrangement = f"{bar_count}-{bar_dia}φ ({bar_layers}L×{bar_count//bar_layers})"

        steps.append(CalculationStep(
            step_number=step_num,
            description="Reinforcement provided",
            formula=f"Provide {arrangement}",
            substitution=f"Ast,provided = {Ast_provided:.0f} mm²",
            result=round(Ast_provided, 0),
            unit="mm²",
            code_reference=""
        ))
        step_num += 1

        # Checks
        fctm = self.code.get_fctm(fck)
        Ast_min = max(0.26 * fctm / fy * bw * d, 0.0013 * bw * d)
        Ast_max = 0.04 * bw * D
        min_ok = Ast_provided >= Ast_min
        max_ok = Ast_provided <= Ast_max
        pt_provided = 100 * Ast_provided / (bw * d)

        # Recalculate with provided steel
        xu_actual = 0.87 * fy * Ast_provided / (0.36 * fck * bf)
        z_actual = d - 0.42 * xu_actual
        Mu_capacity = 0.87 * fy * Ast_provided * z_actual / 1e6

        xu_d_ratio = xu_actual / d
        is_under_reinforced = xu_actual <= xu_max

        return FlexureResult(
            status=DesignStatus.PASS if (min_ok and max_ok and is_under_reinforced) else DesignStatus.WARNING,
            design_moment=Mu_kNm,
            moment_capacity_limit=Mu_lim,
            is_tbeam=True,
            flange_width=bf,
            flange_depth=Df,
            moment_if_na_at_flange=Mf,
            na_in_flange=True,
            required_ast=Ast_required,
            provided_ast=Ast_provided,
            tension_bar_dia=bar_dia,
            tension_bar_count=bar_count,
            tension_bar_layers=bar_layers,
            tension_arrangement=arrangement,
            is_doubly=False,
            required_asc=None,
            provided_asc=None,
            compression_bar_dia=None,
            compression_bar_count=None,
            compression_arrangement=None,
            xu=xu_actual,
            xu_max=xu_max,
            xu_d_ratio=xu_d_ratio,
            lever_arm=z_actual,
            moment_capacity=Mu_capacity,
            is_under_reinforced=is_under_reinforced,
            min_ast_ok=min_ok,
            max_ast_ok=max_ok,
            pt_provided=pt_provided,
            steps=steps
        )

    def _design_tbeam_na_in_web(
        self,
        Mu: float,  # N.mm
        bw: float,
        d: float,
        D: float,
        d_prime: float,
        fck: float,
        fy: float,
        xu_max: float,
        xu_max_d: float,
        Mf: float,
        bf: float,
        Df: float,
        allow_doubly: bool,
        preferred_bars: List[int],
        steps: List[CalculationStep],
        step_num: int,
    ) -> FlexureResult:
        """Design T-beam with NA in web."""

        Mu_kNm = Mu / 1e6

        # Force in flange overhang
        yf = Df / 2  # CG of flange from top
        Cf = 0.447 * fck * (bf - bw) * Df  # N
        Muf = Cf * (d - yf) / 1e6  # kNm (moment from flange)

        steps.append(CalculationStep(
            step_number=step_num,
            description="Compression force in flange overhang (Cf)",
            formula="Cf = 0.447 × fck × (bf - bw) × Df",
            substitution=f"= 0.447 × {fck} × ({bf:.0f} - {bw:.0f}) × {Df:.0f}",
            result=round(Cf/1000, 1),
            unit="kN",
            code_reference=""
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Moment from flange (Muf)",
            formula="Muf = Cf × (d - Df/2)",
            substitution=f"= {Cf/1000:.1f} × ({d:.0f} - {Df/2:.0f})",
            result=round(Muf, 2),
            unit="kNm",
            code_reference=""
        ))
        step_num += 1

        # Moment to be resisted by web
        Muw = Mu_kNm - Muf  # kNm

        steps.append(CalculationStep(
            step_number=step_num,
            description="Moment for web (Muw)",
            formula="Muw = Mu - Muf",
            substitution=f"= {Mu_kNm:.2f} - {Muf:.2f}",
            result=round(Muw, 2),
            unit="kNm",
            code_reference=""
        ))
        step_num += 1

        # Steel for flange
        Ast_flange = Cf / (0.87 * fy)

        # Calculate limiting moment for web
        Mu_lim_web = 0.36 * fck * bw * xu_max * (d - 0.42 * xu_max) / 1e6

        if Muw * 1e6 <= Mu_lim_web * 1e6:
            # Web portion is singly reinforced
            # Solve for Ast_web from Muw
            Muw_Nmm = Muw * 1e6

            A = 0.42 * (0.87 * fy) ** 2 / (0.36 * fck * bw)
            B = -0.87 * fy * d
            C = Muw_Nmm

            discriminant = B**2 - 4*A*C
            if discriminant >= 0:
                Ast_web = (-B - math.sqrt(discriminant)) / (2*A)
            else:
                # Section inadequate
                Ast_web = Muw_Nmm / (0.87 * fy * 0.9 * d)

            Ast_required = Ast_flange + Ast_web

            # Calculate xu in web
            xu_web = 0.87 * fy * Ast_web / (0.36 * fck * bw)
            xu = Df + xu_web  # Total NA depth from top

            steps.append(CalculationStep(
                step_number=step_num,
                description="Steel for flange (Ast,f)",
                formula="Ast,f = Cf / (0.87×fy)",
                substitution=f"= {Cf:.0f} / (0.87×{fy})",
                result=round(Ast_flange, 0),
                unit="mm²",
                code_reference=""
            ))
            step_num += 1

            steps.append(CalculationStep(
                step_number=step_num,
                description="Steel for web (Ast,w)",
                formula="From Muw equilibrium",
                substitution=f"Ast,w = {Ast_web:.0f} mm²",
                result=round(Ast_web, 0),
                unit="mm²",
                code_reference=""
            ))
            step_num += 1

            steps.append(CalculationStep(
                step_number=step_num,
                description="Total tension steel (Ast)",
                formula="Ast = Ast,f + Ast,w",
                substitution=f"= {Ast_flange:.0f} + {Ast_web:.0f}",
                result=round(Ast_required, 0),
                unit="mm²",
                code_reference=""
            ))
            step_num += 1

            # Select bars
            bar_dia, bar_count, bar_layers, Ast_provided = self._select_bars_with_layers(
                Ast_required, bw, preferred_bars
            )
            arrangement = f"{bar_count}-{bar_dia}φ ({bar_layers}L×{bar_count//bar_layers})"

            # Checks
            fctm = self.code.get_fctm(fck)
            Ast_min = max(0.26 * fctm / fy * bw * d, 0.0013 * bw * d)
            Ast_max = 0.04 * bw * D
            min_ok = Ast_provided >= Ast_min
            max_ok = Ast_provided <= Ast_max
            pt_provided = 100 * Ast_provided / (bw * d)

            z = d - 0.42 * xu
            Mu_capacity = (Cf * (d - yf) + 0.87 * fy * Ast_web * (d - Df - 0.42*xu_web)) / 1e6
            xu_d_ratio = xu / d
            is_under_reinforced = xu <= xu_max

            Mu_lim = Muf + Mu_lim_web

            return FlexureResult(
                status=DesignStatus.PASS if (min_ok and max_ok and is_under_reinforced) else DesignStatus.WARNING,
                design_moment=Mu_kNm,
                moment_capacity_limit=Mu_lim,
                is_tbeam=True,
                flange_width=bf,
                flange_depth=Df,
                moment_if_na_at_flange=Mf,
                na_in_flange=False,
                required_ast=Ast_required,
                provided_ast=Ast_provided,
                tension_bar_dia=bar_dia,
                tension_bar_count=bar_count,
                tension_bar_layers=bar_layers,
                tension_arrangement=arrangement,
                is_doubly=False,
                required_asc=None,
                provided_asc=None,
                compression_bar_dia=None,
                compression_bar_count=None,
                compression_arrangement=None,
                xu=xu,
                xu_max=xu_max,
                xu_d_ratio=xu_d_ratio,
                lever_arm=z,
                moment_capacity=Mu_capacity,
                is_under_reinforced=is_under_reinforced,
                min_ast_ok=min_ok,
                max_ast_ok=max_ok,
                pt_provided=pt_provided,
                steps=steps
            )
        else:
            # Web portion needs doubly reinforced design
            if allow_doubly:
                Mu_lim = Muf + Mu_lim_web
                return self._design_doubly_reinforced(
                    Mu, bw, d, D, d_prime, fck, fy, xu_max, xu_max_d, Mu_lim,
                    preferred_bars, steps, step_num, is_tbeam=True, bw=bw, bf=bf, Df=Df, Mf=Mf
                )
            else:
                Mu_lim = Muf + Mu_lim_web
                return self._return_fail_result(
                    Mu_kNm, Mu_lim, xu_max, steps, step_num,
                    is_tbeam=True, bf=bf, Df=Df, Mf=Mf
                )

    def _design_rectangular(
        self,
        Mu: float,  # N.mm
        b: float,   # mm
        d: float,   # mm
        D: float,   # mm
        d_prime: float,  # mm
        fck: float,
        fy: float,
        xu_max: float,
        xu_max_d: float,
        allow_doubly: bool,
        preferred_bars: List[int],
        steps: List[CalculationStep],
        step_num: int,
    ) -> FlexureResult:
        """Design rectangular beam section."""

        Mu_kNm = Mu / 1e6

        # Limiting moment
        Mu_lim = 0.36 * fck * b * xu_max * (d - 0.42 * xu_max) / 1e6

        steps.append(CalculationStep(
            step_number=step_num,
            description="Limiting moment (Mu,lim)",
            formula="Mu,lim = 0.36 × fck × b × xu,max × (d - 0.42×xu,max)",
            substitution=f"= 0.36 × {fck} × {b:.0f} × {xu_max:.0f} × ({d:.0f} - 0.42×{xu_max:.0f})",
            result=round(Mu_lim, 2),
            unit="kNm",
            code_reference="IRC 112:2020"
        ))
        step_num += 1

        if Mu_kNm <= Mu_lim:
            # Singly reinforced
            steps.append(CalculationStep(
                step_number=step_num,
                description="Singly reinforced design applicable",
                formula="Mu ≤ Mu,lim",
                substitution=f"{Mu_kNm:.2f} kNm ≤ {Mu_lim:.2f} kNm ✓",
                result=Mu_kNm,
                unit="kNm",
                code_reference=""
            ))
            step_num += 1

            # Calculate Ast
            A = 0.42 * (0.87 * fy) ** 2 / (0.36 * fck * b)
            B = -0.87 * fy * d
            C = Mu

            discriminant = B**2 - 4*A*C
            Ast_required = (-B - math.sqrt(discriminant)) / (2*A)

            xu = 0.87 * fy * Ast_required / (0.36 * fck * b)
            z = d - 0.42 * xu

            steps.append(CalculationStep(
                step_number=step_num,
                description="Neutral axis depth (xu)",
                formula="xu = 0.87×fy×Ast / (0.36×fck×b)",
                substitution=f"xu = {xu:.0f} mm",
                result=round(xu, 0),
                unit="mm",
                code_reference=""
            ))
            step_num += 1

            steps.append(CalculationStep(
                step_number=step_num,
                description="Required tension steel (Ast)",
                formula="Ast = Mu / (0.87×fy×z)",
                substitution=f"= {Mu_kNm:.2f}×10⁶ / (0.87×{fy}×{z:.0f})",
                result=round(Ast_required, 0),
                unit="mm²",
                code_reference=""
            ))
            step_num += 1

            # Select bars
            bar_dia, bar_count, bar_layers, Ast_provided = self._select_bars_with_layers(
                Ast_required, b, preferred_bars
            )
            arrangement = f"{bar_count}-{bar_dia}φ ({bar_layers}L×{bar_count//bar_layers})"

            # Checks
            fctm = self.code.get_fctm(fck)
            Ast_min = max(0.26 * fctm / fy * b * d, 0.0013 * b * d)
            Ast_max = 0.04 * b * D
            min_ok = Ast_provided >= Ast_min
            max_ok = Ast_provided <= Ast_max
            pt_provided = 100 * Ast_provided / (b * d)

            xu_actual = 0.87 * fy * Ast_provided / (0.36 * fck * b)
            z_actual = d - 0.42 * xu_actual
            Mu_capacity = 0.87 * fy * Ast_provided * z_actual / 1e6
            xu_d_ratio = xu_actual / d
            is_under_reinforced = xu_actual <= xu_max

            return FlexureResult(
                status=DesignStatus.PASS if (min_ok and max_ok and is_under_reinforced) else DesignStatus.WARNING,
                design_moment=Mu_kNm,
                moment_capacity_limit=Mu_lim,
                is_tbeam=False,
                flange_width=None,
                flange_depth=None,
                moment_if_na_at_flange=None,
                na_in_flange=None,
                required_ast=Ast_required,
                provided_ast=Ast_provided,
                tension_bar_dia=bar_dia,
                tension_bar_count=bar_count,
                tension_bar_layers=bar_layers,
                tension_arrangement=arrangement,
                is_doubly=False,
                required_asc=None,
                provided_asc=None,
                compression_bar_dia=None,
                compression_bar_count=None,
                compression_arrangement=None,
                xu=xu_actual,
                xu_max=xu_max,
                xu_d_ratio=xu_d_ratio,
                lever_arm=z_actual,
                moment_capacity=Mu_capacity,
                is_under_reinforced=is_under_reinforced,
                min_ast_ok=min_ok,
                max_ast_ok=max_ok,
                pt_provided=pt_provided,
                steps=steps
            )
        else:
            if allow_doubly:
                return self._design_doubly_reinforced(
                    Mu, b, d, D, d_prime, fck, fy, xu_max, xu_max_d, Mu_lim,
                    preferred_bars, steps, step_num, is_tbeam=False
                )
            else:
                return self._return_fail_result(
                    Mu_kNm, Mu_lim, xu_max, steps, step_num, is_tbeam=False
                )

    def _design_doubly_reinforced(
        self,
        Mu: float,  # N.mm
        b: float,   # Width to use for compression (bf for T-beam NA in flange, bw otherwise)
        d: float,
        D: float,
        d_prime: float,
        fck: float,
        fy: float,
        xu_max: float,
        xu_max_d: float,
        Mu_lim: float,  # kNm
        preferred_bars: List[int],
        steps: List[CalculationStep],
        step_num: int,
        is_tbeam: bool = False,
        bw: float = None,
        bf: float = None,
        Df: float = None,
        Mf: float = None,
    ) -> FlexureResult:
        """Design doubly reinforced section."""

        Mu_kNm = Mu / 1e6
        Mu_lim_Nmm = Mu_lim * 1e6

        steps.append(CalculationStep(
            step_number=step_num,
            description="Doubly reinforced design required",
            formula="Mu > Mu,lim",
            substitution=f"{Mu_kNm:.2f} kNm > {Mu_lim:.2f} kNm",
            result=Mu_kNm,
            unit="kNm",
            code_reference=""
        ))
        step_num += 1

        # Additional moment
        Mu2 = Mu - Mu_lim_Nmm  # N.mm
        Mu2_kNm = Mu2 / 1e6

        steps.append(CalculationStep(
            step_number=step_num,
            description="Additional moment for compression steel",
            formula="Mu2 = Mu - Mu,lim",
            substitution=f"= {Mu_kNm:.2f} - {Mu_lim:.2f}",
            result=round(Mu2_kNm, 2),
            unit="kNm",
            code_reference=""
        ))
        step_num += 1

        # Tension steel for Mu_lim
        Ast1 = (0.36 * fck * b * xu_max) / (0.87 * fy)

        # Compression steel stress
        fyd = fy / 1.15
        Es = 200000
        epsilon_cu = 0.0035
        epsilon_sc = epsilon_cu * (xu_max - d_prime) / xu_max
        epsilon_y = fyd / Es

        if epsilon_sc >= epsilon_y:
            fsc = fyd
        else:
            fsc = Es * epsilon_sc

        # Compression steel area
        fcc = 0.446 * fck
        Asc_required = Mu2 / ((fsc - fcc) * (d - d_prime))

        # Additional tension steel
        Ast2 = Asc_required * (fsc - fcc) / (0.87 * fy)
        Ast_required = Ast1 + Ast2

        steps.append(CalculationStep(
            step_number=step_num,
            description="Tension steel for Mu,lim (Ast1)",
            formula="Ast1 = 0.36×fck×b×xu,max / (0.87×fy)",
            substitution=f"= 0.36×{fck}×{b:.0f}×{xu_max:.0f} / (0.87×{fy})",
            result=round(Ast1, 0),
            unit="mm²",
            code_reference=""
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Compression steel stress (fsc)",
            formula="fsc = fy/γs or Es×εsc",
            substitution=f"fsc = {fsc:.0f} MPa",
            result=round(fsc, 0),
            unit="MPa",
            code_reference=""
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Compression steel (Asc)",
            formula="Asc = Mu2 / [(fsc - fcc) × (d - d')]",
            substitution=f"= {Mu2/1e6:.2f}×10⁶ / [({fsc:.0f} - {fcc:.1f}) × ({d:.0f} - {d_prime:.0f})]",
            result=round(Asc_required, 0),
            unit="mm²",
            code_reference=""
        ))
        step_num += 1

        steps.append(CalculationStep(
            step_number=step_num,
            description="Total tension steel (Ast)",
            formula="Ast = Ast1 + Ast2",
            substitution=f"= {Ast1:.0f} + {Ast2:.0f}",
            result=round(Ast_required, 0),
            unit="mm²",
            code_reference=""
        ))
        step_num += 1

        # Select bars
        width_for_bars = bw if bw else b
        t_bar_dia, t_bar_count, t_bar_layers, Ast_provided = self._select_bars_with_layers(
            Ast_required, width_for_bars, preferred_bars
        )
        t_arrangement = f"{t_bar_count}-{t_bar_dia}φ ({t_bar_layers}L×{t_bar_count//t_bar_layers})"

        c_bar_dia, c_bar_count, _, Asc_provided = self._select_bars_with_layers(
            Asc_required, width_for_bars, [20, 25, 28]
        )
        c_arrangement = f"{c_bar_count}-{c_bar_dia}φ"

        # Checks
        fctm = self.code.get_fctm(fck)
        Ast_min = max(0.26 * fctm / fy * width_for_bars * d, 0.0013 * width_for_bars * d)
        Ast_max = 0.04 * width_for_bars * D
        min_ok = Ast_provided >= Ast_min
        max_ok = Ast_provided <= Ast_max
        pt_provided = 100 * Ast_provided / (width_for_bars * d)

        z = d - 0.42 * xu_max
        Mu_capacity = (0.87 * fy * Ast1 * z + Asc_provided * (fsc - fcc) * (d - d_prime)) / 1e6

        return FlexureResult(
            status=DesignStatus.PASS if (min_ok and max_ok) else DesignStatus.WARNING,
            design_moment=Mu_kNm,
            moment_capacity_limit=Mu_lim,
            is_tbeam=is_tbeam,
            flange_width=bf,
            flange_depth=Df,
            moment_if_na_at_flange=Mf,
            na_in_flange=True if (is_tbeam and Mf and Mu_kNm <= Mf) else False,
            required_ast=Ast_required,
            provided_ast=Ast_provided,
            tension_bar_dia=t_bar_dia,
            tension_bar_count=t_bar_count,
            tension_bar_layers=t_bar_layers,
            tension_arrangement=t_arrangement,
            is_doubly=True,
            required_asc=Asc_required,
            provided_asc=Asc_provided,
            compression_bar_dia=c_bar_dia,
            compression_bar_count=c_bar_count,
            compression_arrangement=c_arrangement,
            xu=xu_max,
            xu_max=xu_max,
            xu_d_ratio=xu_max_d,
            lever_arm=z,
            moment_capacity=Mu_capacity,
            is_under_reinforced=True,
            min_ast_ok=min_ok,
            max_ast_ok=max_ok,
            pt_provided=pt_provided,
            steps=steps
        )

    def _return_fail_result(
        self,
        Mu_kNm: float,
        Mu_lim: float,
        xu_max: float,
        steps: List[CalculationStep],
        step_num: int,
        is_tbeam: bool = False,
        bf: float = None,
        Df: float = None,
        Mf: float = None,
    ) -> FlexureResult:
        """Return a failed design result."""

        steps.append(CalculationStep(
            step_number=step_num,
            description="Section capacity check",
            formula="Mu > Mu,lim (Section undersized)",
            substitution=f"{Mu_kNm:.2f} kNm > {Mu_lim:.2f} kNm",
            result=0,
            unit="",
            code_reference="Increase section size or allow doubly reinforced"
        ))

        return FlexureResult(
            status=DesignStatus.FAIL,
            design_moment=Mu_kNm,
            moment_capacity_limit=Mu_lim,
            is_tbeam=is_tbeam,
            flange_width=bf,
            flange_depth=Df,
            moment_if_na_at_flange=Mf,
            na_in_flange=None,
            required_ast=0,
            provided_ast=0,
            tension_bar_dia=0,
            tension_bar_count=0,
            tension_bar_layers=0,
            tension_arrangement="N/A",
            is_doubly=False,
            required_asc=None,
            provided_asc=None,
            compression_bar_dia=None,
            compression_bar_count=None,
            compression_arrangement=None,
            xu=0,
            xu_max=xu_max,
            xu_d_ratio=0,
            lever_arm=0,
            moment_capacity=Mu_lim,
            is_under_reinforced=False,
            min_ast_ok=False,
            max_ast_ok=False,
            pt_provided=0,
            steps=steps
        )

    def _select_bars_with_layers(
        self,
        ast_required: float,
        width: float,
        preferred_bars: List[int],
    ) -> Tuple[int, int, int, float]:
        """
        Select bar diameter, count, and layers for bridge girders.

        Returns:
            Tuple of (bar_diameter, total_bar_count, num_layers, area_provided)
        """
        best_selection = None
        min_excess = float('inf')

        for dia in preferred_bars:
            if dia not in BAR_AREAS:
                continue

            bar_area = BAR_AREAS[dia]
            total_count = math.ceil(ast_required / bar_area)

            # Check how many bars fit per layer
            min_spacing = max(dia, 25)
            cover = 50
            bars_per_layer = int((width - 2*cover - 2*10 + min_spacing) / (dia + min_spacing))
            bars_per_layer = max(2, bars_per_layer)

            layers = math.ceil(total_count / bars_per_layer)
            # Round up total count to fill layers
            total_count = layers * bars_per_layer

            area_provided = total_count * bar_area
            excess = area_provided - ast_required

            if excess >= 0 and excess < min_excess:
                min_excess = excess
                best_selection = (dia, total_count, layers, area_provided)

        if best_selection is None:
            # Fallback
            dia = max(preferred_bars)
            bar_area = BAR_AREAS.get(dia, BAR_AREAS[32])
            total_count = math.ceil(ast_required / bar_area)
            layers = max(1, total_count // 4)
            area_provided = total_count * bar_area
            best_selection = (dia, total_count, layers, area_provided)

        return best_selection

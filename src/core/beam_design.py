"""
Main bridge girder design orchestrator per IRC 112:2020.

Coordinates the complete design workflow:
1. Dead load calculation (girder + deck + wearing coat)
2. IRC 6 live load analysis (Class A, 70R Wheeled, 70R Tracked)
3. Load combinations (ULS and SLS)
4. Preliminary sizing
5. Flexural design (T-beam or rectangular)
6. Shear design
7. SLS stress checks (critical for bridges!)
8. Deflection check
"""

import math
from typing import Tuple, List
from datetime import datetime

from src.codes.base_code import DesignCode
from src.codes.irc112 import (
    IRC112, irc_class_a_train, irc_class_70r_wheeled, irc_class_70r_tracked,
    impact_factor_class_a, impact_factor_70r_wheeled, impact_factor_70r_tracked,
    max_bm_simply_supported, max_sf_simply_supported,
    max_bm_tracked_vehicle, max_sf_tracked_vehicle
)
from src.models.inputs import BridgeGirderInput
from src.models.outputs import (
    BridgeGirderOutput, FlexuralDesignOutput, ShearDesignOutput,
    SLSStressCheckOutput, CrackWidthCheckOutput, DeflectionCheckOutput,
    LoadAnalysisOutput, LiveLoadResult, SideFaceReinforcementOutput,
    DesignStatus
)
from src.core.flexure import FlexureDesigner, FlexureResult
from src.core.shear import ShearDesigner, ShearResult
from src.core.serviceability import ServiceabilityChecker, SLSStressResult, CrackWidthResult, DeflectionResult


class BridgeGirderDesignEngine:
    """
    Main calculation engine for IRC 112 bridge girder design.

    Key features:
    - IRC 6 vehicle loading (Class A, 70R Wheeled, 70R Tracked)
    - T-beam design with deck slab as flange
    - ULS and SLS checks per IRC 112
    - SLS stress checks (critical for bridge design)
    """

    UNIT_WT_CONCRETE = 25  # kN/m³
    UNIT_WT_WEARING_COAT = 22  # kN/m³

    def __init__(self, code: DesignCode = None):
        self.code = code or IRC112()
        self.flexure_designer = FlexureDesigner(self.code)
        self.shear_designer = ShearDesigner(self.code)
        self.serviceability_checker = ServiceabilityChecker(self.code)

    def design(self, inputs: BridgeGirderInput) -> BridgeGirderOutput:
        """
        Execute complete bridge girder design workflow.

        Args:
            inputs: BridgeGirderInput with all parameters

        Returns:
            BridgeGirderOutput with complete design results
        """
        warnings = []
        notes = []

        # Material properties from simplified input
        fck = inputs.fck
        fy = inputs.fy
        fcd = 0.67 * fck / 1.5
        fyd = fy / 1.15
        ecm = self.code.get_ecm(fck)

        # Geometry from simplified input
        span_m = inputs.span_length
        girder_type = inputs.girder_type
        if isinstance(girder_type, str):
            is_tbeam = girder_type == "t_beam"
        else:
            is_tbeam = girder_type.value == "t_beam"

        # Get cover from input property
        cover = inputs.cover

        # Get dimensions from input properties
        bw = inputs.web_width
        D = inputs.overall_depth
        bf = inputs.flange_width if is_tbeam else bw
        Df = inputs.flange_depth if is_tbeam else 0

        # Effective depth
        bar_dia_main = self._select_main_bar_size(span_m)
        stirrup_dia = inputs.preferred_stirrup_size
        d = D - cover - stirrup_dia - bar_dia_main / 2

        # Dead load calculation
        load_analysis = self._calculate_loads(
            inputs, span_m, bw, D, bf, Df, is_tbeam
        )

        # Design moment and shear
        Mu = load_analysis.factored_moment
        Vu = load_analysis.factored_shear
        Ms = load_analysis.service_moment
        Mll = load_analysis.governing_live_moment

        # Flexural design
        reinf_pref = inputs.reinforcement_preference
        if isinstance(reinf_pref, str):
            allow_doubly = reinf_pref == "allow_doubly"
        else:
            allow_doubly = reinf_pref.value == "allow_doubly"

        flexure_result = self.flexure_designer.design(
            moment=Mu,
            web_width=bw,
            depth=D,
            effective_depth=d,
            fck=fck,
            fy=fy,
            is_tbeam=is_tbeam,
            flange_width=bf if is_tbeam else None,
            flange_depth=Df if is_tbeam else None,
            effective_cover=cover + stirrup_dia + bar_dia_main / 2,
            allow_doubly=allow_doubly,
            preferred_bars=inputs.preferred_bar_sizes,
        )

        # Update effective depth for multiple layers
        if flexure_result.tension_bar_layers > 1:
            layer_spacing = flexure_result.tension_bar_dia + 25
            d_eff = D - cover - stirrup_dia - flexure_result.tension_bar_dia / 2 - \
                    (flexure_result.tension_bar_layers - 1) * layer_spacing / 2
        else:
            d_eff = d

        # Shear design
        shear_result = self.shear_designer.design(
            shear_force=Vu,
            span_m=span_m,
            web_width=bw,
            effective_depth=d_eff,
            fck=fck,
            fy=fy,
            ast_provided=flexure_result.provided_ast,
            factored_udl=1.35 * load_analysis.total_dead_load,
            stirrup_dia=stirrup_dia,
        )

        # SLS stress check (CRITICAL for bridges)
        Mdl = load_analysis.dead_load_moment
        sls_result = self.serviceability_checker.check_sls_stresses(
            service_moment=Ms,
            web_width=bw,
            effective_depth=d_eff,
            ast_provided=flexure_result.provided_ast,
            fck=fck,
            fy=fy,
            ecm=ecm,
            is_tbeam=is_tbeam,
            flange_width=bf if is_tbeam else None,
            flange_depth=Df if is_tbeam else None,
            dead_load_moment=Mdl,  # For quasi-permanent calculation
        )

        if not sls_result.sigma_c_ok:
            warnings.append(
                f"SLS concrete stress ({sls_result.sigma_c:.1f} MPa) exceeds limit "
                f"({sls_result.sigma_c_limit:.1f} MPa). Consider higher concrete grade."
            )

        if not sls_result.sigma_c_qp_ok:
            warnings.append(
                f"Quasi-permanent concrete stress ({sls_result.sigma_c_qp:.1f} MPa) exceeds limit "
                f"({sls_result.sigma_c_qp_limit:.1f} MPa). Non-linear creep may occur."
            )

        # Crack width check
        # Get bar spacing from arrangement
        exposure = inputs.exposure_condition
        if not isinstance(exposure, str):
            exposure = exposure.value

        # Calculate bar spacing
        n_bars_per_layer = min(flexure_result.tension_bar_count,
                               int((bw - 2*cover - 2*stirrup_dia - flexure_result.tension_bar_dia) /
                                   (flexure_result.tension_bar_dia + 25)) + 1)
        bar_spacing = ((bw - 2*cover - 2*stirrup_dia - flexure_result.tension_bar_dia) /
                       max(1, n_bars_per_layer - 1)) if n_bars_per_layer > 1 else 0

        crack_width_result = self.serviceability_checker.check_crack_width(
            service_moment=sls_result.qp_moment,  # Use quasi-permanent moment
            web_width=bw,
            effective_depth=d_eff,
            overall_depth=D,
            ast_provided=flexure_result.provided_ast,
            bar_dia=flexure_result.tension_bar_dia,
            bar_spacing=bar_spacing if bar_spacing > 0 else flexure_result.tension_bar_dia + 25,
            clear_cover=cover,
            fck=fck,
            fy=fy,
            ecm=ecm,
            exposure=exposure,
        )

        if not crack_width_result.crack_width_ok:
            warnings.append(
                f"Crack width ({crack_width_result.wk_calculated:.3f} mm) exceeds limit "
                f"({crack_width_result.wk_limit:.1f} mm). Consider more bars or smaller diameter."
            )

        # Deflection check
        deflection_result = self.serviceability_checker.check_deflection(
            span_m=span_m,
            web_width=bw,
            overall_depth=D,
            effective_depth=d_eff,
            ast_provided=flexure_result.provided_ast,
            service_moment=Ms,
            total_dl=load_analysis.total_dead_load,
            live_moment=Mll,
            fck=fck,
            ecm=ecm,
            is_tbeam=is_tbeam,
            flange_width=bf if is_tbeam else None,
            flange_depth=Df if is_tbeam else None,
        )

        # Side face reinforcement for deep beams
        side_face = self._check_side_face_reinforcement(D, bw)

        # Overall status
        statuses = [
            flexure_result.status,
            shear_result.status,
            sls_result.status,
            crack_width_result.status,
            deflection_result.status,
        ]

        if DesignStatus.FAIL in statuses:
            overall_status = DesignStatus.FAIL
        elif DesignStatus.WARNING in statuses:
            overall_status = DesignStatus.WARNING
        else:
            overall_status = DesignStatus.PASS

        # Generate summary
        summary = self._generate_summary(
            bw, D, bf, Df, is_tbeam,
            flexure_result, shear_result, sls_result,
            overall_status, load_analysis
        )

        # Notes
        if is_tbeam:
            notes.append(f"T-beam design with flange width {bf:.0f}mm and depth {Df:.0f}mm")
        notes.append(f"Governing live load: {load_analysis.governing_bm_case}")

        # Convert to output models
        flexure_output = self._convert_flexure_output(flexure_result, Mu, bf, Df, is_tbeam)
        shear_output = self._convert_shear_output(shear_result, Vu)
        sls_output = self._convert_sls_output(sls_result)
        crack_width_output = self._convert_crack_width_output(crack_width_result)
        deflection_output = self._convert_deflection_output(deflection_result)

        # Get grade strings
        concrete_grade_str = inputs.concrete_grade
        steel_grade_str = inputs.steel_grade
        if not isinstance(concrete_grade_str, str):
            concrete_grade_str = concrete_grade_str.value
        if not isinstance(steel_grade_str, str):
            steel_grade_str = steel_grade_str.value

        return BridgeGirderOutput(
            overall_status=overall_status,
            design_summary=summary,
            span_length=span_m,
            girder_type="T-beam" if is_tbeam else "Rectangular",
            girder_width=bw,
            girder_depth=D,
            effective_depth=d_eff,
            clear_cover=cover,
            flange_width=bf if is_tbeam else None,
            flange_depth=Df if is_tbeam else None,
            concrete_grade=concrete_grade_str,
            steel_grade=steel_grade_str,
            fck=fck,
            fy=fy,
            fcd=fcd,
            fyd=fyd,
            ecm=ecm,
            loads=load_analysis,
            flexure=flexure_output,
            shear=shear_output,
            sls_stress=sls_output,
            crack_width=crack_width_output,
            deflection=deflection_output,
            side_face_reinforcement=side_face,
            design_code="IRC 112:2020",
            computed_at=datetime.now(),
            warnings=warnings,
            notes=notes,
        )

    def _select_main_bar_size(self, span_m: float) -> int:
        """Select main bar diameter based on span."""
        if span_m <= 10:
            return 25
        elif span_m <= 15:
            return 28
        elif span_m <= 20:
            return 32
        else:
            return 36

    def _calculate_loads(
        self,
        inputs: BridgeGirderInput,
        span_m: float,
        bw: float,
        D: float,
        bf: float,
        Df: float,
        is_tbeam: bool,
    ) -> LoadAnalysisOutput:
        """
        Calculate dead loads and live loads with IRC 6 analysis.
        """
        # Get deck parameters
        deck_width_m = inputs.deck.width  # Already in meters per girder
        deck_thickness_m = inputs.deck.thickness / 1000  # mm to m
        wc_thickness_m = inputs.deck.wearing_coat / 1000  # mm to m

        # Self weight of girder
        if is_tbeam:
            # Web portion below deck
            girder_depth_below_deck = (D - Df) / 1000
            sw_girder = (bw / 1000) * girder_depth_below_deck * self.UNIT_WT_CONCRETE
        else:
            sw_girder = (bw / 1000) * (D / 1000) * self.UNIT_WT_CONCRETE

        # Weight of deck slab (tributary width = deck width per girder)
        sw_deck = deck_width_m * deck_thickness_m * self.UNIT_WT_CONCRETE

        # Weight of wearing coat
        sw_wc = deck_width_m * wc_thickness_m * self.UNIT_WT_WEARING_COAT

        # Total dead load
        total_dl = sw_girder + sw_deck + sw_wc

        # Dead load BM and SF
        Mdl = total_dl * span_m ** 2 / 8
        Vdl = total_dl * span_m / 2

        # Live load analysis
        live_load_results = self._analyze_live_loads(span_m)

        # Find governing cases
        gov_bm_case = max(live_load_results, key=lambda x: x.max_bm_with_impact)
        gov_sf_case = max(live_load_results, key=lambda x: x.max_sf_with_impact)

        for result in live_load_results:
            if result.load_type == gov_bm_case.load_type:
                result.is_governing_bm = True
            if result.load_type == gov_sf_case.load_type:
                result.is_governing_sf = True

        Mll = gov_bm_case.max_bm_with_impact
        Vll = gov_sf_case.max_sf_with_impact

        # Load combinations (IRC 112 Table B.1)
        gamma_dl = 1.35
        gamma_ll = 1.50

        Mu = gamma_dl * Mdl + gamma_ll * Mll
        Vu = gamma_dl * Vdl + gamma_ll * Vll
        Ms = Mdl + Mll  # Service moment (rare combination)

        return LoadAnalysisOutput(
            sw_girder=round(sw_girder, 2),
            sw_deck=round(sw_deck, 2),
            sw_wearing_coat=round(sw_wc, 2),
            total_dead_load=round(total_dl, 2),
            dead_load_moment=round(Mdl, 2),
            dead_load_shear=round(Vdl, 2),
            live_load_results=live_load_results,
            governing_bm_case=gov_bm_case.load_type,
            governing_sf_case=gov_sf_case.load_type,
            governing_live_moment=round(Mll, 2),
            governing_live_shear=round(Vll, 2),
            factored_moment=round(Mu, 2),
            factored_shear=round(Vu, 2),
            service_moment=round(Ms, 2),
        )

    def _analyze_live_loads(
        self,
        span_m: float,
    ) -> List[LiveLoadResult]:
        """Analyze IRC 6 vehicle loads."""
        results = []

        # Class A
        class_a_loads = irc_class_a_train()
        Ma_class_a, _, _ = max_bm_simply_supported(span_m, class_a_loads)
        Va_class_a, _ = max_sf_simply_supported(span_m, class_a_loads)
        if_class_a = impact_factor_class_a(span_m)

        results.append(LiveLoadResult(
            load_type="Class A",
            max_bm_without_impact=round(Ma_class_a, 2),
            max_sf_without_impact=round(Va_class_a, 2),
            impact_factor=round(if_class_a, 3),
            max_bm_with_impact=round(Ma_class_a * (1 + if_class_a), 2),
            max_sf_with_impact=round(Va_class_a * (1 + if_class_a), 2),
        ))

        # 70R Wheeled
        class_70r_wh = irc_class_70r_wheeled()
        Ma_70r_wh, _, _ = max_bm_simply_supported(span_m, class_70r_wh)
        Va_70r_wh, _ = max_sf_simply_supported(span_m, class_70r_wh)
        if_70r_wh = impact_factor_70r_wheeled(span_m)

        results.append(LiveLoadResult(
            load_type="70R Wheeled",
            max_bm_without_impact=round(Ma_70r_wh, 2),
            max_sf_without_impact=round(Va_70r_wh, 2),
            impact_factor=round(if_70r_wh, 3),
            max_bm_with_impact=round(Ma_70r_wh * (1 + if_70r_wh), 2),
            max_sf_with_impact=round(Va_70r_wh * (1 + if_70r_wh), 2),
        ))

        # 70R Tracked
        tracked_load, contact_length = irc_class_70r_tracked()
        Ma_70r_tr = max_bm_tracked_vehicle(span_m, tracked_load, contact_length)
        Va_70r_tr = max_sf_tracked_vehicle(span_m, tracked_load, contact_length)
        if_70r_tr = impact_factor_70r_tracked(span_m)

        results.append(LiveLoadResult(
            load_type="70R Tracked",
            max_bm_without_impact=round(Ma_70r_tr, 2),
            max_sf_without_impact=round(Va_70r_tr, 2),
            impact_factor=round(if_70r_tr, 3),
            max_bm_with_impact=round(Ma_70r_tr * (1 + if_70r_tr), 2),
            max_sf_with_impact=round(Va_70r_tr * (1 + if_70r_tr), 2),
        ))

        return results

    def _check_side_face_reinforcement(
        self,
        D: float,
        bw: float,
    ) -> SideFaceReinforcementOutput:
        """Check if side face reinforcement is required (D > 750mm)."""
        if D <= 750:
            return SideFaceReinforcementOutput(required=False)

        # Minimum area per face: 0.05% of web area
        min_area = 0.0005 * bw * D
        bar_dia = 12
        bar_area = math.pi * bar_dia ** 2 / 4
        n_bars = max(2, math.ceil(min_area / bar_area))

        return SideFaceReinforcementOutput(
            required=True,
            min_area_per_face=round(min_area, 0),
            bar_diameter=bar_dia,
            bar_count_per_face=n_bars,
            arrangement=f"{n_bars}-{bar_dia}dia each face",
        )

    def _generate_summary(
        self,
        bw: float,
        D: float,
        bf: float,
        Df: float,
        is_tbeam: bool,
        flexure: FlexureResult,
        shear: ShearResult,
        sls: SLSStressResult,
        status: DesignStatus,
        loads: LoadAnalysisOutput,
    ) -> str:
        """Generate design summary."""
        if status == DesignStatus.PASS:
            status_text = "DESIGN SAFE"
        elif status == DesignStatus.WARNING:
            status_text = "DESIGN OK (WITH WARNINGS)"
        else:
            status_text = "DESIGN INADEQUATE"

        section = f"{bw:.0f}x{D:.0f}mm"
        if is_tbeam:
            section += f" (T-beam, bf={bf:.0f}mm)"

        summary = f"{status_text}\n"
        summary += f"Section: {section}\n"
        summary += f"Governing Live Load: {loads.governing_bm_case}\n"
        summary += f"Bottom Steel: {flexure.tension_arrangement}\n"
        summary += f"Stirrups: {shear.stirrup_legs}L-{shear.stirrup_dia}dia @ {shear.spacing_provided:.0f}mm\n"
        summary += f"SLS Concrete Stress: {sls.sigma_c:.1f}/{sls.sigma_c_limit:.1f} MPa "
        summary += f"({'OK' if sls.sigma_c_ok else 'FAILS'})"

        return summary

    def _convert_flexure_output(
        self,
        result: FlexureResult,
        Mu: float,
        bf: float,
        Df: float,
        is_tbeam: bool,
    ) -> FlexuralDesignOutput:
        """Convert flexure result to output model."""
        return FlexuralDesignOutput(
            status=result.status,
            design_moment=Mu,
            moment_capacity_limit=result.moment_capacity_limit,
            is_tbeam=is_tbeam,
            flange_width=bf if is_tbeam else None,
            flange_depth=Df if is_tbeam else None,
            moment_if_na_at_flange=result.moment_if_na_at_flange,
            na_in_flange=result.na_in_flange,
            required_ast=result.required_ast,
            provided_ast=result.provided_ast,
            tension_bar_diameter=result.tension_bar_dia,
            tension_bar_count=result.tension_bar_count,
            tension_bar_layers=result.tension_bar_layers,
            tension_bar_arrangement=result.tension_arrangement,
            is_doubly_reinforced=result.is_doubly,
            required_asc=result.required_asc,
            provided_asc=result.provided_asc,
            compression_bar_diameter=result.compression_bar_dia,
            compression_bar_count=result.compression_bar_count,
            compression_bar_arrangement=result.compression_arrangement,
            neutral_axis_depth=result.xu,
            xu_d_ratio=result.xu_d_ratio,
            xu_max=result.xu_max,
            lever_arm=result.lever_arm,
            moment_capacity=result.moment_capacity,
            utilization_ratio=Mu / result.moment_capacity if result.moment_capacity > 0 else 0,
            is_under_reinforced=result.is_under_reinforced,
            min_ast_satisfied=result.min_ast_ok,
            max_ast_satisfied=result.max_ast_ok,
            pt_provided=result.pt_provided,
            calculation_steps=result.steps,
        )

    def _convert_shear_output(
        self,
        result: ShearResult,
        Vu: float,
    ) -> ShearDesignOutput:
        """Convert shear result to output model."""
        return ShearDesignOutput(
            status=result.status,
            design_shear=Vu,
            shear_at_d_from_support=result.shear_at_d,
            size_effect_factor_k=result.size_effect_k,
            reinforcement_ratio_rho=result.rho_l,
            concrete_shear_capacity=result.VRdc,
            min_concrete_shear_capacity=result.VRdc_min,
            max_shear_capacity=result.VRd_max,
            shear_reinforcement_required=result.shear_reinf_required,
            stirrup_diameter=result.stirrup_dia,
            stirrup_legs=result.stirrup_legs,
            stirrup_area=result.Asv,
            spacing_required=result.spacing_required,
            spacing_provided=result.spacing_provided,
            spacing_max=result.spacing_max,
            shear_adequacy_check=result.shear_adequacy,
            min_shear_reinf_check=result.min_reinf_ok,
            calculation_steps=result.steps,
        )

    def _convert_sls_output(
        self,
        result: SLSStressResult,
    ) -> SLSStressCheckOutput:
        """Convert SLS result to output model."""
        return SLSStressCheckOutput(
            status=result.status,
            rare_combination_moment=result.rare_moment,
            modular_ratio=result.modular_ratio,
            elastic_na_depth=result.elastic_na_depth,
            cracked_moment_of_inertia=result.cracked_mi,
            concrete_stress=result.sigma_c,
            concrete_stress_limit=result.sigma_c_limit,
            concrete_stress_check=result.sigma_c_ok,
            concrete_utilization=result.sigma_c_utilization,
            steel_stress=result.sigma_s,
            steel_stress_limit=result.sigma_s_limit,
            steel_stress_check=result.sigma_s_ok,
            steel_utilization=result.sigma_s_utilization,
            qp_moment=result.qp_moment,
            qp_concrete_stress=result.sigma_c_qp,
            qp_concrete_stress_limit=result.sigma_c_qp_limit,
            qp_concrete_stress_check=result.sigma_c_qp_ok,
            calculation_steps=result.steps,
        )

    def _convert_crack_width_output(
        self,
        result: CrackWidthResult,
    ) -> CrackWidthCheckOutput:
        """Convert crack width result to output model."""
        return CrackWidthCheckOutput(
            status=result.status,
            steel_stress_qp=result.sigma_s,
            clear_cover=result.cover,
            bar_diameter=result.bar_dia,
            bar_spacing=result.spacing,
            effective_reinforcement_ratio=result.rho_eff,
            max_crack_spacing=result.sr_max,
            strain_difference=result.eps_sm_minus_eps_cm,
            calculated_crack_width=result.wk_calculated,
            allowable_crack_width=result.wk_limit,
            crack_width_check=result.crack_width_ok,
            calculation_steps=result.steps,
        )

    def _convert_deflection_output(
        self,
        result: DeflectionResult,
    ) -> DeflectionCheckOutput:
        """Convert deflection result to output model."""
        return DeflectionCheckOutput(
            status=result.status,
            gross_moment_of_inertia=result.gross_mi,
            centroid_from_bottom=result.centroid,
            cracking_moment=result.cracking_moment,
            effective_moment_of_inertia=result.effective_mi,
            calculated_deflection=result.calculated,
            allowable_deflection_total=result.allowable_total,
            allowable_deflection_live=result.allowable_live,
            deflection_check_passed=result.deflection_ok,
            calculation_steps=result.steps,
        )


# Backward compatibility alias
BeamDesignEngine = BridgeGirderDesignEngine

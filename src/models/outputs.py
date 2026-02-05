"""
Output data models for IRC 112 bridge girder design results.
"""

from pydantic import BaseModel
from typing import Optional
from enum import Enum
from datetime import datetime


class DesignStatus(str, Enum):
    """Status of a design check."""
    PASS = "pass"
    FAIL = "fail"
    WARNING = "warning"


class CalculationStep(BaseModel):
    """Single calculation step for transparency."""
    step_number: int
    description: str
    formula: str
    substitution: str
    result: float
    unit: str
    code_reference: Optional[str] = None


class LiveLoadResult(BaseModel):
    """Results from IRC vehicle load analysis."""
    load_type: str  # 'Class A', '70R Wheeled', '70R Tracked'
    max_bm_without_impact: float  # kNm
    max_sf_without_impact: float  # kN
    impact_factor: float
    max_bm_with_impact: float  # kNm
    max_sf_with_impact: float  # kN
    is_governing_bm: bool = False
    is_governing_sf: bool = False


class LoadAnalysisOutput(BaseModel):
    """Complete load analysis results."""
    # Dead loads
    sw_girder: float  # Self weight of girder (kN/m)
    sw_deck: float  # Weight of deck slab (kN/m)
    sw_wearing_coat: float  # Weight of wearing coat (kN/m)
    total_dead_load: float  # Total DL (kN/m)

    # Dead load effects
    dead_load_moment: float  # Mdl (kNm)
    dead_load_shear: float  # Vdl (kN)

    # Live load analysis
    live_load_results: list[LiveLoadResult]
    governing_bm_case: str
    governing_sf_case: str
    governing_live_moment: float  # Mll (kNm)
    governing_live_shear: float  # Vll (kN)

    # Factored loads (ULS)
    factored_moment: float  # Mu = 1.35*Mdl + 1.50*Mll (kNm)
    factored_shear: float  # Vu = 1.35*Vdl + 1.50*Vll (kN)

    # Service loads (SLS)
    service_moment: float  # Ms = Mdl + Mll (kNm)


class FlexuralDesignOutput(BaseModel):
    """Flexural design results for IRC 112."""
    status: DesignStatus

    # Design forces
    design_moment: float  # Mu in kNm
    moment_capacity_limit: float  # Mu_lim in kNm

    # T-beam analysis
    is_tbeam: bool = False
    flange_width: Optional[float] = None  # bf in mm
    flange_depth: Optional[float] = None  # Df in mm
    moment_if_na_at_flange: Optional[float] = None  # Mf in kNm
    na_in_flange: Optional[bool] = None

    # Reinforcement - Tension
    required_ast: float  # Required tension steel area (mm²)
    provided_ast: float  # Provided tension steel area (mm²)
    tension_bar_diameter: int  # mm
    tension_bar_count: int
    tension_bar_layers: int
    tension_bar_arrangement: str  # e.g., "8-32φ (2L×4)"

    # Reinforcement - Compression (for doubly reinforced)
    is_doubly_reinforced: bool = False
    required_asc: Optional[float] = None  # mm²
    provided_asc: Optional[float] = None  # mm²
    compression_bar_diameter: Optional[int] = None
    compression_bar_count: Optional[int] = None
    compression_bar_arrangement: Optional[str] = None

    # Section analysis
    neutral_axis_depth: float  # xu in mm
    xu_d_ratio: float
    xu_max: float  # xu,max in mm
    lever_arm: float  # z in mm
    moment_capacity: float  # Mu_provided in kNm
    utilization_ratio: float  # Mu_applied / Mu_capacity

    # Code checks
    is_under_reinforced: bool
    min_ast_satisfied: bool
    max_ast_satisfied: bool
    pt_provided: float  # Tension steel percentage

    # Calculation breakdown
    calculation_steps: list[CalculationStep]


class ShearDesignOutput(BaseModel):
    """Shear design results per IRC 112."""
    status: DesignStatus

    # Design forces
    design_shear: float  # Vu in kN
    shear_at_d_from_support: float  # Vu at d from support in kN

    # IRC 112 shear capacity
    size_effect_factor_k: float  # k = 1 + sqrt(200/d)
    reinforcement_ratio_rho: float  # ρl
    concrete_shear_capacity: float  # VRd,c in kN
    min_concrete_shear_capacity: float  # VRd,c,min in kN
    max_shear_capacity: float  # VRd,max in kN

    # Shear reinforcement
    shear_reinforcement_required: bool
    stirrup_diameter: int  # mm
    stirrup_legs: int
    stirrup_area: float  # Asv in mm²
    spacing_required: float  # mm
    spacing_provided: float  # mm
    spacing_max: float  # mm

    # Code checks
    shear_adequacy_check: bool  # Vu <= VRd,max
    min_shear_reinf_check: bool

    # Calculation breakdown
    calculation_steps: list[CalculationStep]


class SLSStressCheckOutput(BaseModel):
    """Service Limit State stress check results per IRC 112 Cl. 12.2.1."""
    status: DesignStatus

    # Service moment
    rare_combination_moment: float  # M_rare = Mdl + Mll (kNm)

    # Elastic analysis
    modular_ratio: float  # m = Es/Ecm
    elastic_na_depth: float  # x in mm
    cracked_moment_of_inertia: float  # Icr in mm⁴

    # Concrete stress (Rare combination)
    concrete_stress: float  # σc in MPa
    concrete_stress_limit: float  # 0.48*fck in MPa
    concrete_stress_check: bool
    concrete_utilization: float  # σc / σc_allow

    # Steel stress (Rare combination)
    steel_stress: float  # σs in MPa
    steel_stress_limit: float  # 0.8*fy in MPa
    steel_stress_check: bool
    steel_utilization: float  # σs / σs_allow

    # Quasi-permanent stress check (to avoid non-linear creep)
    qp_moment: Optional[float] = None  # M_qp = Mdl + 0.3*Mll (kNm)
    qp_concrete_stress: Optional[float] = None  # σc,qp in MPa
    qp_concrete_stress_limit: Optional[float] = None  # 0.36*fck in MPa
    qp_concrete_stress_check: Optional[bool] = None

    # Calculation breakdown
    calculation_steps: list[CalculationStep]


class CrackWidthCheckOutput(BaseModel):
    """Crack width check results per IRC 112 Cl. 12.3.4."""
    status: DesignStatus

    # Steel stress under quasi-permanent
    steel_stress_qp: float  # σs under QP loading (MPa)

    # Crack spacing parameters
    clear_cover: float  # c (mm)
    bar_diameter: float  # φ (mm)
    bar_spacing: float  # s (mm)
    effective_reinforcement_ratio: float  # ρp,eff

    # Maximum crack spacing
    max_crack_spacing: float  # sr,max (mm)

    # Strain difference
    strain_difference: float  # εsm - εcm

    # Crack width
    calculated_crack_width: float  # wk (mm)
    allowable_crack_width: float  # wk,max (mm)
    crack_width_check: bool

    # Calculation breakdown
    calculation_steps: list[CalculationStep]


class DeflectionCheckOutput(BaseModel):
    """Deflection check results per IRC 112."""
    status: DesignStatus

    # Section properties
    gross_moment_of_inertia: float  # Ig in mm⁴
    centroid_from_bottom: float  # mm
    cracking_moment: float  # Mcr in kNm
    effective_moment_of_inertia: float  # Ie in mm⁴

    # Deflection
    calculated_deflection: float  # mm
    allowable_deflection_total: float  # L/250 mm
    allowable_deflection_live: float  # L/800 mm
    deflection_check_passed: bool

    # Calculation breakdown
    calculation_steps: list[CalculationStep]


class SideFaceReinforcementOutput(BaseModel):
    """Side face reinforcement for deep beams (D > 750mm)."""
    required: bool
    min_area_per_face: Optional[float] = None  # mm²
    bar_diameter: Optional[int] = None  # mm
    bar_count_per_face: Optional[int] = None
    arrangement: Optional[str] = None


class BridgeGirderOutput(BaseModel):
    """Complete bridge girder design output per IRC 112."""

    # Overall status
    overall_status: DesignStatus
    design_summary: str

    # Input echo
    span_length: float  # m
    girder_type: str  # 'T-beam' or 'Rectangular'

    # Final dimensions
    girder_width: float  # bw in mm
    girder_depth: float  # D in mm
    effective_depth: float  # d in mm
    clear_cover: float  # mm
    flange_width: Optional[float] = None  # bf for T-beam
    flange_depth: Optional[float] = None  # Df for T-beam

    # Material properties used
    concrete_grade: str
    steel_grade: str
    fck: float  # MPa
    fy: float  # MPa
    fcd: float  # Design compressive strength (MPa)
    fyd: float  # Design yield strength (MPa)
    ecm: float  # Modulus of elasticity (MPa)

    # Load analysis
    loads: LoadAnalysisOutput

    # Design results
    flexure: FlexuralDesignOutput
    shear: ShearDesignOutput
    sls_stress: SLSStressCheckOutput
    crack_width: Optional[CrackWidthCheckOutput] = None
    deflection: DeflectionCheckOutput
    side_face_reinforcement: SideFaceReinforcementOutput

    # Metadata
    design_code: str = "IRC 112:2020"
    computed_at: datetime

    # Warnings and notes
    warnings: list[str]
    notes: list[str]

    @property
    def reinforcement_summary(self) -> str:
        """Quick summary of reinforcement."""
        summary = f"Bottom: {self.flexure.tension_bar_arrangement}"
        if self.flexure.is_doubly_reinforced:
            summary += f" | Top: {self.flexure.compression_bar_arrangement}"
        summary += f" | Stirrups: {self.shear.stirrup_legs}L-{self.shear.stirrup_diameter}φ @ {self.shear.spacing_provided:.0f}mm"
        if self.side_face_reinforcement.required:
            summary += f" | Side: {self.side_face_reinforcement.arrangement}"
        return summary

    @property
    def is_safe(self) -> bool:
        """Check if all design checks pass."""
        crack_ok = self.crack_width is None or self.crack_width.status == DesignStatus.PASS
        return (
            self.flexure.status == DesignStatus.PASS and
            self.shear.status == DesignStatus.PASS and
            self.sls_stress.status == DesignStatus.PASS and
            crack_ok and
            self.deflection.status == DesignStatus.PASS
        )


# Backward compatibility alias
BeamDesignOutput = BridgeGirderOutput

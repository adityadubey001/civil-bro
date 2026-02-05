"""
Input data models for bridge girder design using Pydantic for validation.
Updated for IRC 112:2020 bridge design.
"""

from pydantic import BaseModel, Field, field_validator
from typing import Optional, Literal, List
from enum import Enum


class ConcreteGrade(str, Enum):
    """Available concrete grades per IRC 112."""
    M25 = "M25"
    M30 = "M30"
    M35 = "M35"
    M40 = "M40"
    M45 = "M45"
    M50 = "M50"
    M55 = "M55"
    M60 = "M60"


class SteelGrade(str, Enum):
    """Available steel grades."""
    FE415 = "Fe415"
    FE500 = "Fe500"
    FE550 = "Fe550"


class ExposureCondition(str, Enum):
    """Exposure conditions per IRC 112 Table 14.2."""
    MODERATE = "moderate"
    SEVERE = "severe"
    VERY_SEVERE = "very_severe"
    EXTREME = "extreme"


class ReinforcementPreference(str, Enum):
    """User preference for beam reinforcement type."""
    SINGLY_ONLY = "singly_only"  # Increase section if Mu > Mu_lim
    ALLOW_DOUBLY = "allow_doubly"  # Use compression steel if needed


class LiveLoadType(str, Enum):
    """IRC 6 vehicle loading types."""
    CLASS_A = "class_a"
    CLASS_70R_WHEELED = "70r_wheeled"
    CLASS_70R_TRACKED = "70r_tracked"


class GirderType(str, Enum):
    """Type of girder cross-section."""
    RECTANGULAR = "rectangular"
    T_BEAM = "t_beam"


class DeckInput(BaseModel):
    """Deck slab parameters - simplified for UI."""
    width: float = Field(
        ...,
        gt=0,
        le=15,
        description="Effective deck width per girder in meters"
    )
    thickness: float = Field(
        default=240,
        ge=150,
        le=400,
        description="Deck slab thickness in mm"
    )
    wearing_coat: float = Field(
        default=75,
        ge=40,
        le=150,
        description="Wearing coat thickness in mm"
    )


class BridgeGirderInput(BaseModel):
    """Complete input model for bridge girder design per IRC 112.

    Simplified flat structure for UI integration.
    """
    # Span
    span_length: float = Field(
        ...,
        gt=0,
        le=40,
        description="Effective span in meters"
    )

    # Deck parameters
    deck: DeckInput
    num_girders: int = Field(
        default=2,
        ge=2,
        le=10,
        description="Number of main girders"
    )

    # Girder type
    girder_type: GirderType = GirderType.T_BEAM

    # Materials
    concrete_grade: ConcreteGrade = ConcreteGrade.M40
    steel_grade: SteelGrade = SteelGrade.FE500

    # Live load
    live_load_type: LiveLoadType = LiveLoadType.CLASS_70R_WHEELED

    # Exposure
    exposure_condition: ExposureCondition = ExposureCondition.MODERATE

    # Design preferences
    reinforcement_preference: ReinforcementPreference = ReinforcementPreference.ALLOW_DOUBLY
    preferred_bar_sizes: Optional[List[int]] = Field(
        default=[20, 25, 32],
        description="Preferred main bar diameters in mm"
    )
    preferred_stirrup_size: int = Field(
        default=12,
        ge=8,
        le=16,
        description="Preferred stirrup diameter in mm"
    )

    # Optional custom dimensions
    custom_web_width: Optional[float] = Field(
        None,
        ge=250,
        le=1000,
        description="Custom girder web width in mm"
    )
    custom_depth: Optional[float] = Field(
        None,
        ge=400,
        le=3000,
        description="Custom overall girder depth in mm"
    )

    class Config:
        use_enum_values = True

    @property
    def fck(self) -> float:
        """Characteristic compressive strength in MPa."""
        grade = self.concrete_grade
        if isinstance(grade, str):
            return float(grade[1:])
        return float(grade.value[1:])

    @property
    def fy(self) -> float:
        """Characteristic yield strength in MPa."""
        grade = self.steel_grade
        if isinstance(grade, str):
            return float(grade[2:])
        return float(grade.value[2:])

    @property
    def cover(self) -> float:
        """Get nominal cover based on exposure condition."""
        exposure = self.exposure_condition
        if isinstance(exposure, str):
            exp_val = exposure
        else:
            exp_val = exposure.value

        cover_map = {
            "moderate": 40,
            "severe": 50,
            "very_severe": 60,
            "extreme": 75
        }
        return cover_map.get(exp_val, 40)

    @property
    def web_width(self) -> float:
        """Get web width (custom or auto-calculated)."""
        if self.custom_web_width:
            return self.custom_web_width
        # Auto: D/2.5, min 300mm
        depth = self.overall_depth
        width = max(300, depth / 2.5)
        return round(width / 25) * 25  # Round to 25mm

    @property
    def overall_depth(self) -> float:
        """Get overall depth (custom or auto-calculated)."""
        if self.custom_depth:
            return self.custom_depth
        # Auto: L/12 for T-beam, L/10 for rectangular
        girder = self.girder_type
        if isinstance(girder, str):
            is_tbeam = girder == "t_beam"
        else:
            is_tbeam = girder == GirderType.T_BEAM

        ratio = 12 if is_tbeam else 10
        depth = (self.span_length * 1000) / ratio
        return round(depth / 50) * 50  # Round to 50mm

    @property
    def flange_width(self) -> float:
        """Get effective flange width for T-beam."""
        # Effective width = deck width per girder (in mm)
        return self.deck.width * 1000

    @property
    def flange_depth(self) -> float:
        """Get flange depth (deck slab thickness)."""
        return self.deck.thickness


# Keep backward compatibility
BeamDesignInput = BridgeGirderInput

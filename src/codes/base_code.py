"""
Abstract base class for design code provisions.
Enables extensibility for different codes (IS 456, IRC 112, etc.)
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Dict


@dataclass
class CoverRequirements:
    """Minimum cover requirements based on exposure."""
    nominal_cover: float  # mm
    exposure_class: str
    min_grade: str  # Minimum concrete grade


class DesignCode(ABC):
    """
    Abstract base class for structural design codes.

    Purpose:
    - Define interface for code-specific provisions
    - Enable switching between IS 456, IRC 112, etc.
    - Centralize code clause references
    """

    @property
    @abstractmethod
    def code_name(self) -> str:
        """Return the code name/version."""
        pass

    @abstractmethod
    def get_partial_safety_factors(self) -> Dict[str, float]:
        """Return partial safety factors for materials and loads."""
        pass

    @abstractmethod
    def get_span_depth_ratio(self, support_type: str) -> float:
        """Return basic span/effective depth ratio per code."""
        pass

    @abstractmethod
    def get_minimum_cover(self, exposure_class: str) -> CoverRequirements:
        """Return minimum cover requirement for exposure condition."""
        pass

    @abstractmethod
    def get_minimum_reinforcement_ratio(self, fy: float) -> float:
        """Return minimum reinforcement percentage."""
        pass

    @abstractmethod
    def get_maximum_reinforcement_ratio(self) -> float:
        """Return maximum reinforcement percentage."""
        pass

    @abstractmethod
    def get_shear_strength_concrete(self, pt: float, fck: float) -> float:
        """Return design shear strength of concrete (τc)."""
        pass

    @abstractmethod
    def get_maximum_shear_stress(self, fck: float) -> float:
        """Return maximum shear stress limit (τc,max)."""
        pass

    @abstractmethod
    def get_xu_max_ratio(self, fy: float) -> float:
        """Return xu_max/d ratio for balanced section."""
        pass

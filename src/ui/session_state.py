"""
Session state management for IRC 112 bridge girder design Q&A flow.
"""

import streamlit as st
from dataclasses import dataclass, field
from typing import Optional, Any, List
from enum import IntEnum


class DesignStep(IntEnum):
    """Steps in the Q&A design flow."""
    BRIDGE_GEOMETRY = 1
    DECK_PARAMETERS = 2
    GIRDER_TYPE = 3
    MATERIALS = 4
    LIVE_LOAD = 5
    EXPOSURE = 6
    DESIGN_OPTIONS = 7
    REVIEW = 8
    RESULTS = 9


@dataclass
class DesignState:
    """Manages application state across the Q&A flow."""

    current_step: DesignStep = DesignStep.BRIDGE_GEOMETRY

    # Step 1: Bridge Geometry
    span_length: Optional[float] = None  # meters
    support_type: str = "simply_supported"

    # Step 2: Deck Parameters
    deck_width: Optional[float] = None  # meters (effective width per girder)
    deck_thickness: Optional[float] = None  # mm
    wearing_coat_thickness: Optional[float] = None  # mm
    num_girders: int = 2

    # Step 3: Girder Type
    girder_type: str = "t_beam"  # "t_beam" or "rectangular"
    custom_web_width: Optional[float] = None  # mm
    custom_depth: Optional[float] = None  # mm

    # Step 4: Materials
    concrete_grade: str = "M40"
    steel_grade: str = "Fe500"

    # Step 5: Live Load
    live_load_type: str = "70r_wheeled"  # "class_a", "70r_wheeled", "70r_tracked"
    num_lanes: int = 1

    # Step 6: Exposure
    exposure_condition: str = "moderate"

    # Step 7: Design Options
    reinforcement_preference: str = "allow_doubly"
    preferred_bar_sizes: List[int] = field(default_factory=lambda: [20, 25, 32])
    preferred_stirrup_size: int = 12

    # Results
    design_output: Optional[Any] = None

    def is_step_complete(self, step: DesignStep) -> bool:
        """Check if a step has all required inputs."""
        validators = {
            DesignStep.BRIDGE_GEOMETRY: lambda: (
                self.span_length is not None and
                self.span_length > 0
            ),
            DesignStep.DECK_PARAMETERS: lambda: (
                self.deck_width is not None and
                self.deck_thickness is not None and
                self.deck_width > 0 and
                self.deck_thickness > 0
            ),
            DesignStep.GIRDER_TYPE: lambda: self.girder_type in ["t_beam", "rectangular"],
            DesignStep.MATERIALS: lambda: (
                self.concrete_grade is not None and
                self.steel_grade is not None
            ),
            DesignStep.LIVE_LOAD: lambda: (
                self.live_load_type in ["class_a", "70r_wheeled", "70r_tracked"]
            ),
            DesignStep.EXPOSURE: lambda: self.exposure_condition is not None,
            DesignStep.DESIGN_OPTIONS: lambda: True,  # Has defaults
            DesignStep.REVIEW: lambda: True,
            DesignStep.RESULTS: lambda: self.design_output is not None,
        }
        return validators.get(step, lambda: True)()

    def can_proceed_to(self, step: DesignStep) -> bool:
        """Check if user can navigate to a step."""
        for s in DesignStep:
            if s < step and not self.is_step_complete(s):
                return False
        return True

    def get_progress(self) -> float:
        """Get progress percentage (0-1)."""
        completed = sum(1 for s in DesignStep if s < DesignStep.RESULTS and self.is_step_complete(s))
        total = len(DesignStep) - 1  # Exclude RESULTS
        return completed / total

    def reset(self):
        """Reset all inputs to defaults."""
        self.current_step = DesignStep.BRIDGE_GEOMETRY
        self.span_length = None
        self.support_type = "simply_supported"
        self.deck_width = None
        self.deck_thickness = None
        self.wearing_coat_thickness = None
        self.num_girders = 2
        self.girder_type = "t_beam"
        self.custom_web_width = None
        self.custom_depth = None
        self.concrete_grade = "M40"
        self.steel_grade = "Fe500"
        self.live_load_type = "70r_wheeled"
        self.num_lanes = 1
        self.exposure_condition = "moderate"
        self.reinforcement_preference = "allow_doubly"
        self.preferred_bar_sizes = [20, 25, 32]
        self.preferred_stirrup_size = 12
        self.design_output = None


def get_state() -> DesignState:
    """Get or initialize session state."""
    if 'design_state' not in st.session_state:
        st.session_state.design_state = DesignState()
    return st.session_state.design_state


def set_step(step: DesignStep):
    """Set the current step."""
    state = get_state()
    state.current_step = step

"""Session state management for the Substructure Design app.

Manages wizard step flow, stores all input parameters with defaults,
and provides helpers for building the config dict consumed by the engine.
"""

from __future__ import annotations

import streamlit as st

# Total wizard steps
NUM_STEPS = 8

STEP_NAMES = {
    1: "Project & Levels",
    2: "Superstructure",
    3: "Pier",
    4: "Pier Cap & Bearings",
    5: "Foundation",
    6: "Materials & Creep",
    7: "Seismic & Wind",
    8: "Review & Run",
}

# ---------------------------------------------------------------------------
# Default values (from sample_input.yaml)
# ---------------------------------------------------------------------------

_DEFAULTS = {
    # -- control --
    "input_mode": "yaml",
    "current_step": 1,
    "results": None,
    "config": None,

    # -- Step 1: Project & Levels --
    "project_name": "MY-PROJECT-001",
    "pier_id": "P1",
    "designer": "",
    "checker": "",
    "project_date": "2025-01-01",
    "frl": 315.685,
    "gl": 303.236,
    "gwt": 303.236,
    "hfl": 303.236,

    # -- Step 2: Superstructure --
    "deck_width": 12.0,
    "super_type": "U Girder",
    "left_span": 58.35,
    "right_span": 58.35,
    "num_girders": 2,
    "girder_spacing": 6.0,
    "cg_below_deck_top": 1.023,
    "depth_incl_slab": 3.23,
    "deck_slab_thickness": 0.25,
    "girder_flange_width": 0.3,
    "girder_inertia": 3.245,
    "continuity": "Deck Continuity",
    "num_continuous_spans": 2,
    "camber_superelevation": 2.5,
    "bearing_pedestal_depth": 0.5,
    "skew_angle": 0.0,
    "radius_of_curvature": 2000.0,
    "wearing_coat_thickness": 45.0,
    "crash_barrier_width": 0.5,
    "median_barrier_width": 0.0,
    "crash_barrier_height": 1.1,
    "noise_barrier_height": 1.5,
    "dl_reaction_lhs": 0.0,
    "dl_reaction_rhs": 0.0,
    "sidl_reaction_lhs": 0.0,
    "sidl_reaction_rhs": 0.0,

    # -- Step 3: Pier --
    "pier_type": "Circular",
    "pier_diameter_bottom": 2.0,
    "pier_diameter_top": 2.0,
    "pier_end_condition": "fixed-guided",
    "pier_n_bars": 90,
    "pier_bar_dia": 20,

    # -- Step 4: Pier Cap & Bearings --
    "piercap_width_long": 2.85,
    "piercap_length_trans": 9.2,
    "piercap_depth_max": 1.75,
    "piercap_depth_min": 0.75,
    "bearing_offset_long": 0.825,
    "bearing_type": "Elastomeric",
    "bearings_lhs": [[-0.825, -3.7], [-0.825, -2.3], [-0.825, 2.3], [-0.825, 3.7]],
    "bearings_rhs": [[0.825, -3.7], [0.825, -2.3], [0.825, 2.3], [0.825, 3.7]],
    "pedestal_length": 0.8,
    "pedestal_width": 0.8,
    "pedestal_beyond": 0.2,

    # -- Step 5: Foundation --
    "pile_diameter": 1.0,
    "piles_long": 2,
    "piles_trans": 2,
    "pile_spacing_factor": 3,
    "pilecap_edge_clearance": 150,
    "pilecap_thickness": 1.8,
    "pilecap_top_below_gl": 500,
    "total_pile_depth": 15.0,
    "pile_capacity_top": 550.0,
    "pile_uplift_capacity": -187.0,
    "soil_density": 20.0,
    "friction_angle": 30.0,
    "subgrade_modulus": 1.4,
    "fixity_ratio": 2.3,
    "pile_cover": 75,
    "pile_n_bars": 18,
    "pile_bar_dia": 22.6,

    # -- Step 6: Materials & Creep --
    "fck_pier": 50,
    "aggregate_pier": "Quartzite",
    "fck_piercap": 45,
    "aggregate_piercap": "Quartzite",
    "fck_pedestal": 45,
    "fck_pile_pilecap": 35,
    "fyk": 550,
    "Es": 200000,
    "concrete_density": 25,
    "exposure": "Severe",
    "water_density": 10,
    "age_of_loading": 7,
    "relative_humidity": 54.5,
    "design_life": 36500,

    # -- Step 7: Seismic & Wind --
    "seismic_zone": "IV",
    "importance_factor": 1.2,
    "response_reduction": 4.0,
    "soil_type": "medium",
    "wind_speed": 44.0,
    "terrain_category": 2,

    # -- Step 8: BOQ rates --
    "concrete_rate": 8000,
    "steel_rate": 70000,
}


def init_state() -> None:
    """Initialize session state with defaults if not already set."""
    for key, val in _DEFAULTS.items():
        if key not in st.session_state:
            # Deep-copy lists to avoid shared references
            if isinstance(val, list):
                st.session_state[key] = [row[:] for row in val]
            else:
                st.session_state[key] = val


def reset_state() -> None:
    """Clear all state and reset to defaults."""
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    init_state()


def get_progress() -> float:
    """Return wizard progress as 0.0-1.0."""
    step = st.session_state.get("current_step", 1)
    return (step - 1) / NUM_STEPS


def build_config_from_state() -> dict:
    """Assemble a YAML-compatible config dict from session state values.

    This produces the same structure as sample_input.yaml so the engine
    modules can consume it directly.
    """
    s = st.session_state
    return {
        "project": {
            "name": s.get("project_name", ""),
            "pier_id": s.get("pier_id", ""),
            "designer": s.get("designer", ""),
            "checker": s.get("checker", ""),
            "date": s.get("project_date", ""),
        },
        "superstructure": {
            "deck_width": s.get("deck_width"),
            "type": s.get("super_type"),
            "left_span_cc_bearing": s.get("left_span"),
            "right_span_cc_bearing": s.get("right_span"),
            "num_girders": s.get("num_girders"),
            "girder_spacing": s.get("girder_spacing"),
            "cg_below_deck_top": s.get("cg_below_deck_top"),
            "depth_incl_slab": s.get("depth_incl_slab"),
            "deck_slab_thickness": s.get("deck_slab_thickness"),
            "girder_flange_width": s.get("girder_flange_width"),
            "girder_inertia": s.get("girder_inertia"),
            "continuity": s.get("continuity"),
            "num_continuous_spans": s.get("num_continuous_spans"),
            "camber_superelevation": s.get("camber_superelevation"),
            "bearing_pedestal_depth": s.get("bearing_pedestal_depth"),
            "skew_angle": s.get("skew_angle"),
            "radius_of_curvature": s.get("radius_of_curvature"),
            "wearing_coat_thickness": s.get("wearing_coat_thickness"),
            "crash_barrier_width": s.get("crash_barrier_width"),
            "median_barrier_width": s.get("median_barrier_width"),
            "crash_barrier_height": s.get("crash_barrier_height"),
            "noise_barrier_height": s.get("noise_barrier_height"),
            "dl_reaction_lhs": s.get("dl_reaction_lhs"),
            "dl_reaction_rhs": s.get("dl_reaction_rhs"),
            "sidl_reaction_lhs": s.get("sidl_reaction_lhs"),
            "sidl_reaction_rhs": s.get("sidl_reaction_rhs"),
        },
        "pier": {
            "type": s.get("pier_type"),
            "diameter_bottom": s.get("pier_diameter_bottom"),
            "diameter_top": s.get("pier_diameter_top"),
            "end_condition": s.get("pier_end_condition"),
            "n_bars": s.get("pier_n_bars"),
            "bar_dia": s.get("pier_bar_dia"),
        },
        "pier_cap": {
            "width_long": s.get("piercap_width_long"),
            "length_trans": s.get("piercap_length_trans"),
            "depth_max": s.get("piercap_depth_max"),
            "depth_min": s.get("piercap_depth_min"),
            "bearing_offset_long": s.get("bearing_offset_long"),
        },
        "bearings": {
            "type": s.get("bearing_type"),
            "coordinates_lhs": s.get("bearings_lhs"),
            "coordinates_rhs": s.get("bearings_rhs"),
            "pedestal_length": s.get("pedestal_length"),
            "pedestal_width": s.get("pedestal_width"),
            "pedestal_beyond": s.get("pedestal_beyond"),
        },
        "levels": {
            "frl": s.get("frl"),
            "gl": s.get("gl"),
            "gwt": s.get("gwt"),
            "hfl": s.get("hfl"),
        },
        "materials": {
            "pier": {
                "fck": s.get("fck_pier"),
                "aggregate": s.get("aggregate_pier"),
            },
            "pier_cap": {
                "fck": s.get("fck_piercap"),
                "aggregate": s.get("aggregate_piercap"),
            },
            "pedestal": {
                "fck": s.get("fck_pedestal"),
            },
            "pile_pilecap": {
                "fck": s.get("fck_pile_pilecap"),
            },
            "steel": {
                "fyk": s.get("fyk"),
                "Es": s.get("Es"),
            },
            "concrete_density": s.get("concrete_density"),
            "exposure": s.get("exposure"),
            "water_density": s.get("water_density"),
        },
        "creep": {
            "age_of_loading": s.get("age_of_loading"),
            "relative_humidity": s.get("relative_humidity"),
            "design_life": s.get("design_life"),
        },
        "foundation": {
            "pile_diameter": s.get("pile_diameter"),
            "piles_long": s.get("piles_long"),
            "piles_trans": s.get("piles_trans"),
            "pile_spacing_factor": s.get("pile_spacing_factor"),
            "pilecap_edge_clearance": s.get("pilecap_edge_clearance"),
            "pilecap_thickness": s.get("pilecap_thickness"),
            "pilecap_top_below_gl": s.get("pilecap_top_below_gl"),
            "total_pile_depth": s.get("total_pile_depth"),
            "pile_capacity_top": s.get("pile_capacity_top"),
            "pile_uplift_capacity": s.get("pile_uplift_capacity"),
            "soil_density": s.get("soil_density"),
            "friction_angle": s.get("friction_angle"),
            "subgrade_modulus": s.get("subgrade_modulus"),
            "fixity_ratio": s.get("fixity_ratio"),
            "cover": s.get("pile_cover"),
            "pile_n_bars": s.get("pile_n_bars"),
            "pile_bar_dia": s.get("pile_bar_dia"),
        },
        "seismic": {
            "zone": s.get("seismic_zone"),
            "importance_factor": s.get("importance_factor"),
            "response_reduction": s.get("response_reduction"),
            "soil_type": s.get("soil_type"),
        },
        "wind": {
            "basic_speed": s.get("wind_speed"),
            "terrain_category": s.get("terrain_category"),
        },
        "boq": {
            "concrete_rate": s.get("concrete_rate"),
            "steel_rate": s.get("steel_rate"),
        },
    }

"""Step-by-step wizard input mode for the Substructure Design app.

Guides users through 8 steps to collect all design parameters,
then assembles a config dict and runs the design.
"""

from __future__ import annotations

import streamlit as st

from ui.state import NUM_STEPS, STEP_NAMES, build_config_from_state
from ui.runner import run_design


# ── Main entry point ─────────────────────────────────────────────────────────

def render_wizard() -> None:
    """Render the current wizard step."""
    step = st.session_state.get("current_step", 1)
    st.markdown(f"## Step {step} of {NUM_STEPS}: {STEP_NAMES.get(step, '')}")

    renderers = {
        1: _step_project,
        2: _step_superstructure,
        3: _step_pier,
        4: _step_piercap_bearings,
        5: _step_foundation,
        6: _step_materials,
        7: _step_seismic_wind,
        8: _step_review,
    }

    renderer = renderers.get(step, _step_project)
    renderer()


# ── Navigation helpers ───────────────────────────────────────────────────────

def _nav_buttons(show_back: bool = True, show_next: bool = True,
                 next_label: str = "Next \u2192") -> None:
    """Render Back / Next navigation buttons."""
    cols = st.columns([1, 1]) if show_back else [st.columns(1)[0]]
    if show_back:
        with cols[0]:
            if st.button("\u2190 Back", use_container_width=True):
                st.session_state["current_step"] = max(1, st.session_state["current_step"] - 1)
                st.rerun()
    if show_next:
        with cols[-1]:
            if st.button(next_label, type="primary", use_container_width=True):
                st.session_state["current_step"] = min(NUM_STEPS, st.session_state["current_step"] + 1)
                st.rerun()


# ── Step 1: Project & Levels ─────────────────────────────────────────────────

def _step_project() -> None:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.session_state["project_name"] = st.text_input(
            "Project Name", value=st.session_state.get("project_name", ""))
        st.session_state["pier_id"] = st.text_input(
            "Pier ID", value=st.session_state.get("pier_id", ""))

        c1, c2 = st.columns(2)
        with c1:
            st.session_state["designer"] = st.text_input(
                "Designer (optional)", value=st.session_state.get("designer", ""))
        with c2:
            st.session_state["checker"] = st.text_input(
                "Checker (optional)", value=st.session_state.get("checker", ""))

        st.markdown("#### Levels (m)")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["frl"] = st.number_input(
                "FRL (Finished Road Level)", value=st.session_state.get("frl", 315.685),
                format="%.3f", step=0.1)
            st.session_state["gwt"] = st.number_input(
                "GWT (Ground Water Table)", value=st.session_state.get("gwt", 303.236),
                format="%.3f", step=0.1)
        with c2:
            st.session_state["gl"] = st.number_input(
                "GL (Ground Level)", value=st.session_state.get("gl", 303.236),
                format="%.3f", step=0.1)
            st.session_state["hfl"] = st.number_input(
                "HFL (High Flood Level)", value=st.session_state.get("hfl", 303.236),
                format="%.3f", step=0.1)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Levels** are used to compute pier height, "
            "submersion depth for hydrodynamic loads, "
            "and pile cap position.\n\n"
            "FRL must be above GL."
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons(show_back=False)


# ── Step 2: Superstructure ───────────────────────────────────────────────────

def _step_superstructure() -> None:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("#### Primary Parameters")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["deck_width"] = st.number_input(
                "Deck Width (m)", value=st.session_state.get("deck_width", 12.0),
                min_value=3.0, max_value=30.0, step=0.5)
            st.session_state["left_span"] = st.number_input(
                "Left Span c/c Bearing (m)", value=st.session_state.get("left_span", 58.35),
                min_value=5.0, step=0.5)
            st.session_state["depth_incl_slab"] = st.number_input(
                "Depth incl. Slab (m)", value=st.session_state.get("depth_incl_slab", 3.23),
                min_value=0.3, step=0.1)
        with c2:
            st.session_state["super_type"] = st.selectbox(
                "Girder Type",
                ["U Girder", "Box Girder", "Slab", "I Girder"],
                index=["U Girder", "Box Girder", "Slab", "I Girder"].index(
                    st.session_state.get("super_type", "U Girder")))
            st.session_state["right_span"] = st.number_input(
                "Right Span c/c Bearing (m)", value=st.session_state.get("right_span", 58.35),
                min_value=5.0, step=0.5)
            st.session_state["num_girders"] = st.number_input(
                "Number of Girders", value=st.session_state.get("num_girders", 2),
                min_value=1, max_value=10, step=1)

        with st.expander("Advanced Superstructure Parameters"):
            c1, c2, c3 = st.columns(3)
            with c1:
                st.session_state["girder_spacing"] = st.number_input(
                    "Girder Spacing (m)", value=st.session_state.get("girder_spacing", 6.0), step=0.5)
                st.session_state["cg_below_deck_top"] = st.number_input(
                    "CG below Deck Top (m)", value=st.session_state.get("cg_below_deck_top", 1.023), step=0.1)
                st.session_state["deck_slab_thickness"] = st.number_input(
                    "Deck Slab Thickness (m)", value=st.session_state.get("deck_slab_thickness", 0.25), step=0.05)
                st.session_state["girder_flange_width"] = st.number_input(
                    "Girder Flange Width (m)", value=st.session_state.get("girder_flange_width", 0.3), step=0.05)
                st.session_state["girder_inertia"] = st.number_input(
                    "Girder Inertia (m\u2074)", value=st.session_state.get("girder_inertia", 3.245), step=0.1, format="%.3f")
            with c2:
                st.session_state["continuity"] = st.selectbox(
                    "Continuity",
                    ["Deck Continuity", "Simply Supported"],
                    index=0 if st.session_state.get("continuity") == "Deck Continuity" else 1)
                st.session_state["num_continuous_spans"] = st.number_input(
                    "Continuous Spans", value=st.session_state.get("num_continuous_spans", 2),
                    min_value=1, step=1)
                st.session_state["camber_superelevation"] = st.number_input(
                    "Camber / Superelevation (%)", value=st.session_state.get("camber_superelevation", 2.5), step=0.5)
                st.session_state["bearing_pedestal_depth"] = st.number_input(
                    "Bearing Pedestal Depth (m)", value=st.session_state.get("bearing_pedestal_depth", 0.5), step=0.1)
                st.session_state["skew_angle"] = st.number_input(
                    "Skew Angle (\u00b0)", value=st.session_state.get("skew_angle", 0.0),
                    min_value=0.0, max_value=90.0, step=5.0)
            with c3:
                st.session_state["radius_of_curvature"] = st.number_input(
                    "Radius of Curvature (m)", value=st.session_state.get("radius_of_curvature", 2000.0), step=100.0,
                    help="Use 0 for straight bridges")
                st.session_state["wearing_coat_thickness"] = st.number_input(
                    "Wearing Coat (mm)", value=st.session_state.get("wearing_coat_thickness", 45.0), step=5.0)
                st.session_state["crash_barrier_width"] = st.number_input(
                    "Crash Barrier Width (m)", value=st.session_state.get("crash_barrier_width", 0.5), step=0.1)
                st.session_state["crash_barrier_height"] = st.number_input(
                    "Crash Barrier Height (m)", value=st.session_state.get("crash_barrier_height", 1.1), step=0.1)
                st.session_state["median_barrier_width"] = st.number_input(
                    "Median Barrier Width (m)", value=st.session_state.get("median_barrier_width", 0.0), step=0.1)
                st.session_state["noise_barrier_height"] = st.number_input(
                    "Noise Barrier Height (m)", value=st.session_state.get("noise_barrier_height", 1.5), step=0.1)

        with st.expander("DL / SIDL Reaction Overrides (set 0 for auto)"):
            c1, c2 = st.columns(2)
            with c1:
                st.session_state["dl_reaction_lhs"] = st.number_input(
                    "DL Reaction LHS (kN)", value=st.session_state.get("dl_reaction_lhs", 0.0), step=10.0)
                st.session_state["sidl_reaction_lhs"] = st.number_input(
                    "SIDL Reaction LHS (kN)", value=st.session_state.get("sidl_reaction_lhs", 0.0), step=10.0)
            with c2:
                st.session_state["dl_reaction_rhs"] = st.number_input(
                    "DL Reaction RHS (kN)", value=st.session_state.get("dl_reaction_rhs", 0.0), step=10.0)
                st.session_state["sidl_reaction_rhs"] = st.number_input(
                    "SIDL Reaction RHS (kN)", value=st.session_state.get("sidl_reaction_rhs", 0.0), step=10.0)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Superstructure** parameters define the deck geometry "
            "and load characteristics.\n\n"
            "The tool auto-calculates DL from girder type if "
            "reaction overrides are set to 0."
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons()


# ── Step 3: Pier ─────────────────────────────────────────────────────────────

def _step_pier() -> None:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.session_state["pier_type"] = st.selectbox(
            "Pier Type",
            ["Circular", "Rectangular"],
            index=0 if st.session_state.get("pier_type") == "Circular" else 1)

        c1, c2 = st.columns(2)
        with c1:
            label_bottom = "Diameter Bottom (m)" if st.session_state["pier_type"] == "Circular" else "Width Long (m)"
            st.session_state["pier_diameter_bottom"] = st.number_input(
                label_bottom, value=st.session_state.get("pier_diameter_bottom", 2.0),
                min_value=0.5, step=0.1)
        with c2:
            label_top = "Diameter Top (m)" if st.session_state["pier_type"] == "Circular" else "Length Trans (m)"
            st.session_state["pier_diameter_top"] = st.number_input(
                label_top, value=st.session_state.get("pier_diameter_top", 2.0),
                min_value=0.5, step=0.1)

        st.session_state["pier_end_condition"] = st.selectbox(
            "End Condition",
            ["fixed-free", "fixed-guided", "fixed-fixed"],
            index=["fixed-free", "fixed-guided", "fixed-fixed"].index(
                st.session_state.get("pier_end_condition", "fixed-guided")))

        st.markdown("#### Reinforcement")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["pier_n_bars"] = st.number_input(
                "Number of Bars", value=st.session_state.get("pier_n_bars", 90),
                min_value=8, step=2)
        with c2:
            st.session_state["pier_bar_dia"] = st.number_input(
                "Bar Diameter (mm)", value=st.session_state.get("pier_bar_dia", 20),
                min_value=10, max_value=40, step=2)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Pier** is the vertical column supporting the superstructure.\n\n"
            "**End conditions** affect effective length:\n"
            "- Fixed-free (cantilever): l0 = 2.0\u00d7L\n"
            "- Fixed-guided: l0 = 1.0\u00d7L\n"
            "- Fixed-fixed: l0 = 0.7\u00d7L"
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons()


# ── Step 4: Pier Cap & Bearings ──────────────────────────────────────────────

def _step_piercap_bearings() -> None:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("#### Pier Cap Dimensions")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["piercap_width_long"] = st.number_input(
                "Width Longitudinal (m)", value=st.session_state.get("piercap_width_long", 2.85),
                min_value=0.5, step=0.1)
            st.session_state["piercap_depth_max"] = st.number_input(
                "Depth Max (m)", value=st.session_state.get("piercap_depth_max", 1.75),
                min_value=0.3, step=0.1)
        with c2:
            st.session_state["piercap_length_trans"] = st.number_input(
                "Length Transverse (m)", value=st.session_state.get("piercap_length_trans", 9.2),
                min_value=1.0, step=0.1)
            st.session_state["piercap_depth_min"] = st.number_input(
                "Depth Min (m)", value=st.session_state.get("piercap_depth_min", 0.75),
                min_value=0.2, step=0.1)

        st.session_state["bearing_offset_long"] = st.number_input(
            "Bearing Offset Longitudinal (m)", value=st.session_state.get("bearing_offset_long", 0.825),
            min_value=0.0, step=0.05)

        st.markdown("#### Bearings")
        st.session_state["bearing_type"] = st.selectbox(
            "Bearing Type", ["Elastomeric", "Pot"],
            index=0 if st.session_state.get("bearing_type") == "Elastomeric" else 1)

        st.markdown("**LHS Bearing Coordinates** [along_traffic, across_traffic]")
        lhs = st.session_state.get("bearings_lhs", [[-0.825, -3.7], [-0.825, -2.3], [-0.825, 2.3], [-0.825, 3.7]])
        new_lhs = []
        for i, coord in enumerate(lhs):
            c1, c2 = st.columns(2)
            with c1:
                x = st.number_input(f"LHS Brg {i+1}: x", value=float(coord[0]),
                                    format="%.3f", step=0.1, key=f"lhs_x_{i}")
            with c2:
                y = st.number_input(f"LHS Brg {i+1}: y", value=float(coord[1]),
                                    format="%.3f", step=0.1, key=f"lhs_y_{i}")
            new_lhs.append([x, y])
        st.session_state["bearings_lhs"] = new_lhs

        st.markdown("**RHS Bearing Coordinates** [along_traffic, across_traffic]")
        rhs = st.session_state.get("bearings_rhs", [[0.825, -3.7], [0.825, -2.3], [0.825, 2.3], [0.825, 3.7]])
        new_rhs = []
        for i, coord in enumerate(rhs):
            c1, c2 = st.columns(2)
            with c1:
                x = st.number_input(f"RHS Brg {i+1}: x", value=float(coord[0]),
                                    format="%.3f", step=0.1, key=f"rhs_x_{i}")
            with c2:
                y = st.number_input(f"RHS Brg {i+1}: y", value=float(coord[1]),
                                    format="%.3f", step=0.1, key=f"rhs_y_{i}")
            new_rhs.append([x, y])
        st.session_state["bearings_rhs"] = new_rhs

        with st.expander("Pedestal Dimensions"):
            c1, c2, c3 = st.columns(3)
            with c1:
                st.session_state["pedestal_length"] = st.number_input(
                    "Pedestal Length (m)", value=st.session_state.get("pedestal_length", 0.8), step=0.1)
            with c2:
                st.session_state["pedestal_width"] = st.number_input(
                    "Pedestal Width (m)", value=st.session_state.get("pedestal_width", 0.8), step=0.1)
            with c3:
                st.session_state["pedestal_beyond"] = st.number_input(
                    "Pedestal Beyond (m)", value=st.session_state.get("pedestal_beyond", 0.2), step=0.05)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Pier Cap** is the transverse beam on top of the pier that "
            "distributes bearing loads.\n\n"
            "**Bearing coordinates** are relative to pier center:\n"
            "- x = along traffic direction\n"
            "- y = across traffic direction\n\n"
            "Tapered cap: depth varies linearly from max (at pier) to min (at tip)."
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons()


# ── Step 5: Foundation ───────────────────────────────────────────────────────

def _step_foundation() -> None:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("#### Pile Configuration")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.session_state["pile_diameter"] = st.number_input(
                "Pile Diameter (m)", value=st.session_state.get("pile_diameter", 1.0),
                min_value=0.3, step=0.1)
            st.session_state["piles_long"] = st.number_input(
                "Piles Longitudinal", value=st.session_state.get("piles_long", 2),
                min_value=1, max_value=10, step=1)
        with c2:
            st.session_state["pile_spacing_factor"] = st.number_input(
                "Spacing Factor (x dia)", value=st.session_state.get("pile_spacing_factor", 3),
                min_value=2, max_value=6, step=1)
            st.session_state["piles_trans"] = st.number_input(
                "Piles Transverse", value=st.session_state.get("piles_trans", 2),
                min_value=1, max_value=10, step=1)
        with c3:
            st.session_state["total_pile_depth"] = st.number_input(
                "Total Pile Depth (m)", value=st.session_state.get("total_pile_depth", 15.0),
                min_value=3.0, step=1.0)

        # Show pile grid preview
        n_long = st.session_state.get("piles_long", 2)
        n_trans = st.session_state.get("piles_trans", 2)
        total = n_long * n_trans
        st.caption(f"Pile grid: {n_long} x {n_trans} = {total} piles")

        st.markdown("#### Pile Cap")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.session_state["pilecap_thickness"] = st.number_input(
                "Pile Cap Thickness (m)", value=st.session_state.get("pilecap_thickness", 1.8),
                min_value=0.5, step=0.1)
        with c2:
            st.session_state["pilecap_edge_clearance"] = st.number_input(
                "Edge Clearance (mm)", value=st.session_state.get("pilecap_edge_clearance", 150),
                min_value=50, step=25)
        with c3:
            st.session_state["pilecap_top_below_gl"] = st.number_input(
                "Cap Top Below GL (mm)", value=st.session_state.get("pilecap_top_below_gl", 500),
                min_value=0, step=100)

        st.markdown("#### Pile Capacity & Soil")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.session_state["pile_capacity_top"] = st.number_input(
                "Safe Capacity (tonnes)", value=st.session_state.get("pile_capacity_top", 550.0),
                min_value=50.0, step=25.0)
            st.session_state["soil_density"] = st.number_input(
                "Soil Density (kN/m\u00b3)", value=st.session_state.get("soil_density", 20.0), step=1.0)
        with c2:
            st.session_state["pile_uplift_capacity"] = st.number_input(
                "Uplift Capacity (tonnes)", value=st.session_state.get("pile_uplift_capacity", -187.0),
                max_value=0.0, step=10.0)
            st.session_state["friction_angle"] = st.number_input(
                "Friction Angle (\u00b0)", value=st.session_state.get("friction_angle", 30.0),
                min_value=10.0, max_value=45.0, step=1.0)
        with c3:
            st.session_state["subgrade_modulus"] = st.number_input(
                "Subgrade Modulus (MN/m\u00b3)", value=st.session_state.get("subgrade_modulus", 1.4),
                min_value=0.1, step=0.1, format="%.1f")
            st.session_state["fixity_ratio"] = st.number_input(
                "Fixity Ratio (Lf/T)", value=st.session_state.get("fixity_ratio", 2.3),
                min_value=1.0, max_value=5.0, step=0.1, format="%.1f")

        st.markdown("#### Pile Reinforcement")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.session_state["pile_n_bars"] = st.number_input(
                "Pile Bars (nos)", value=st.session_state.get("pile_n_bars", 18),
                min_value=6, step=2)
        with c2:
            st.session_state["pile_bar_dia"] = st.number_input(
                "Bar Diameter (mm)", value=st.session_state.get("pile_bar_dia", 22.6),
                min_value=10.0, step=0.1, format="%.1f")
        with c3:
            st.session_state["pile_cover"] = st.number_input(
                "Cover (mm)", value=st.session_state.get("pile_cover", 75),
                min_value=40, max_value=100, step=5)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Foundation** defines the pile group and pile cap.\n\n"
            "**Spacing factor** = c/c spacing / pile diameter\n"
            "(typically 3\u20135)\n\n"
            "**Fixity ratio** Lf/T from IS:2911 Fig 4.\n\n"
            "**Safe capacity** includes FOS per IS:2911. "
            "IRC 78 adds 25% bonus for wind/seismic combinations."
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons()


# ── Step 6: Materials & Creep ────────────────────────────────────────────────

def _step_materials() -> None:
    col1, col2 = st.columns([2, 1])

    fck_options = [20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90]
    agg_options = ["Quartzite", "Granite", "Limestone", "Sandstone", "Basalt"]

    with col1:
        st.markdown("#### Concrete Grades (fck, MPa)")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.session_state["fck_pier"] = st.selectbox(
                "Pier", fck_options,
                index=fck_options.index(st.session_state.get("fck_pier", 50)))
            st.session_state["aggregate_pier"] = st.selectbox(
                "Pier Aggregate", agg_options,
                index=agg_options.index(st.session_state.get("aggregate_pier", "Quartzite")))
        with c2:
            st.session_state["fck_piercap"] = st.selectbox(
                "Pier Cap", fck_options,
                index=fck_options.index(st.session_state.get("fck_piercap", 45)))
            st.session_state["aggregate_piercap"] = st.selectbox(
                "Pier Cap Aggregate", agg_options,
                index=agg_options.index(st.session_state.get("aggregate_piercap", "Quartzite")))
        with c3:
            st.session_state["fck_pedestal"] = st.selectbox(
                "Pedestal", fck_options,
                index=fck_options.index(st.session_state.get("fck_pedestal", 45)),
                key="fck_ped_sel")
            st.session_state["fck_pile_pilecap"] = st.selectbox(
                "Pile / Pile Cap", fck_options,
                index=fck_options.index(st.session_state.get("fck_pile_pilecap", 35)))

        st.markdown("#### Steel")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["fyk"] = st.number_input(
                "fyk (MPa)", value=st.session_state.get("fyk", 550),
                min_value=400, max_value=600, step=50)
        with c2:
            st.session_state["Es"] = st.number_input(
                "Es (MPa)", value=st.session_state.get("Es", 200000),
                min_value=190000, max_value=210000, step=5000)

        st.markdown("#### General")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["exposure"] = st.selectbox(
                "Exposure Condition",
                ["Mild", "Moderate", "Severe", "Very Severe", "Extreme"],
                index=["Mild", "Moderate", "Severe", "Very Severe", "Extreme"].index(
                    st.session_state.get("exposure", "Severe")))
            st.session_state["concrete_density"] = st.number_input(
                "Concrete Density (kN/m\u00b3)", value=st.session_state.get("concrete_density", 25), step=1)
        with c2:
            st.session_state["water_density"] = st.number_input(
                "Water Density (kN/m\u00b3)", value=st.session_state.get("water_density", 10), step=1)

        st.markdown("#### Creep Parameters (IRC 112 Annex A-2)")
        c1, c2, c3 = st.columns(3)
        with c1:
            st.session_state["age_of_loading"] = st.number_input(
                "Age of Loading (days)", value=st.session_state.get("age_of_loading", 7),
                min_value=1, step=1)
        with c2:
            st.session_state["relative_humidity"] = st.number_input(
                "Relative Humidity (%)", value=st.session_state.get("relative_humidity", 54.5),
                min_value=20.0, max_value=100.0, step=5.0)
        with c3:
            st.session_state["design_life"] = st.number_input(
                "Design Life (days)", value=st.session_state.get("design_life", 36500),
                min_value=3650, step=3650)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Materials** can differ by component:\n"
            "- Pier: typically M45\u2013M50\n"
            "- Pier cap: M40\u2013M45\n"
            "- Pile/pile cap: M30\u2013M40\n\n"
            "**Exposure** affects clear cover per IRC 112 Table 14.1 "
            "and crack width limits."
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons()


# ── Step 7: Seismic & Wind ───────────────────────────────────────────────────

def _step_seismic_wind() -> None:
    col1, col2 = st.columns([2, 1])

    with col1:
        st.markdown("#### Seismic (IRC 6 Cl. 211)")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["seismic_zone"] = st.selectbox(
                "Seismic Zone",
                ["II", "III", "IV", "V"],
                index=["II", "III", "IV", "V"].index(
                    st.session_state.get("seismic_zone", "IV")))
            st.session_state["importance_factor"] = st.number_input(
                "Importance Factor (I)", value=st.session_state.get("importance_factor", 1.2),
                min_value=1.0, max_value=1.5, step=0.1)
        with c2:
            st.session_state["response_reduction"] = st.number_input(
                "Response Reduction (R)", value=st.session_state.get("response_reduction", 4.0),
                min_value=1.0, max_value=5.0, step=0.5)
            st.session_state["soil_type"] = st.selectbox(
                "Soil Type",
                ["hard", "medium", "soft"],
                index=["hard", "medium", "soft"].index(
                    st.session_state.get("soil_type", "medium")))

        st.markdown("#### Wind (IRC 6 Cl. 209)")
        c1, c2 = st.columns(2)
        with c1:
            st.session_state["wind_speed"] = st.number_input(
                "Basic Wind Speed (m/s)", value=st.session_state.get("wind_speed", 44.0),
                min_value=20.0, max_value=60.0, step=1.0)
        with c2:
            st.session_state["terrain_category"] = st.number_input(
                "Terrain Category (1\u20134)", value=st.session_state.get("terrain_category", 2),
                min_value=1, max_value=4, step=1)

    with col2:
        st.markdown('<div class="info-box">', unsafe_allow_html=True)
        st.markdown(
            "**Seismic zones (IRC 6):**\n"
            "- Zone II: Z=0.10 (Low)\n"
            "- Zone III: Z=0.16 (Moderate)\n"
            "- Zone IV: Z=0.24 (Severe)\n"
            "- Zone V: Z=0.36 (Very Severe)\n\n"
            "**R factor:**\n"
            "- Ductile: R=4.0\n"
            "- Limited ductile: R=3.0\n"
            "- Brittle: R=1.0"
        )
        st.markdown('</div>', unsafe_allow_html=True)

    st.markdown("---")
    _nav_buttons()


# ── Step 8: Review & Run ─────────────────────────────────────────────────────

def _step_review() -> None:
    st.markdown("### Review Your Inputs")

    s = st.session_state

    # Summary tables
    st.markdown("#### Project")
    st.markdown(
        f"| Parameter | Value |\n|---|---|\n"
        f"| Project | {s.get('project_name')} |\n"
        f"| Pier ID | {s.get('pier_id')} |\n"
        f"| FRL | {s.get('frl')} m |\n"
        f"| GL | {s.get('gl')} m |"
    )

    st.markdown("#### Superstructure")
    st.markdown(
        f"| Parameter | Value |\n|---|---|\n"
        f"| Deck Width | {s.get('deck_width')} m |\n"
        f"| Type | {s.get('super_type')} |\n"
        f"| Left Span | {s.get('left_span')} m |\n"
        f"| Right Span | {s.get('right_span')} m |\n"
        f"| Depth | {s.get('depth_incl_slab')} m |\n"
        f"| Continuity | {s.get('continuity')} |"
    )

    st.markdown("#### Pier")
    st.markdown(
        f"| Parameter | Value |\n|---|---|\n"
        f"| Type | {s.get('pier_type')} |\n"
        f"| Diameter / Width | {s.get('pier_diameter_bottom')} m |\n"
        f"| End Condition | {s.get('pier_end_condition')} |\n"
        f"| Reinforcement | {s.get('pier_n_bars')} nos x {s.get('pier_bar_dia')} mm |"
    )

    st.markdown("#### Foundation")
    n_piles = s.get("piles_long", 2) * s.get("piles_trans", 2)
    st.markdown(
        f"| Parameter | Value |\n|---|---|\n"
        f"| Pile Diameter | {s.get('pile_diameter')} m |\n"
        f"| Pile Grid | {s.get('piles_long')} x {s.get('piles_trans')} = {n_piles} piles |\n"
        f"| Pile Depth | {s.get('total_pile_depth')} m |\n"
        f"| Safe Capacity | {s.get('pile_capacity_top')} tonnes |\n"
        f"| Pile Cap Thickness | {s.get('pilecap_thickness')} m |"
    )

    st.markdown("#### Materials")
    st.markdown(
        f"| Component | fck (MPa) |\n|---|---|\n"
        f"| Pier | {s.get('fck_pier')} |\n"
        f"| Pier Cap | {s.get('fck_piercap')} |\n"
        f"| Pile / Pile Cap | {s.get('fck_pile_pilecap')} |\n"
        f"| Steel fyk | {s.get('fyk')} MPa |\n"
        f"| Exposure | {s.get('exposure')} |"
    )

    st.markdown("#### Loading")
    st.markdown(
        f"| Parameter | Value |\n|---|---|\n"
        f"| Seismic Zone | {s.get('seismic_zone')} |\n"
        f"| Importance Factor | {s.get('importance_factor')} |\n"
        f"| R Factor | {s.get('response_reduction')} |\n"
        f"| Wind Speed | {s.get('wind_speed')} m/s |"
    )

    st.markdown("---")

    # BOQ rates
    c1, c2 = st.columns(2)
    with c1:
        st.session_state["concrete_rate"] = st.number_input(
            "Concrete Rate (INR/m\u00b3)", value=s.get("concrete_rate", 8000), step=500)
    with c2:
        st.session_state["steel_rate"] = st.number_input(
            "Steel Rate (INR/tonne)", value=s.get("steel_rate", 70000), step=5000)

    st.markdown("---")

    # Navigation
    col1, col2 = st.columns(2)
    with col1:
        if st.button("\u2190 Back", use_container_width=True):
            st.session_state["current_step"] = 7
            st.rerun()
    with col2:
        if st.button("Run Design", type="primary", use_container_width=True):
            config = build_config_from_state()
            with st.spinner("Running substructure design..."):
                results = run_design(config)
            st.session_state["results"] = results
            st.session_state["config"] = config
            st.rerun()

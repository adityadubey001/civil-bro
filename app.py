"""
IRC 112 Bridge Girder Design Application

Streamlined design for simply supported rectangular girders
per IRC 112:2020 with IRC 6 vehicle loading.
"""

import streamlit as st
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from src.models.inputs import (
    BridgeGirderInput, DeckInput,
    ConcreteGrade, SteelGrade, ExposureCondition,
    ReinforcementPreference, LiveLoadType, GirderType
)
from src.core.beam_design import BridgeGirderDesignEngine
from src.models.outputs import DesignStatus

st.set_page_config(
    page_title="IRC 112 Bridge Girder Designer",
    page_icon="Bridge",
    layout="wide",
    initial_sidebar_state="collapsed"
)

st.markdown("""
<style>
    .info-box {
        background-color: #e8f4f8;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3498db;
    }
    .result-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
        border: 1px solid #dee2e6;
    }
</style>
""", unsafe_allow_html=True)


def init_state():
    """Initialize session state with defaults."""
    defaults = {
        'step': 1,
        'span_length': 15.0,
        'deck_width': 2.75,
        'result': None,
        'deck_thickness': 240,
        'wearing_coat': 75,
        'concrete_grade': "M40",
        'steel_grade': "Fe500",
    }
    for key, val in defaults.items():
        if key not in st.session_state:
            st.session_state[key] = val


def run_design():
    """Execute the bridge girder design."""
    deck = DeckInput(
        width=st.session_state.deck_width,
        thickness=st.session_state.deck_thickness,
        wearing_coat=st.session_state.wearing_coat,
    )

    inputs = BridgeGirderInput(
        span_length=st.session_state.span_length,
        deck=deck,
        num_girders=2,
        girder_type=GirderType.RECTANGULAR,  # Rectangular - deck only for loading
        concrete_grade=ConcreteGrade(st.session_state.concrete_grade),
        steel_grade=SteelGrade(st.session_state.steel_grade),
        live_load_type=LiveLoadType.CLASS_70R_WHEELED,  # All cases analyzed, this is just default
        exposure_condition=ExposureCondition.MODERATE,
        reinforcement_preference=ReinforcementPreference.ALLOW_DOUBLY,
        preferred_bar_sizes=[20, 25, 32],
        preferred_stirrup_size=12,
    )

    engine = BridgeGirderDesignEngine()
    st.session_state.result = engine.design(inputs)
    st.session_state.step = 3


def render_step1():
    """Step 1: Span Length"""
    st.markdown("## Step 1: What is the span of your bridge?")

    col1, col2 = st.columns([2, 1])

    with col1:
        span = st.number_input(
            "Effective Span (meters)",
            min_value=5.0,
            max_value=30.0,
            value=st.session_state.span_length,
            step=0.5,
            help="Clear span between supports"
        )
        st.session_state.span_length = span

        if span <= 10:
            st.info("Short span - Rectangular girders are economical")
        elif span <= 20:
            st.info("Medium span - Standard design applies")
        else:
            st.warning("Long span - Consider deeper sections")

    with col2:
        st.code(f"""
    <-- {span:.1f} m -->

    ================
         |  |
         |  | Girders
    ----+--+----
    Support  Support
        """)

    st.markdown("---")
    if st.button("Next ->", type="primary", use_container_width=True):
        st.session_state.step = 2
        st.rerun()


def render_step2():
    """Step 2: Deck Width and Run Design"""
    st.markdown("## Step 2: What is the deck width per girder?")

    col1, col2 = st.columns([2, 1])

    with col1:
        deck_width = st.number_input(
            "Deck width per girder (meters)",
            min_value=1.0,
            max_value=6.0,
            value=st.session_state.deck_width,
            step=0.25,
            help="Tributary width of deck slab supported by each girder"
        )
        st.session_state.deck_width = deck_width

        st.markdown(f"""
        <div class="info-box">
        This determines the <b>dead load</b> from deck slab on each girder.<br>
        Typical values: 2.5m - 4.0m for highway bridges
        </div>
        """, unsafe_allow_html=True)

    with col2:
        st.markdown("**Example:**")
        st.markdown("""
        - Total carriage width: 7.5m
        - Number of girders: 3
        - Width per girder: 2.5m
        """)

    # Advanced options (collapsed)
    with st.expander("Advanced Options"):
        col1, col2 = st.columns(2)
        with col1:
            st.session_state.deck_thickness = st.number_input(
                "Deck thickness (mm)", 150, 400, st.session_state.deck_thickness, 10
            )
            st.session_state.wearing_coat = st.number_input(
                "Wearing coat (mm)", 40, 100, st.session_state.wearing_coat, 5
            )
        with col2:
            st.session_state.concrete_grade = st.selectbox(
                "Concrete grade",
                ["M30", "M35", "M40", "M45", "M50"],
                index=["M30", "M35", "M40", "M45", "M50"].index(st.session_state.concrete_grade)
            )
            st.session_state.steel_grade = st.selectbox(
                "Steel grade",
                ["Fe415", "Fe500", "Fe550"],
                index=["Fe415", "Fe500", "Fe550"].index(st.session_state.steel_grade)
            )

    # Summary
    st.markdown("### Design Parameters")
    st.markdown(f"""
    | Parameter | Value |
    |-----------|-------|
    | Span | {st.session_state.span_length} m |
    | Deck width/girder | {st.session_state.deck_width} m |
    | Girder type | Rectangular |
    | Live load | IRC 6 (all cases analyzed) |
    | Concrete | {st.session_state.concrete_grade} |
    | Steel | {st.session_state.steel_grade} |
    """)

    st.markdown("---")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("<- Back", use_container_width=True):
            st.session_state.step = 1
            st.rerun()
    with col2:
        if st.button("Run Design", type="primary", use_container_width=True):
            with st.spinner("Designing girder per IRC 112..."):
                try:
                    run_design()
                    st.rerun()
                except Exception as e:
                    st.error(f"Design failed: {e}")
                    import traceback
                    st.code(traceback.format_exc())


def render_results():
    """Step 3: Results"""
    result = st.session_state.result

    if result is None:
        st.warning("No results. Please complete the design steps.")
        if st.button("Start Over"):
            st.session_state.step = 1
            st.rerun()
        return

    # Status banner
    if result.overall_status == DesignStatus.PASS:
        st.success("DESIGN ADEQUATE - All ULS and SLS checks passed")
    elif result.overall_status == DesignStatus.WARNING:
        st.warning("DESIGN OK WITH WARNINGS")
    else:
        st.error("DESIGN INADEQUATE")

    st.markdown("## Design Results")

    # Key results
    col1, col2, col3 = st.columns(3)

    with col1:
        st.markdown("### Section")
        st.metric("Width", f"{result.girder_width:.0f} mm")
        st.metric("Depth", f"{result.girder_depth:.0f} mm")
        st.metric("Effective Depth", f"{result.effective_depth:.0f} mm")

    with col2:
        st.markdown("### Reinforcement")
        st.metric("Bottom Steel", result.flexure.tension_bar_arrangement)
        st.metric("Area", f"{result.flexure.provided_ast:.0f} mm²")
        st.metric("Stirrups", f"{result.shear.stirrup_legs}L-{result.shear.stirrup_diameter}φ @ {result.shear.spacing_provided:.0f}mm")

    with col3:
        st.markdown("### Design Forces")
        st.metric("Mu (factored)", f"{result.loads.factored_moment:.0f} kNm")
        st.metric("Vu (factored)", f"{result.loads.factored_shear:.0f} kN")
        st.metric("Governing LL", result.loads.governing_bm_case)

    # ULS Checks Section
    st.markdown("### ULS Checks (IRC 112)")

    # Flexure Check
    st.markdown("#### Flexure Check (Cl. 10.3)")
    col1, col2, col3 = st.columns(3)

    with col1:
        flex_status = "PASS" if result.flexure.status.value == "pass" else "FAIL"
        flex_util = result.flexure.utilization_ratio * 100
        st.metric(
            "Moment Capacity (Mu ≤ MuR)",
            f"{result.loads.factored_moment:.0f} / {result.flexure.moment_capacity:.0f} kNm",
            delta=f"{flex_status} ({flex_util:.0f}%)"
        )

    with col2:
        xu_status = "PASS" if result.flexure.is_under_reinforced else "FAIL"
        st.metric(
            "Ductility (xu/d ≤ xu,max/d)",
            f"{result.flexure.xu_d_ratio:.3f} / 0.456",
            delta=f"{xu_status} ({'Under' if result.flexure.is_under_reinforced else 'Over'}-reinforced)"
        )

    with col3:
        min_status = "PASS" if result.flexure.min_ast_satisfied else "FAIL"
        st.metric(
            "Min Steel Check",
            f"{result.flexure.pt_provided:.2f}%",
            delta=f"{min_status}"
        )

    # Shear Check
    st.markdown("#### Shear Check (Cl. 10.4)")
    col1, col2, col3 = st.columns(3)

    with col1:
        shear_status = "PASS" if result.shear.status.value == "pass" else "FAIL"
        shear_util = (result.shear.shear_at_d_from_support / result.shear.max_shear_capacity) * 100
        st.metric(
            "Shear Capacity (Vu ≤ VRd,max)",
            f"{result.shear.shear_at_d_from_support:.0f} / {result.shear.max_shear_capacity:.0f} kN",
            delta=f"{shear_status} ({shear_util:.0f}%)"
        )

    with col2:
        vrdc_status = "Yes" if result.shear.shear_reinforcement_required else "No"
        st.metric(
            "Concrete Shear (VRd,c)",
            f"{result.shear.concrete_shear_capacity:.0f} kN",
            delta=f"Reinf required: {vrdc_status}"
        )

    with col3:
        st.metric(
            "Stirrup Spacing",
            f"{result.shear.spacing_provided:.0f} / {result.shear.spacing_max:.0f} mm",
            delta="OK" if result.shear.spacing_provided <= result.shear.spacing_max else "FAIL"
        )

    # SLS Checks Section
    st.markdown("### SLS Checks (IRC 112)")

    # Rare Combination Stress Check
    st.markdown("#### Stress Check - Rare Combination (Cl. 12.2.1)")
    col1, col2 = st.columns(2)

    with col1:
        c_status = "PASS" if result.sls_stress.concrete_stress_check else "FAIL"
        c_util = result.sls_stress.concrete_utilization * 100
        st.metric(
            "Concrete Stress (σc ≤ 0.48fck)",
            f"{result.sls_stress.concrete_stress:.1f} / {result.sls_stress.concrete_stress_limit:.1f} MPa",
            delta=f"{c_status} ({c_util:.0f}%)"
        )

    with col2:
        s_status = "PASS" if result.sls_stress.steel_stress_check else "FAIL"
        s_util = result.sls_stress.steel_utilization * 100
        st.metric(
            "Steel Stress (σs ≤ 0.8fy)",
            f"{result.sls_stress.steel_stress:.0f} / {result.sls_stress.steel_stress_limit:.0f} MPa",
            delta=f"{s_status} ({s_util:.0f}%)"
        )

    # Quasi-permanent Stress Check
    st.markdown("#### Stress Check - Quasi-Permanent (Cl. 12.2.1(3))")
    col1, col2 = st.columns(2)

    with col1:
        if result.sls_stress.qp_concrete_stress is not None:
            qp_status = "PASS" if result.sls_stress.qp_concrete_stress_check else "FAIL"
            qp_util = (result.sls_stress.qp_concrete_stress / result.sls_stress.qp_concrete_stress_limit) * 100
            st.metric(
                "Concrete Stress (σc,qp ≤ 0.36fck)",
                f"{result.sls_stress.qp_concrete_stress:.1f} / {result.sls_stress.qp_concrete_stress_limit:.1f} MPa",
                delta=f"{qp_status} ({qp_util:.0f}%)"
            )
        else:
            st.info("QP check not available")

    with col2:
        if result.sls_stress.qp_moment is not None:
            st.metric(
                "Quasi-Permanent Moment",
                f"{result.sls_stress.qp_moment:.0f} kNm",
                delta="Mdl only (ψ2=0)"
            )

    # Crack Width Check
    st.markdown("#### Crack Width Check (Cl. 12.3.4)")
    if result.crack_width is not None:
        col1, col2, col3 = st.columns(3)

        with col1:
            cw_status = "PASS" if result.crack_width.crack_width_check else "FAIL"
            st.metric(
                "Crack Width (wk)",
                f"{result.crack_width.calculated_crack_width:.3f} / {result.crack_width.allowable_crack_width:.1f} mm",
                delta=cw_status
            )

        with col2:
            st.metric(
                "Max Crack Spacing (sr,max)",
                f"{result.crack_width.max_crack_spacing:.0f} mm"
            )

        with col3:
            st.metric(
                "Steel Stress (QP)",
                f"{result.crack_width.steel_stress_qp:.0f} MPa"
            )
    else:
        st.info("Crack width check not available")

    # Detailed tabs
    tab1, tab2, tab3, tab4 = st.tabs(["Load Analysis", "Flexure", "Shear", "SLS Details"])

    with tab1:
        st.markdown("#### Dead Loads")
        col1, col2, col3 = st.columns(3)
        col1.metric("Girder SW", f"{result.loads.sw_girder:.2f} kN/m")
        col2.metric("Deck Slab", f"{result.loads.sw_deck:.2f} kN/m")
        col3.metric("Wearing Coat", f"{result.loads.sw_wearing_coat:.2f} kN/m")
        st.metric("Total DL", f"{result.loads.total_dead_load:.2f} kN/m")

        st.markdown("#### Live Load Analysis (IRC 6) - All Cases")
        for ll in result.loads.live_load_results:
            gov = " **(GOVERNS)**" if ll.is_governing_bm else ""
            st.markdown(f"**{ll.load_type}**{gov}: BM = {ll.max_bm_with_impact:.0f} kNm, SF = {ll.max_sf_with_impact:.0f} kN (IF = {ll.impact_factor:.3f})")

        st.markdown("#### Load Combination (ULS: 1.35 DL + 1.50 LL)")
        st.markdown(f"**Mu** = 1.35 × {result.loads.dead_load_moment:.0f} + 1.50 × {result.loads.governing_live_moment:.0f} = **{result.loads.factored_moment:.0f} kNm**")
        st.markdown(f"**Vu** = 1.35 × {result.loads.dead_load_shear:.0f} + 1.50 × {result.loads.governing_live_shear:.0f} = **{result.loads.factored_shear:.0f} kN**")

    with tab2:
        st.markdown(f"**Design Moment:** {result.flexure.design_moment:.0f} kNm")
        st.markdown(f"**Moment Capacity:** {result.flexure.moment_capacity:.0f} kNm")
        st.markdown(f"**Utilization:** {result.flexure.utilization_ratio*100:.1f}%")
        st.markdown(f"**xu/d ratio:** {result.flexure.xu_d_ratio:.3f} (limit: 0.456 for Fe500)")
        st.markdown(f"**Under-reinforced:** {'Yes' if result.flexure.is_under_reinforced else 'No'}")
        if result.flexure.is_doubly_reinforced:
            st.markdown(f"**Compression Steel:** {result.flexure.compression_bar_arrangement}")

    with tab3:
        st.markdown(f"**Design Shear:** {result.shear.design_shear:.0f} kN")
        st.markdown(f"**Shear at d from support:** {result.shear.shear_at_d_from_support:.0f} kN")
        st.markdown(f"**VRd,c (concrete):** {result.shear.concrete_shear_capacity:.0f} kN")
        st.markdown(f"**VRd,max:** {result.shear.max_shear_capacity:.0f} kN")
        st.markdown(f"**Size Effect k:** {result.shear.size_effect_factor_k:.3f}")
        st.markdown(f"**Shear reinf required:** {'Yes' if result.shear.shear_reinforcement_required else 'No (minimum provided)'}")

    with tab4:
        st.markdown("#### SLS Stress Analysis Details")
        st.markdown(f"**Modular Ratio (m):** {result.sls_stress.modular_ratio:.2f}")
        st.markdown(f"**Elastic NA Depth (x):** {result.sls_stress.elastic_na_depth:.0f} mm")
        st.markdown(f"**Cracked MI (Icr):** {result.sls_stress.cracked_moment_of_inertia/1e9:.4f} × 10⁹ mm⁴")

        st.markdown("---")
        st.markdown("##### Rare Combination (Mdl + Mll)")
        st.markdown(f"**Service Moment:** {result.sls_stress.rare_combination_moment:.0f} kNm")
        st.markdown(f"**Concrete Stress:** {result.sls_stress.concrete_stress:.2f} MPa ≤ {result.sls_stress.concrete_stress_limit:.2f} MPa (0.48×fck)")
        st.markdown(f"**Steel Stress:** {result.sls_stress.steel_stress:.0f} MPa ≤ {result.sls_stress.steel_stress_limit:.0f} MPa (0.8×fy)")

        if result.sls_stress.qp_moment is not None:
            st.markdown("---")
            st.markdown("##### Quasi-Permanent Combination (Mdl only, ψ2=0 per Table B.3)")
            st.markdown(f"**QP Moment:** {result.sls_stress.qp_moment:.0f} kNm (dead load only)")
            st.markdown(f"**Concrete Stress (QP):** {result.sls_stress.qp_concrete_stress:.2f} MPa ≤ {result.sls_stress.qp_concrete_stress_limit:.2f} MPa (0.36×fck)")
            st.markdown("*Note: ψ2=0 for traffic per IRC 112 Table B.3. Limit 0.36×fck to avoid non-linear creep per Cl. 12.2.1(3)*")

        if result.crack_width is not None:
            st.markdown("---")
            st.markdown("##### Crack Width Calculation (Cl. 12.3.4)")
            st.markdown(f"**Effective Rebar Ratio (ρp,eff):** {result.crack_width.effective_reinforcement_ratio:.4f}")
            st.markdown(f"**Max Crack Spacing (sr,max):** {result.crack_width.max_crack_spacing:.0f} mm")
            st.markdown(f"**Strain Difference (εsm - εcm):** {result.crack_width.strain_difference*1000:.4f} × 10⁻³")
            st.markdown(f"**Calculated Crack Width:** {result.crack_width.calculated_crack_width:.3f} mm")
            st.markdown(f"**Allowable Crack Width:** {result.crack_width.allowable_crack_width:.1f} mm")

        st.markdown("---")
        st.markdown("##### Deflection Check")
        st.markdown(f"**Calculated:** {result.deflection.calculated_deflection:.1f} mm")
        st.markdown(f"**Allowable (L/250):** {result.deflection.allowable_deflection_total:.1f} mm")
        st.markdown(f"**Status:** {'PASS' if result.deflection.deflection_check_passed else 'FAIL'}")

    # Warnings
    if result.warnings:
        st.markdown("### Warnings")
        for w in result.warnings:
            st.warning(w)

    st.markdown("---")
    if st.button("Start New Design", use_container_width=True):
        st.session_state.step = 1
        st.session_state.result = None
        st.rerun()


def main():
    init_state()

    st.markdown("# IRC 112 Bridge Girder Designer")
    st.markdown("*Rectangular girder design with IRC 6 vehicle loading (auto-analyzed)*")
    st.markdown("---")

    if st.session_state.step == 1:
        render_step1()
    elif st.session_state.step == 2:
        render_step2()
    else:
        render_results()


if __name__ == "__main__":
    main()

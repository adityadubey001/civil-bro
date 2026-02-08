"""Results display for the Substructure Design app.

Renders the design dashboard and detail tabs after a design run completes.
"""

from __future__ import annotations

import streamlit as st

from ui.runner import safe_attr, get_status, all_checks_ok


# ── Formatting helpers ───────────────────────────────────────────────────────

def _fnum(val, dp: int = 2, unit: str = "") -> str:
    """Format a number safely."""
    if val is None:
        return "—"
    try:
        s = f"{float(val):.{dp}f}"
        return f"{s} {unit}".strip() if unit else s
    except (TypeError, ValueError):
        return str(val)


def _status_delta(obj) -> str | None:
    """Return a delta string for st.metric based on status."""
    s = get_status(obj)
    if s == "OK":
        return "OK"
    elif s == "NOT OK":
        return "FAIL"
    return None


def _delta_color(obj) -> str:
    """Return 'normal' for OK (green up arrow) or 'inverse' for FAIL."""
    s = get_status(obj)
    if s == "OK":
        return "normal"
    return "inverse"


# ── Main render function ─────────────────────────────────────────────────────

def render_results(results: dict) -> None:
    """Render the full results page."""

    st.markdown("## Substructure Design Results")

    # ── Status banner ────────────────────────────────────────────────────
    if all_checks_ok(results):
        st.success("All design checks passed.")
    else:
        errors = results.get("errors", [])
        if errors:
            st.error("Design encountered errors:")
            for e in errors:
                st.markdown(f"- {e}")
        else:
            st.warning("Some design checks did not pass — review details below.")

    # Show warnings
    for w in results.get("warnings", []):
        st.warning(w)

    # ── Dashboard metrics ────────────────────────────────────────────────
    _render_dashboard(results)

    st.markdown("---")

    # ── Detail tabs ──────────────────────────────────────────────────────
    _render_detail_tabs(results)


# ── Dashboard ────────────────────────────────────────────────────────────────

def _render_dashboard(results: dict) -> None:
    """Render the 6-card summary dashboard."""

    pc = results.get("pier_cap")
    pier = results.get("pier")
    pile_cap = results.get("pile_capacity")
    pile_d = results.get("pile_design")
    pilecap = results.get("pilecap")
    boq = results.get("boq")

    cols = st.columns(6)

    with cols[0]:
        _design_metric(
            "Pier Cap",
            pc, "flexure_util_pier",
            sub_attrs=["shear_util"],
            sub_labels=["Shear"],
        )

    with cols[1]:
        _design_metric(
            "Pier",
            pier, "util_biaxial",
        )

    with cols[2]:
        _design_metric(
            "Pile Capacity",
            pile_cap, "compression_util",
        )

    with cols[3]:
        _design_metric(
            "Pile Design",
            pile_d, "util_PM",
        )

    with cols[4]:
        _design_metric(
            "Pile Cap",
            pilecap, "punch_util_pile",
            sub_attrs=["punch_util_pier"],
            sub_labels=["Pier punch"],
        )

    with cols[5]:
        if boq is not None:
            cost = safe_attr(boq, "total_cost", 0)
            cost_lakhs = cost / 1e5
            st.metric("Cost", f"INR {cost_lakhs:.1f}L")
        else:
            st.metric("Cost", "—")


def _design_metric(label: str, obj, primary_attr: str,
                   sub_attrs: list[str] | None = None,
                   sub_labels: list[str] | None = None) -> None:
    """Render a design check metric card."""
    val = safe_attr(obj, primary_attr)
    status = get_status(obj)

    if val is not None:
        st.metric(
            label,
            _fnum(val, 3),
            delta=status,
            delta_color="normal" if status == "OK" else "inverse",
        )
    else:
        st.metric(label, "—", delta="Not computed")

    # Show sub-attributes as small text
    if sub_attrs and obj is not None:
        parts = []
        for attr, lbl in zip(sub_attrs, sub_labels or sub_attrs):
            v = safe_attr(obj, attr)
            if v is not None:
                parts.append(f"{lbl}: {_fnum(v, 3)}")
        if parts:
            st.caption(" | ".join(parts))


# ── Detail tabs ──────────────────────────────────────────────────────────────

def _render_detail_tabs(results: dict) -> None:
    """Render the 6 detail tabs."""

    tabs = st.tabs(["Pier Cap", "Pier", "Piles", "Pile Cap", "Loads", "BOQ"])

    with tabs[0]:
        _render_pier_cap_tab(results)
    with tabs[1]:
        _render_pier_tab(results)
    with tabs[2]:
        _render_piles_tab(results)
    with tabs[3]:
        _render_pilecap_tab(results)
    with tabs[4]:
        _render_loads_tab(results)
    with tabs[5]:
        _render_boq_tab(results)


# ── Pier Cap tab ─────────────────────────────────────────────────────────────

def _render_pier_cap_tab(results: dict) -> None:
    pc = results.get("pier_cap")
    if pc is None:
        st.info("Pier cap design was not computed.")
        return

    st.markdown("### Pier Cap Design")

    # Geometry
    st.markdown("#### Section Geometry")
    c1, c2, c3 = st.columns(3)
    c1.metric("Cantilever Length", _fnum(safe_attr(pc, "cantilever_length_mm"), 0, "mm"))
    c2.metric("Depth at Pier Face", _fnum(safe_attr(pc, "depth_pier_mm"), 0, "mm"))
    c3.metric("Depth at Curtailment", _fnum(safe_attr(pc, "depth_curt_mm"), 0, "mm"))

    # Flexure
    st.markdown("#### Flexure (IRC 112 Cl. 8.2)")
    c1, c2, c3 = st.columns(3)
    c1.metric("M_Ed (Pier Face)", _fnum(safe_attr(pc, "MEd_pier"), 1, "kN.m"))
    c2.metric("Utilization (Pier)", _fnum(safe_attr(pc, "flexure_util_pier"), 3),
              delta="OK" if safe_attr(pc, "flexure_util_pier", 999) <= 1.0 else "FAIL",
              delta_color="normal" if safe_attr(pc, "flexure_util_pier", 999) <= 1.0 else "inverse")
    c3.metric("Utilization (Curtailment)", _fnum(safe_attr(pc, "flexure_util_curt"), 3),
              delta="OK" if safe_attr(pc, "flexure_util_curt", 999) <= 1.0 else "FAIL",
              delta_color="normal" if safe_attr(pc, "flexure_util_curt", 999) <= 1.0 else "inverse")

    # Reinforcement
    c1, c2 = st.columns(2)
    c1.metric("Ast Provided (Pier)", _fnum(safe_attr(pc, "Ast_pier_provided"), 0, "mm\u00b2"))
    c2.metric("Ast Provided (Curt)", _fnum(safe_attr(pc, "Ast_curt_provided"), 0, "mm\u00b2"))

    # Shear
    st.markdown("#### Shear (IRC 112 Cl. 10.3)")
    c1, c2, c3 = st.columns(3)
    c1.metric("V_Ed", _fnum(safe_attr(pc, "VEd_pier"), 1, "kN"))
    c2.metric("V_Rd,s", _fnum(safe_attr(pc, "VRds"), 1, "kN"))
    c3.metric("Shear Utilization", _fnum(safe_attr(pc, "shear_util"), 3),
              delta="OK" if safe_attr(pc, "shear_util", 999) <= 1.0 else "FAIL",
              delta_color="normal" if safe_attr(pc, "shear_util", 999) <= 1.0 else "inverse")

    # Crack width
    st.markdown("#### Crack Width (IRC 112 Cl. 12.3.4)")
    c1, c2 = st.columns(2)
    c1.metric("Crack Width (w_k)", _fnum(safe_attr(pc, "crack_width"), 3, "mm"))
    c2.metric("Limit", _fnum(safe_attr(pc, "crack_width_limit"), 1, "mm"),
              delta="OK" if safe_attr(pc, "crack_width", 999) <= safe_attr(pc, "crack_width_limit", 0.3) else "FAIL",
              delta_color="normal" if safe_attr(pc, "crack_width", 999) <= safe_attr(pc, "crack_width_limit", 0.3) else "inverse")


# ── Pier tab ─────────────────────────────────────────────────────────────────

def _render_pier_tab(results: dict) -> None:
    pier = results.get("pier")
    if pier is None:
        st.info("Pier design was not computed.")
        return

    st.markdown("### Pier Design")

    # Section properties
    st.markdown("#### Section Properties")
    c1, c2, c3 = st.columns(3)
    c1.metric("Diameter / Width", _fnum(safe_attr(pier, "diameter_mm"), 0, "mm"))
    c2.metric("Reinforcement", f"{safe_attr(pier, 'n_bars', '—')} nos x {safe_attr(pier, 'bar_dia_mm', '—')} mm")
    c3.metric("Steel Ratio", _fnum(safe_attr(pier, "steel_pct"), 2, "%"))

    # Slenderness
    st.markdown("#### Slenderness Check (IRC 112 Cl. 8.3)")
    c1, c2, c3 = st.columns(3)
    c1.metric("Slenderness (\u03bb)", _fnum(safe_attr(pier, "lambda_"), 1))
    c2.metric("Limit (\u03bb_lim)", _fnum(safe_attr(pier, "lambda_lim"), 1))
    c3.metric("2nd Order Effects",
              "Required" if safe_attr(pier, "second_order_needed") else "Not required")

    # ULS Biaxial
    st.markdown("#### ULS Biaxial Check")
    c1, c2, c3 = st.columns(3)
    c1.metric("Biaxial Utilization", _fnum(safe_attr(pier, "util_biaxial"), 3),
              delta="OK" if safe_attr(pier, "util_biaxial", 999) <= 1.0 else "FAIL",
              delta_color="normal" if safe_attr(pier, "util_biaxial", 999) <= 1.0 else "inverse")
    c2.metric("Governing Combo", safe_attr(pier, "governing_combo_name", "—"))
    c3.metric("Axial Load (P)", _fnum(safe_attr(pier, "P_governing"), 1, "kN"))

    # SLS
    st.markdown("#### SLS Checks")
    c1, c2, c3 = st.columns(3)
    c1.metric("Concrete Stress (Rare)", _fnum(safe_attr(pier, "sigma_c_rare"), 2, "MPa"))
    c2.metric("Steel Stress (Rare)", _fnum(safe_attr(pier, "sigma_s_rare"), 1, "MPa"))
    c3.metric("Crack Width", _fnum(safe_attr(pier, "crack_width"), 3, "mm"))

    # Ductile detailing
    st.markdown("#### Ductile Detailing (IRC 112 Cl. 17.2.1)")
    c1, c2, c3 = st.columns(3)
    c1.metric("Spiral Diameter", _fnum(safe_attr(pier, "spiral_dia"), 0, "mm"))
    c2.metric("Spiral Spacing", _fnum(safe_attr(pier, "spiral_spacing"), 0, "mm"))
    c3.metric("\u03c9_wd Provided", _fnum(safe_attr(pier, "omega_wd"), 3))


# ── Piles tab ────────────────────────────────────────────────────────────────

def _render_piles_tab(results: dict) -> None:
    pile_cap = results.get("pile_capacity")
    pile_d = results.get("pile_design")

    st.markdown("### Pile Design")

    # Capacity
    st.markdown("#### Pile Capacity (IRC 78)")
    if pile_cap is not None:
        c1, c2, c3 = st.columns(3)
        c1.metric("Max Compression", _fnum(safe_attr(pile_cap, "max_compression_pile.P",
                   safe_attr(safe_attr(pile_cap, "max_compression_pile"), "P")), 1, "kN"))
        c2.metric("Compression Util", _fnum(safe_attr(pile_cap, "compression_util"), 3),
                  delta=get_status(pile_cap),
                  delta_color="normal" if get_status(pile_cap) == "OK" else "inverse")
        c3.metric("Uplift Util", _fnum(safe_attr(pile_cap, "uplift_util"), 3))

        # Show governing combo
        gov = safe_attr(pile_cap, "governing_compression_combo")
        if gov:
            st.caption(f"Governing: {gov}")
    else:
        st.info("Pile capacity was not computed.")

    # Structural design
    st.markdown("#### Pile Structural Design")
    if pile_d is not None:
        c1, c2, c3 = st.columns(3)
        c1.metric("P-M Utilization", _fnum(safe_attr(pile_d, "util_PM"), 3),
                  delta="OK" if safe_attr(pile_d, "util_PM", 999) <= 1.0 else "FAIL",
                  delta_color="normal" if safe_attr(pile_d, "util_PM", 999) <= 1.0 else "inverse")
        c2.metric("Shear Utilization", _fnum(safe_attr(pile_d, "shear_util"), 3))
        c3.metric("Crack Width", _fnum(safe_attr(pile_d, "crack_width"), 3, "mm"))
    else:
        st.info("Pile structural design was not computed.")


# ── Pile Cap tab ─────────────────────────────────────────────────────────────

def _render_pilecap_tab(results: dict) -> None:
    pc = results.get("pilecap")
    if pc is None:
        st.info("Pile cap design was not computed.")
        return

    st.markdown("### Pile Cap Design")

    # Flexure
    st.markdown("#### Flexure (IRC 112 Cl. 8.2)")
    c1, c2 = st.columns(2)
    c1.metric("Longitudinal Util", _fnum(safe_attr(pc, "flexure_util_long"), 3),
              delta="OK" if safe_attr(pc, "flexure_util_long", 999) <= 1.0 else "FAIL",
              delta_color="normal" if safe_attr(pc, "flexure_util_long", 999) <= 1.0 else "inverse")
    c2.metric("Transverse Util", _fnum(safe_attr(pc, "flexure_util_trans"), 3),
              delta="OK" if safe_attr(pc, "flexure_util_trans", 999) <= 1.0 else "FAIL",
              delta_color="normal" if safe_attr(pc, "flexure_util_trans", 999) <= 1.0 else "inverse")

    # Punching shear
    st.markdown("#### Punching Shear (IRC 112 Cl. 10.4)")
    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Pier Face", _fnum(safe_attr(pc, "punch_util_pier_face"), 3))
    c2.metric("Pier at 2d", _fnum(safe_attr(pc, "punch_util_pier"), 3))
    c3.metric("Pile Face", _fnum(safe_attr(pc, "punch_util_pile_face"), 3))
    c4.metric("Pile at c", _fnum(safe_attr(pc, "punch_util_pile"), 3))

    # Crack width
    st.markdown("#### Crack Width")
    c1, c2 = st.columns(2)
    c1.metric("Crack Width (w_k)", _fnum(safe_attr(pc, "crack_width"), 3, "mm"))
    c2.metric("Limit", _fnum(safe_attr(pc, "crack_width_limit"), 1, "mm"))


# ── Loads tab ────────────────────────────────────────────────────────────────

def _render_loads_tab(results: dict) -> None:
    loads = results.get("loads")
    combos = results.get("combinations")

    st.markdown("### Load Analysis")

    if loads is None:
        st.info("Loads were not computed.")
        return

    # Dead loads
    st.markdown("#### Dead Loads")
    dl_data = {
        "Component": ["Pier Cap", "Pier", "Superstructure (LHS)", "Superstructure (RHS)"],
        "P (kN)": [
            _fnum(safe_attr(loads.dl_piercap, "P"), 1),
            _fnum(safe_attr(loads.dl_pier, "P"), 1),
            _fnum(safe_attr(loads.dl_super_lhs, "P"), 1),
            _fnum(safe_attr(loads.dl_super_rhs, "P"), 1),
        ],
    }
    st.table(dl_data)

    # Live load cases
    st.markdown("#### Live Load Cases")
    if loads.ll_cases:
        ll_data = {"Case": [], "P (kN)": [], "FL (kN)": [], "FT (kN)": [],
                   "ML (kN.m)": [], "MT (kN.m)": []}
        for case in loads.ll_cases:
            fv = case.force_at_pier_base
            ll_data["Case"].append(case.name)
            ll_data["P (kN)"].append(_fnum(fv.P, 1))
            ll_data["FL (kN)"].append(_fnum(fv.FL, 1))
            ll_data["FT (kN)"].append(_fnum(fv.FT, 1))
            ll_data["ML (kN.m)"].append(_fnum(fv.ML, 1))
            ll_data["MT (kN.m)"].append(_fnum(fv.MT, 1))
        st.table(ll_data)

    # Wind
    st.markdown("#### Wind Loads")
    c1, c2 = st.columns(2)
    wind_super = safe_attr(loads, "wind_on_super")
    wind_pier = safe_attr(loads, "wind_on_pier")
    if wind_super:
        c1.metric("Wind on Super (FT)", _fnum(safe_attr(wind_super, "FT"), 1, "kN"))
    if wind_pier:
        c2.metric("Wind on Pier (FT)", _fnum(safe_attr(wind_pier, "FT"), 1, "kN"))

    # Seismic
    st.markdown("#### Seismic")
    c1, c2, c3 = st.columns(3)
    c1.metric("Ah (Long)", _fnum(safe_attr(loads, "Ah_long"), 4))
    c2.metric("Ah (Trans)", _fnum(safe_attr(loads, "Ah_trans"), 4))
    eq_super_long = safe_attr(loads, "eq_super_long")
    if eq_super_long:
        c3.metric("EQ Super (Long, FL)", _fnum(safe_attr(eq_super_long, "FL"), 1, "kN"))

    # Combinations summary
    if combos is not None:
        st.markdown("#### Governing Load Combinations")
        gov_data = {"Envelope": [], "Combination": [], "P (kN)": [],
                    "ML (kN.m)": [], "MT (kN.m)": []}
        for label, combo_attr in [
            ("Max P", "max_P"),
            ("Min P", "min_P"),
            ("Max ML", "max_ML"),
            ("Max MT", "max_MT"),
            ("Gov ULS Basic", "governing_uls_basic"),
            ("Gov ULS Seismic", "governing_uls_seismic"),
        ]:
            combo = safe_attr(combos, combo_attr)
            if combo is not None:
                fv = safe_attr(combo, "forces_pier_base")
                gov_data["Envelope"].append(label)
                gov_data["Combination"].append(safe_attr(combo, "name", "—"))
                gov_data["P (kN)"].append(_fnum(safe_attr(fv, "P"), 1))
                gov_data["ML (kN.m)"].append(_fnum(safe_attr(fv, "ML"), 1))
                gov_data["MT (kN.m)"].append(_fnum(safe_attr(fv, "MT"), 1))
        if gov_data["Envelope"]:
            st.table(gov_data)


# ── BOQ tab ──────────────────────────────────────────────────────────────────

def _render_boq_tab(results: dict) -> None:
    boq = results.get("boq")
    if boq is None:
        st.info("Bill of quantities was not computed.")
        return

    st.markdown("### Bill of Quantities")

    items = safe_attr(boq, "items", [])
    if items:
        data = {
            "Element": [],
            "Concrete (m\u00b3)": [],
            "Steel (kg)": [],
            "Total Cost (INR)": [],
        }
        for item in items:
            data["Element"].append(safe_attr(item, "element", "—"))
            data["Concrete (m\u00b3)"].append(_fnum(safe_attr(item, "concrete_volume"), 2))
            data["Steel (kg)"].append(_fnum(safe_attr(item, "steel_weight"), 1))
            data["Total Cost (INR)"].append(_fnum(safe_attr(item, "total_cost"), 0))
        st.table(data)

    st.markdown("---")
    c1, c2, c3 = st.columns(3)
    c1.metric("Total Concrete", _fnum(safe_attr(boq, "total_concrete"), 2, "m\u00b3"))
    c2.metric("Total Steel", _fnum(safe_attr(boq, "total_steel"), 1, "kg"))
    c3.metric("Total Cost", f"INR {safe_attr(boq, 'total_cost', 0):,.0f}")

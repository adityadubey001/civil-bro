"""IRC 112 Bridge Substructure Design Application

Interactive design tool for bridge piers, pile caps, and pile foundations
per IRC 6:2017 (loads) and IRC 112:2020 (concrete design).
"""

import sys
from pathlib import Path

# Ensure src/ is on the import path
sys.path.insert(0, str(Path(__file__).parent / "src"))

import streamlit as st

from ui.state import init_state, reset_state, get_progress, NUM_STEPS, STEP_NAMES

# ── Page config ──────────────────────────────────────────────────────────────

st.set_page_config(
    page_title="IRC 112 Substructure Designer",
    page_icon="\U0001F3D7\uFE0F",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
    .info-box {
        background-color: #e8f4f8;
        padding: 1rem;
        border-radius: 0.5rem;
        border-left: 4px solid #3498db;
        margin-bottom: 1rem;
    }
    .result-card {
        background-color: #f8f9fa;
        padding: 1rem;
        border-radius: 0.5rem;
        margin-bottom: 1rem;
        border: 1px solid #dee2e6;
    }
    .status-ok { color: #27ae60; font-weight: bold; }
    .status-fail { color: #e74c3c; font-weight: bold; }
</style>
""", unsafe_allow_html=True)


# ── Initialise state ─────────────────────────────────────────────────────────

init_state()


# ── Sidebar ──────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("### IRC 112 Substructure Designer")
    st.markdown("---")

    mode = st.radio(
        "Input Mode",
        ["YAML Upload", "Wizard"],
        index=0 if st.session_state.get("input_mode") == "yaml" else 1,
        help="Choose how to provide design inputs",
    )
    st.session_state["input_mode"] = "yaml" if mode == "YAML Upload" else "wizard"

    # Show wizard progress if in wizard mode
    if st.session_state["input_mode"] == "wizard":
        step = st.session_state.get("current_step", 1)
        st.progress(get_progress(), text=f"Step {step} of {NUM_STEPS}: {STEP_NAMES.get(step, '')}")

    st.markdown("---")

    # Start new design button (always visible)
    if st.session_state.get("results") is not None:
        if st.button("Start New Design", use_container_width=True):
            reset_state()
            st.rerun()

    st.markdown(
        "<small>Design codes: IRC 6:2017, IRC 112:2020<br>"
        "Foundation: IRC 78, IS:2911</small>",
        unsafe_allow_html=True,
    )


# ── Main content ─────────────────────────────────────────────────────────────

def main():
    # If results exist, show them
    if st.session_state.get("results") is not None:
        from ui.results import render_results
        render_results(st.session_state["results"])
        return

    # Otherwise show input mode
    if st.session_state["input_mode"] == "yaml":
        from ui.yaml_mode import render_yaml_mode
        render_yaml_mode()
    else:
        from ui.wizard import render_wizard
        render_wizard()


main()

"""YAML upload / paste input mode for the Substructure Design app."""

from __future__ import annotations

from pathlib import Path

import streamlit as st
import yaml

from ui.runner import run_design


_SAMPLE_PATH = Path(__file__).resolve().parent.parent.parent / "config" / "sample_input.yaml"


def _load_sample() -> str:
    """Read the bundled sample YAML file."""
    try:
        return _SAMPLE_PATH.read_text(encoding="utf-8")
    except FileNotFoundError:
        return "# sample_input.yaml not found"


def render_yaml_mode() -> None:
    """Render the YAML upload / paste interface."""

    st.markdown("## YAML Input")
    st.markdown(
        "Upload a YAML configuration file or paste its contents below. "
        "Click **Load Sample** to start from the bundled example."
    )

    # --- File upload ---
    uploaded = st.file_uploader(
        "Upload YAML file",
        type=["yaml", "yml"],
        help="Upload a .yaml or .yml configuration file",
    )

    # --- Load-sample button ---
    if st.button("Load Sample"):
        st.session_state["yaml_text"] = _load_sample()
        st.rerun()

    # --- Text area ---
    default_text = st.session_state.get("yaml_text", "")
    if uploaded is not None:
        default_text = uploaded.read().decode("utf-8")
        st.session_state["yaml_text"] = default_text

    yaml_text = st.text_area(
        "YAML configuration",
        value=default_text,
        height=400,
        help="Paste or edit your YAML configuration here",
    )
    st.session_state["yaml_text"] = yaml_text

    # --- Action buttons ---
    col1, col2 = st.columns(2)

    with col1:
        validate_clicked = st.button("Validate", use_container_width=True)

    with col2:
        run_clicked = st.button(
            "Run Design", type="primary", use_container_width=True
        )

    # --- Validate ---
    if validate_clicked:
        _validate(yaml_text)

    # --- Run ---
    if run_clicked:
        config = _parse(yaml_text)
        if config is not None:
            with st.spinner("Running substructure design..."):
                results = run_design(config)
            st.session_state["results"] = results
            st.session_state["config"] = config
            st.rerun()


def _parse(text: str) -> dict | None:
    """Parse YAML text and return config dict, or None on error."""
    if not text.strip():
        st.error("YAML input is empty.")
        return None
    try:
        config = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        st.error(f"YAML syntax error: {exc}")
        return None

    if not isinstance(config, dict):
        st.error("YAML content must be a mapping (key-value pairs).")
        return None

    # Basic section check
    required = [
        "project", "superstructure", "pier", "pier_cap",
        "bearings", "levels", "materials", "foundation",
        "seismic", "wind",
    ]
    missing = [s for s in required if s not in config]
    if missing:
        st.error(f"Missing required sections: {', '.join(missing)}")
        return None

    return config


def _validate(text: str) -> None:
    """Parse and validate YAML, showing results inline."""
    config = _parse(text)
    if config is None:
        return  # errors already shown by _parse

    errors = []

    # Levels check
    levels = config.get("levels", {})
    frl = levels.get("frl", 0)
    gl = levels.get("gl", 0)
    if isinstance(frl, (int, float)) and isinstance(gl, (int, float)) and frl <= gl:
        errors.append(f"FRL ({frl}) must be above GL ({gl})")

    # Pier type
    pier_type = config.get("pier", {}).get("type")
    if pier_type not in ("Circular", "Rectangular"):
        errors.append(f"pier.type must be 'Circular' or 'Rectangular', got '{pier_type}'")

    # Concrete grades
    valid_grades = {20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90}
    mat = config.get("materials", {})
    for comp in ("pier", "pier_cap", "pedestal", "pile_pilecap"):
        comp_data = mat.get(comp, {})
        if isinstance(comp_data, dict):
            fck = comp_data.get("fck")
            if fck is not None and fck not in valid_grades:
                errors.append(f"materials.{comp}.fck = {fck} is not a standard IRC 112 grade")

    # Seismic zone
    zone = config.get("seismic", {}).get("zone")
    if zone not in ("II", "III", "IV", "V"):
        errors.append(f"seismic.zone = '{zone}' is invalid (must be II, III, IV, or V)")

    # Exposure
    exposure = mat.get("exposure")
    valid_exp = {"Mild", "Moderate", "Severe", "Very Severe", "Extreme"}
    if exposure not in valid_exp:
        errors.append(f"materials.exposure = '{exposure}' is invalid")

    if errors:
        st.error("Validation errors found:")
        for e in errors:
            st.markdown(f"- {e}")
    else:
        st.success("YAML is valid. Ready to run design.")

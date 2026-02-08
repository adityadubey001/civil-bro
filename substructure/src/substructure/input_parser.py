"""Parse and validate YAML input for bridge substructure design.

Reads a project YAML file, checks that every required section and field is
present, applies defaults for optional fields, and validates value ranges
against Indian code provisions (IRC 112, IS 2911, IRC 6).
"""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any

import yaml


# ---------------------------------------------------------------------------
# Schema definitions
# ---------------------------------------------------------------------------
# Each leaf entry is a tuple:
#   (type, required, default, validator_or_None)
# A validator is a callable (value) -> bool; True means OK.

_VALID_DECK_TYPES = {"U Girder", "Box Girder", "Slab", "I Girder"}
_VALID_CONTINUITY = {"Deck Continuity", "Simply Supported"}
_VALID_PIER_TYPES = {"Circular", "Rectangular"}
_VALID_END_CONDITIONS = {"fixed-free", "fixed-guided", "fixed-fixed"}
_VALID_AGGREGATES = {"Quartzite", "Granite", "Limestone", "Sandstone", "Basalt"}
_VALID_EXPOSURES = {"Mild", "Moderate", "Severe", "Very Severe", "Extreme"}
_VALID_ZONES = {"II", "III", "IV", "V"}
_VALID_SOIL_TYPES = {"hard", "medium", "soft"}

_positive = lambda v: v > 0  # noqa: E731
_non_negative = lambda v: v >= 0  # noqa: E731
_fck_range = lambda v: 20 <= v <= 90  # noqa: E731
_fyk_range = lambda v: 400 <= v <= 600  # noqa: E731
_angle_range = lambda v: 0 <= v <= 90  # noqa: E731
_camber_range = lambda v: -10 <= v <= 10  # noqa: E731
_humidity_range = lambda v: 0 < v <= 100  # noqa: E731
_terrain_range = lambda v: v in (1, 2, 3, 4)  # noqa: E731


def _in_set(valid: set[str]):
    """Return a validator that checks membership in *valid*."""
    return lambda v: v in valid


# Required sections and their fields.  ``None`` as default means the field
# is mandatory.  When the field is itself a dict the entry is nested.
SCHEMA: dict[str, dict[str, tuple]] = {
    "project": {
        "name":     (str,   True,  None, None),
        "pier_id":  (str,   True,  None, None),
        "designer": (str,   False, "",   None),
        "checker":  (str,   False, "",   None),
        "date":     (str,   False, "",   None),
    },
    "superstructure": {
        "deck_width":             (float, True,  None,  _positive),
        "type":                   (str,   True,  None,  _in_set(_VALID_DECK_TYPES)),
        "left_span_cc_bearing":   (float, True,  None,  _positive),
        "right_span_cc_bearing":  (float, True,  None,  _positive),
        "num_girders":            (int,   True,  None,  _positive),
        "girder_spacing":         (float, True,  None,  _positive),
        "cg_below_deck_top":      (float, True,  None,  _positive),
        "depth_incl_slab":        (float, True,  None,  _positive),
        "deck_slab_thickness":    (float, True,  None,  _positive),
        "girder_flange_width":    (float, False, 0.3,   _positive),
        "girder_inertia":         (float, False, 0.0,   _non_negative),
        "continuity":             (str,   True,  None,  _in_set(_VALID_CONTINUITY)),
        "num_continuous_spans":   (int,   False, 1,     _positive),
        "camber_superelevation":  (float, False, 0.0,   _camber_range),
        "bearing_pedestal_depth": (float, True,  None,  _non_negative),
        "skew_angle":             (float, False, 0.0,   _angle_range),
        "radius_of_curvature":   (float, False, 0.0,   _non_negative),
        "wearing_coat_thickness": (float, True,  None,  _non_negative),
        "crash_barrier_width":    (float, False, 0.5,   _non_negative),
        "median_barrier_width":   (float, False, 0.0,   _non_negative),
        "crash_barrier_height":   (float, False, 1.1,   _non_negative),
        "noise_barrier_height":   (float, False, 1.5,   _non_negative),
        "dl_reaction_lhs":        (float, False, 0.0,   None),
        "dl_reaction_rhs":        (float, False, 0.0,   None),
        "sidl_reaction_lhs":      (float, False, 0.0,   None),
        "sidl_reaction_rhs":      (float, False, 0.0,   None),
    },
    "pier": {
        "type":            (str,   True,  None, _in_set(_VALID_PIER_TYPES)),
        "diameter_bottom": (float, True,  None, _positive),
        "diameter_top":    (float, True,  None, _positive),
        "end_condition":   (str,   False, "fixed-free", _in_set(_VALID_END_CONDITIONS)),
        "n_bars":          (int,   False, None, None),
        "bar_dia":         (float, False, None, None),
    },
    "pier_cap": {
        "width_long":          (float, True, None, _positive),
        "length_trans":        (float, True, None, _positive),
        "depth_max":           (float, True, None, _positive),
        "depth_min":           (float, True, None, _positive),
        "bearing_offset_long": (float, True, None, _non_negative),
    },
    "bearings": {
        "type":              (str,   True,  None, None),
        "coordinates_lhs":   (list,  True,  None, None),
        "coordinates_rhs":   (list,  True,  None, None),
        "pedestal_length":   (float, False, 0.8,  _positive),
        "pedestal_width":    (float, False, 0.8,  _positive),
        "pedestal_beyond":   (float, False, 0.2,  _non_negative),
    },
    "levels": {
        "frl": (float, True, None, None),
        "gl":  (float, True, None, None),
        "gwt": (float, True, None, None),
        "hfl": (float, True, None, None),
    },
    "materials": {
        # Nested sub-sections are validated separately.
        "concrete_density": (float, False, 25.0, _positive),
        "exposure":         (str,   True,  None,  _in_set(_VALID_EXPOSURES)),
        "water_density":    (float, False, 10.0,  _positive),
    },
    "creep": {
        "age_of_loading":    (float, False, 7,     _positive),
        "relative_humidity": (float, False, 60.0,  _humidity_range),
        "design_life":       (float, False, 36500, _positive),
    },
    "foundation": {
        "pile_diameter":          (float, True,  None, _positive),
        "piles_long":             (int,   True,  None, _positive),
        "piles_trans":            (int,   True,  None, _positive),
        "pile_spacing_factor":    (float, True,  None, _positive),
        "pilecap_edge_clearance": (float, False, 150,  _non_negative),
        "pilecap_thickness":      (float, True,  None, _positive),
        "pilecap_top_below_gl":   (float, False, 500,  _non_negative),
        "total_pile_depth":       (float, True,  None, _positive),
        "pile_capacity_top":      (float, True,  None, None),
        "pile_uplift_capacity":   (float, False, 0.0,  None),
        "soil_density":           (float, False, 18.0, _positive),
        "friction_angle":         (float, False, 30.0, _angle_range),
        "subgrade_modulus":       (float, True,  None, _positive),
        "fixity_ratio":           (float, False, 2.3,  _positive),
        "cover":                  (float, False, 75.0, _non_negative),
        "pile_n_bars":            (int,   False, None, None),
        "pile_bar_dia":           (float, False, None, None),
    },
    "seismic": {
        "zone":              (str,   True,  None, _in_set(_VALID_ZONES)),
        "importance_factor": (float, False, 1.2,  _positive),
        "response_reduction": (float, False, 4.0, _positive),
        "soil_type":         (str,   False, "medium", _in_set(_VALID_SOIL_TYPES)),
    },
    "wind": {
        "basic_speed":       (float, True,  None, _positive),
        "terrain_category":  (int,   False, 2,    _terrain_range),
    },
    "boq": {
        "concrete_rate": (float, False, 8000,  _positive),
        "steel_rate":    (float, False, 70000, _positive),
    },
}

# Nested material sub-sections (validated independently of the flat schema).
_MATERIAL_SUBS: dict[str, dict[str, tuple]] = {
    "pier": {
        "fck":       (float, True, None, _fck_range),
        "aggregate": (str,   False, "Quartzite", _in_set(_VALID_AGGREGATES)),
    },
    "pier_cap": {
        "fck":       (float, True, None, _fck_range),
        "aggregate": (str,   False, "Quartzite", _in_set(_VALID_AGGREGATES)),
    },
    "pedestal": {
        "fck": (float, True, None, _fck_range),
    },
    "pile_pilecap": {
        "fck": (float, True, None, _fck_range),
    },
    "steel": {
        "fyk": (float, True, None, _fyk_range),
        "Es":  (float, False, 200_000, _positive),
    },
}


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------

class InputError(Exception):
    """Raised when the YAML input is invalid or incomplete."""


def _coerce(value: Any, expected_type: type) -> Any:
    """Attempt to coerce *value* to *expected_type*.

    YAML often reads ``2`` as ``int`` where a ``float`` is expected.  This
    silently promotes ints to floats when the schema says ``float``.
    """
    if expected_type is float and isinstance(value, (int, float)):
        return float(value)
    if expected_type is int and isinstance(value, int):
        return value
    if isinstance(value, expected_type):
        return value
    raise InputError(
        f"Expected type {expected_type.__name__}, got "
        f"{type(value).__name__} for value {value!r}"
    )


def _validate_section(
    data: dict[str, Any],
    schema: dict[str, tuple],
    section_path: str,
    errors: list[str],
) -> dict[str, Any]:
    """Validate *data* against a flat field *schema*.

    Mutates *data* in place: fills defaults and coerces types.  Appends
    human-readable messages to *errors* for every problem found.
    """
    validated: dict[str, Any] = {}
    for field, (ftype, required, default, validator) in schema.items():
        path_str = f"{section_path}.{field}"
        if field not in data:
            if required and default is None:
                errors.append(f"Missing required field: {path_str}")
                continue
            validated[field] = default
            continue

        raw = data[field]
        try:
            coerced = _coerce(raw, ftype)
        except InputError:
            errors.append(
                f"{path_str}: expected {ftype.__name__}, "
                f"got {type(raw).__name__} ({raw!r})"
            )
            continue

        if validator is not None and not validator(coerced):
            errors.append(f"{path_str}: value {coerced!r} is out of range")
            continue

        validated[field] = coerced

    return validated


def _validate_bearing_coords(
    coords: list,
    label: str,
    errors: list[str],
) -> list[list[float]]:
    """Ensure bearing coordinates are a list of ``[x, y]`` pairs."""
    if not isinstance(coords, list) or len(coords) == 0:
        errors.append(f"bearings.{label}: must be a non-empty list of [x, y] pairs")
        return []
    validated: list[list[float]] = []
    for i, pair in enumerate(coords):
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            errors.append(
                f"bearings.{label}[{i}]: each entry must be [x, y], "
                f"got {pair!r}"
            )
            continue
        try:
            validated.append([float(pair[0]), float(pair[1])])
        except (TypeError, ValueError):
            errors.append(
                f"bearings.{label}[{i}]: coordinates must be numeric, "
                f"got {pair!r}"
            )
    return validated


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_input(yaml_path: str) -> dict[str, Any]:
    """Read and validate a project YAML file.

    Parameters
    ----------
    yaml_path:
        Filesystem path to the YAML input file.

    Returns
    -------
    dict
        A fully validated configuration dictionary with defaults applied.

    Raises
    ------
    FileNotFoundError
        If *yaml_path* does not exist.
    InputError
        If validation fails (the message lists every problem found).
    """
    path = Path(yaml_path)
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {yaml_path}")

    with open(path, "r", encoding="utf-8") as fh:
        raw: dict[str, Any] = yaml.safe_load(fh)

    if not isinstance(raw, dict):
        raise InputError("YAML root must be a mapping (dict)")

    errors: list[str] = []
    config: dict[str, Any] = {}

    # ------------------------------------------------------------------
    # 1. Validate top-level sections (flat fields)
    # ------------------------------------------------------------------
    for section_name, field_schema in SCHEMA.items():
        if section_name not in raw:
            # Check whether any field in the section is required.
            has_required = any(
                req and default is None
                for (_, req, default, _) in field_schema.values()
            )
            if has_required:
                errors.append(f"Missing required section: {section_name}")
            # Fill with defaults even when the section is absent.
            config[section_name] = {
                field: default
                for field, (_, _, default, _) in field_schema.items()
            }
            continue

        section_data = raw[section_name]
        if not isinstance(section_data, dict):
            errors.append(f"Section '{section_name}' must be a mapping")
            continue

        config[section_name] = _validate_section(
            section_data, field_schema, section_name, errors
        )

    # ------------------------------------------------------------------
    # 2. Validate nested materials sub-sections
    # ------------------------------------------------------------------
    mat_raw = raw.get("materials", {})
    if not isinstance(mat_raw, dict):
        errors.append("Section 'materials' must be a mapping")
        mat_raw = {}

    for sub_name, sub_schema in _MATERIAL_SUBS.items():
        sub_path = f"materials.{sub_name}"
        if sub_name not in mat_raw:
            has_required = any(
                req and default is None
                for (_, req, default, _) in sub_schema.values()
            )
            if has_required:
                errors.append(f"Missing required sub-section: {sub_path}")
            config["materials"][sub_name] = {
                field: default
                for field, (_, _, default, _) in sub_schema.items()
            }
            continue

        sub_data = mat_raw[sub_name]
        if not isinstance(sub_data, dict):
            errors.append(f"{sub_path}: must be a mapping")
            continue

        config["materials"][sub_name] = _validate_section(
            sub_data, sub_schema, sub_path, errors
        )

    # ------------------------------------------------------------------
    # 3. Validate bearing coordinate lists
    # ------------------------------------------------------------------
    if "bearings" in config:
        brg = config["bearings"]
        if brg.get("coordinates_lhs") is not None:
            brg["coordinates_lhs"] = _validate_bearing_coords(
                brg["coordinates_lhs"], "coordinates_lhs", errors
            )
        if brg.get("coordinates_rhs") is not None:
            brg["coordinates_rhs"] = _validate_bearing_coords(
                brg["coordinates_rhs"], "coordinates_rhs", errors
            )

    # ------------------------------------------------------------------
    # 4. Cross-field sanity checks
    # ------------------------------------------------------------------
    if not errors:
        lvl = config["levels"]
        if lvl["gwt"] > lvl["gl"]:
            errors.append(
                "levels: GWT cannot be above ground level "
                f"(gwt={lvl['gwt']}, gl={lvl['gl']})"
            )
        if lvl["frl"] < lvl["gl"]:
            errors.append(
                "levels: FRL is below ground level "
                f"(frl={lvl['frl']}, gl={lvl['gl']})"
            )

        pc = config["pier_cap"]
        if pc["depth_min"] > pc["depth_max"]:
            errors.append(
                "pier_cap: depth_min cannot exceed depth_max "
                f"(min={pc['depth_min']}, max={pc['depth_max']})"
            )

        fnd = config["foundation"]
        if fnd["pilecap_thickness"] >= fnd["total_pile_depth"]:
            errors.append(
                "foundation: pilecap_thickness must be less than total_pile_depth"
            )

    # ------------------------------------------------------------------
    # Report
    # ------------------------------------------------------------------
    if errors:
        bullet_list = "\n  - ".join(errors)
        raise InputError(
            f"Input validation failed with {len(errors)} error(s):\n"
            f"  - {bullet_list}"
        )

    return config


# ---------------------------------------------------------------------------
# Template generator
# ---------------------------------------------------------------------------

_TEMPLATE_YAML = """\
# Substructure Design Input File
# ================================
# Fill in all values below. Comments show units and valid options.

project:
  name: "PROJECT_NAME"
  pier_id: "P1"
  designer: "Designer Name"
  checker: "Checker Name"
  date: "2024-01-01"

superstructure:
  deck_width: 12.0              # m - Total width of deck
  type: "U Girder"              # Options: U Girder | Box Girder | Slab | I Girder
  left_span_cc_bearing: 40.0    # m - Left span c/c of bearings
  right_span_cc_bearing: 40.0   # m - Right span c/c of bearings
  num_girders: 2                # Number of girders per span
  girder_spacing: 6.0           # m - C/C spacing of girders
  cg_below_deck_top: 1.0        # m - CG of superstructure below deck top
  depth_incl_slab: 3.0          # m - Depth of superstructure including deck slab
  deck_slab_thickness: 0.25     # m
  girder_flange_width: 0.3      # m - Width of top flange
  girder_inertia: 3.0           # m4 - Moment of inertia about centroid
  continuity: "Deck Continuity" # Options: Deck Continuity | Simply Supported
  num_continuous_spans: 2
  camber_superelevation: 2.5    # %
  bearing_pedestal_depth: 0.5   # m
  skew_angle: 0                 # degrees
  radius_of_curvature: 0        # m (use 0 for straight)
  wearing_coat_thickness: 45    # mm
  crash_barrier_width: 0.5      # m
  median_barrier_width: 0       # m (0 if not applicable)
  crash_barrier_height: 1.1     # m
  noise_barrier_height: 1.5     # m
  # Superstructure dead load reactions per bearing (kN)
  # Set to 0 for auto-calculation from geometry
  dl_reaction_lhs: 0
  dl_reaction_rhs: 0
  sidl_reaction_lhs: 0
  sidl_reaction_rhs: 0

pier:
  type: "Circular"              # Options: Circular | Rectangular
  diameter_bottom: 2.0          # m (or width_long for Rectangular)
  diameter_top: 2.0             # m (or length_trans for Rectangular)

pier_cap:
  width_long: 2.85              # m - Width in longitudinal (traffic) direction
  length_trans: 9.2             # m - Length in transverse direction
  depth_max: 1.75               # m - Depth at pier location
  depth_min: 0.75               # m - Depth at far end (cantilever tip)
  bearing_offset_long: 0.825    # m - Distance from pier center to bearing, longitudinal

bearings:
  type: "Elastomeric"
  # Bearing coordinates relative to pier center [along_traffic, across_traffic]
  coordinates_lhs:
    - [-0.825, -3.0]
    - [-0.825,  3.0]
  coordinates_rhs:
    - [ 0.825, -3.0]
    - [ 0.825,  3.0]
  pedestal_length: 0.8          # m
  pedestal_width: 0.8           # m
  pedestal_beyond: 0.2          # m - Portion beyond pedestal

levels:
  frl: 100.000                  # m - Finished Road Level
  gl: 90.000                    # m - Ground Level
  gwt: 90.000                   # m - Ground Water Table Level
  hfl: 90.000                   # m - High Flood Level (= GL for land bridges)

materials:
  pier:
    fck: 50                     # MPa (20-90)
    aggregate: "Quartzite"      # Options: Quartzite | Granite | Limestone | Sandstone | Basalt
  pier_cap:
    fck: 45
    aggregate: "Quartzite"
  pedestal:
    fck: 45
  pile_pilecap:
    fck: 35
  steel:
    fyk: 550                    # MPa
    Es: 200000                  # MPa
  concrete_density: 25          # kN/m3
  exposure: "Severe"            # Options: Mild | Moderate | Severe | Very Severe | Extreme
  water_density: 10             # kN/m3

creep:
  age_of_loading: 7             # days
  relative_humidity: 60         # %
  design_life: 36500            # days (100 years)

foundation:
  pile_diameter: 1.0            # m
  piles_long: 2                 # Number of piles in longitudinal direction
  piles_trans: 2                # Number of piles in transverse direction
  pile_spacing_factor: 3        # C/C spacing = factor * diameter
  pilecap_edge_clearance: 150   # mm - From pile surface to edge of pile cap
  pilecap_thickness: 1.8        # m
  pilecap_top_below_gl: 500     # mm
  total_pile_depth: 15.0        # m - From ground level to founding level
  pile_capacity_top: 500        # tonnes - Safe bearing capacity at pile top
  pile_uplift_capacity: -180    # tonnes - Uplift capacity (negative)
  soil_density: 18              # kN/m3
  friction_angle: 30            # degrees
  subgrade_modulus: 1.4         # MN/m3 (modulus of subgrade reaction)
  fixity_ratio: 2.3             # Lf/T (from IS:2911-2010 Fig 4)
  cover: 75                     # mm - Clear cover for piles

seismic:
  zone: "IV"                    # Options: II | III | IV | V
  importance_factor: 1.2
  response_reduction: 4.0       # R factor
  soil_type: "medium"           # Options: hard | medium | soft

wind:
  basic_speed: 39.0             # m/s - Basic wind speed for the location
  terrain_category: 2           # 1 to 4

boq:
  concrete_rate: 8000           # INR/m3
  steel_rate: 70000             # INR/tonne
"""


def generate_template() -> str:
    """Return a complete sample YAML input template as a string.

    The returned text is ready to be written to a file and edited by the
    user.
    """
    return _TEMPLATE_YAML

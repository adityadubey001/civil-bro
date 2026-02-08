"""Command-line interface for Bridge Substructure Design Tool.

Usage::

    substructure-design run <input_yaml> [-o output_dir]
    substructure-design template
    substructure-design validate <input_yaml>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import click
import yaml


# ---------------------------------------------------------------------------
# Top-level group
# ---------------------------------------------------------------------------

@click.group()
@click.version_option(package_name="substructure-design")
def main():
    """Bridge Substructure Design Tool - IRC 6 & IRC 112."""


# ---------------------------------------------------------------------------
# run
# ---------------------------------------------------------------------------

@main.command()
@click.argument("input_file", type=click.Path(exists=True))
@click.option(
    "-o", "--output",
    default="./output",
    show_default=True,
    help="Output directory for results.",
)
def run(input_file: str, output: str) -> None:
    """Run the full substructure design for INPUT_FILE."""
    input_path = Path(input_file)
    output_dir = Path(output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # 1. Parse input
    # ------------------------------------------------------------------
    click.echo(f"Reading input file: {input_path}")
    try:
        from substructure.input_parser import parse_input
        data = parse_input(input_path)
    except ImportError:
        click.echo("  [input_parser module not found -- loading YAML directly]")
        with open(input_path, encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except Exception as exc:
        click.secho(f"Error parsing input: {exc}", fg="red", err=True)
        raise SystemExit(1) from exc

    project = data.get("project", {})
    pier = data.get("pier", {})
    pier_cap = data.get("pier_cap", {})
    levels = data.get("levels", {})
    foundation = data.get("foundation", {})
    materials_cfg = data.get("materials", {})
    creep_cfg = data.get("creep", {})

    # ------------------------------------------------------------------
    # 2. Geometry
    # ------------------------------------------------------------------
    click.echo("\nCalculating geometry ...")
    try:
        from substructure.geometry import calculate_geometry
        geometry = calculate_geometry(data)
        click.echo("  Geometry calculated.")
    except ImportError:
        click.echo("  [geometry module not yet implemented]")
        geometry = None
    except Exception as exc:
        click.secho(f"  Geometry error: {exc}", fg="red", err=True)
        geometry = None

    # ------------------------------------------------------------------
    # 3. Materials
    # ------------------------------------------------------------------
    click.echo("Calculating material properties ...")
    material_props = {}
    try:
        from substructure.materials import get_concrete_properties, get_steel_properties

        for comp in ("pier", "pier_cap", "pile_pilecap"):
            comp_cfg = materials_cfg.get(comp, {})
            if isinstance(comp_cfg, dict) and "fck" in comp_cfg:
                cp = get_concrete_properties(
                    comp_cfg["fck"],
                    aggregate=comp_cfg.get("aggregate", "Quartzite"),
                    exposure=materials_cfg.get("exposure"),
                )
                material_props[comp] = cp
                click.echo(f"  {comp}: fck={cp.fck}, fcd={cp.fcd:.2f}, Ecm={cp.Ecm:.2f} GPa")

        steel_cfg = materials_cfg.get("steel", {})
        fyk = steel_cfg.get("fyk", 550)
        sp = get_steel_properties(float(fyk))
        material_props["steel"] = sp
        click.echo(f"  Steel: fyk={sp.fyk}, fyd={sp.fyd:.2f} MPa")
    except Exception as exc:
        click.secho(f"  Materials error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 4. Creep
    # ------------------------------------------------------------------
    click.echo("Calculating creep coefficients ...")
    creep_results = None
    if geometry is not None:
        try:
            from substructure.creep import calculate_creep, calculate_notional_size
            from substructure.materials import get_concrete_properties as _gcp

            pier_cp = _gcp(
                materials_cfg.get("pier", {}).get("fck", 50),
                aggregate=materials_cfg.get("pier", {}).get("aggregate", "Quartzite"),
            )
            ps = geometry.pier_section
            area_mm2 = ps.area * 1e6
            peri_mm = ps.perimeter * 1e3
            h0 = calculate_notional_size(area_mm2, peri_mm)

            creep_results = calculate_creep(
                fcm=pier_cp.fcm,
                Ecm=pier_cp.Ecm,
                Es=200.0,
                h0=h0,
                RH=float(creep_cfg.get("relative_humidity", 70)),
                t0=float(creep_cfg.get("age_of_loading", 7)),
                t=float(creep_cfg.get("design_life", 36500)),
            )
            click.echo(f"  Creep coefficient: phi = {creep_results.phi:.3f}")
            click.echo(f"  Ec_eff = {creep_results.Ec_eff:.2f} GPa")
        except Exception as exc:
            click.secho(f"  Creep error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 5. Loads
    # ------------------------------------------------------------------
    click.echo("Calculating loads ...")
    load_results = None
    if geometry is not None:
        try:
            from substructure.loads import calculate_loads
            load_results = calculate_loads(data, geometry)
            click.echo(f"  DL pier: {load_results.dl_pier.P:.2f} kN")
            click.echo(f"  DL piercap: {load_results.dl_piercap.P:.2f} kN")
            click.echo(f"  Total DL at pier base: {load_results.total_dl_at_pier_base.P:.2f} kN")
            click.echo(f"  LL cases: {len(load_results.ll_cases)}")
        except Exception as exc:
            click.secho(f"  Loads error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 6. Load Combinations
    # ------------------------------------------------------------------
    click.echo("Generating load combinations ...")
    combo_results = None
    if load_results is not None and geometry is not None:
        try:
            from substructure.load_combinations import calculate_combinations, summarise_governing
            combo_results = calculate_combinations(data, geometry, load_results)
            n = len(combo_results.all_combinations)
            click.echo(f"  {n} combinations generated.")
            click.echo(f"  Max P: {combo_results.max_P.forces_pier_base.P:.2f} kN")
            click.echo(f"  Min P: {combo_results.min_P.forces_pier_base.P:.2f} kN")
        except Exception as exc:
            click.secho(f"  Load combinations error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 7. Pier Cap Design
    # ------------------------------------------------------------------
    click.echo("Designing pier cap ...")
    piercap_result = None
    if combo_results is not None:
        try:
            from substructure.pier_cap_design import design_pier_cap
            piercap_result = design_pier_cap(data, geometry, combo_results, creep_results)
            click.echo(f"  Flexure util (pier face): {piercap_result.flexure_util_pier:.3f}")
            click.echo(f"  Shear util: {piercap_result.shear_util:.3f} "
                       f"(VEd={piercap_result.VEd_pier:.1f} kN, "
                       f"VRds={piercap_result.VRds:.1f} kN)")
            click.echo(f"  Crack width: {piercap_result.crack_width:.3f} mm "
                       f"(limit {piercap_result.crack_width_limit:.1f} mm)")
            click.echo(f"  Status: {piercap_result.status}")
        except ImportError:
            click.echo("  [pier_cap_design module not yet implemented]")
        except Exception as exc:
            click.secho(f"  Pier cap design error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 8. Pier Design
    # ------------------------------------------------------------------
    click.echo("Designing pier ...")
    pier_result = None
    if combo_results is not None:
        try:
            from substructure.pier_design import design_pier
            pier_result = design_pier(data, geometry, combo_results)
            click.echo(f"  Biaxial util: {pier_result.util_biaxial:.3f}")
            click.echo(f"  Governing: {pier_result.governing_combo_name}")
            click.echo(f"  Status: {pier_result.status}")
        except ImportError:
            click.echo("  [pier_design module not yet implemented]")
        except Exception as exc:
            click.secho(f"  Pier design error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 9. Pile Capacity
    # ------------------------------------------------------------------
    click.echo("Calculating pile capacity ...")
    pile_cap_result = None
    if combo_results is not None:
        try:
            from substructure.pile_capacity import calculate_pile_capacity
            pile_cap_result = calculate_pile_capacity(data, geometry, combo_results)
            click.echo(f"  Max pile compression: {pile_cap_result.max_compression_pile.P:.1f} kN "
                       f"(util {pile_cap_result.compression_util:.3f})")
            click.echo(f"  Max pile uplift:      {pile_cap_result.max_uplift_pile.P:.1f} kN "
                       f"(util {pile_cap_result.uplift_util:.3f})")
            click.echo(f"  Status: {pile_cap_result.status}")
        except ImportError:
            click.echo("  [pile_capacity module not yet implemented]")
        except Exception as exc:
            click.secho(f"  Pile capacity error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 10. Pile Design
    # ------------------------------------------------------------------
    click.echo("Designing piles ...")
    pile_design_result = None
    if pile_cap_result is not None:
        try:
            from substructure.pile_design import design_pile
            pile_design_result = design_pile(data, geometry, combo_results, pile_cap_result)
            click.echo(f"  P-M util: {pile_design_result.util_PM:.3f}")
            click.echo(f"  Shear util: {pile_design_result.shear_util:.3f}")
            click.echo(f"  Crack width: {pile_design_result.crack_width:.4f} mm")
            click.echo(f"  Status: {pile_design_result.status}")
        except ImportError:
            click.echo("  [pile_design module not yet implemented]")
        except Exception as exc:
            click.secho(f"  Pile design error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 11. Pile Cap Design
    # ------------------------------------------------------------------
    click.echo("Designing pile cap ...")
    pilecap_design_result = None
    if pile_cap_result is not None:
        try:
            from substructure.pilecap_design import design_pilecap
            pilecap_design_result = design_pilecap(data, geometry, combo_results, pile_cap_result, creep_results)
            click.echo(f"  Flexure util (long): {pilecap_design_result.flexure_util_long:.3f}")
            click.echo(f"  Flexure util (trans): {pilecap_design_result.flexure_util_trans:.3f}")
            click.echo(f"  Punching (pier): {pilecap_design_result.punch_util_pier:.3f}")
            click.echo(f"  Punching (pile): {pilecap_design_result.punch_util_pile:.3f}")
            click.echo(f"  Status: {pilecap_design_result.status}")
        except ImportError:
            click.echo("  [pilecap_design module not yet implemented]")
        except Exception as exc:
            click.secho(f"  Pile cap design error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # 12. Bill of Quantities
    # ------------------------------------------------------------------
    click.echo("Calculating BOQ ...")
    boq_result = None
    if piercap_result and pier_result and pile_design_result and pilecap_design_result:
        try:
            from substructure.boq import calculate_boq
            boq_result = calculate_boq(
                data, geometry, piercap_result, pier_result,
                pile_design_result, pilecap_design_result,
            )
            click.echo(f"  Total concrete: {boq_result.total_concrete:.2f} m3")
            click.echo(f"  Total steel: {boq_result.total_steel:.1f} kg")
            click.echo(f"  Total cost: INR {boq_result.total_cost:,.0f}")
        except Exception as exc:
            click.secho(f"  BOQ error: {exc}", fg="red", err=True)

    # ------------------------------------------------------------------
    # Console summary
    # ------------------------------------------------------------------
    click.echo("")
    click.secho("=" * 60, bold=True)
    click.secho("  SUBSTRUCTURE DESIGN SUMMARY", bold=True)
    click.secho("=" * 60, bold=True)

    click.echo(f"\n  Project       : {project.get('name', 'N/A')}")
    click.echo(f"  Pier ID       : {project.get('pier_id', 'N/A')}")

    def _status_line(label, result_obj, status_attr="status"):
        val = getattr(result_obj, status_attr, None) if result_obj else None
        colour = "green" if val == "OK" else "red" if val else "yellow"
        click.echo(f"  {label:<26}: ", nl=False)
        click.secho(val or "not computed", fg=colour)

    click.echo(f"\n  --- Design Status ---")
    _status_line("Pier Cap", piercap_result)
    _status_line("Pier", pier_result)
    _status_line("Pile Capacity", pile_cap_result)
    _status_line("Pile Design", pile_design_result)
    _status_line("Pile Cap Design", pilecap_design_result)

    if boq_result:
        click.echo(f"\n  --- Cost ---")
        click.echo(f"  Total cost              : INR {boq_result.total_cost:,.0f}")

    click.secho("\n" + "=" * 60, bold=True)

    # ------------------------------------------------------------------
    # Save results
    # ------------------------------------------------------------------
    results = {
        "project": project,
        "geometry": geometry,
        "materials": material_props,
        "creep": creep_results,
        "loads": load_results,
        "combinations": combo_results,
        "pier_cap_design": piercap_result,
        "pier_design": pier_result,
        "pile_capacity": pile_cap_result,
        "pile_design": pile_design_result,
        "pilecap_design": pilecap_design_result,
        "boq": boq_result,
    }

    results_file = output_dir / "results.json"
    serialisable = {}
    for key, val in results.items():
        if val is None:
            serialisable[key] = "not yet computed"
        else:
            serialisable[key] = val

    with open(results_file, "w", encoding="utf-8") as fh:
        json.dump(serialisable, fh, indent=2, default=str)

    click.echo(f"\nResults saved to {results_file.resolve()}")

    # ------------------------------------------------------------------
    # 13. PDF Report
    # ------------------------------------------------------------------
    click.echo("Generating PDF report ...")
    try:
        from substructure.report import generate_report
        pdf_path = output_dir / f"{project.get('name', 'report')}.pdf"
        generate_report(str(pdf_path), data, results)
        click.secho(f"Report saved to {pdf_path.resolve()}", fg="green")
    except Exception as exc:
        click.secho(f"  Report error: {exc}", fg="red", err=True)


# ---------------------------------------------------------------------------
# template
# ---------------------------------------------------------------------------

_SAMPLE_YAML = """\
# Substructure Design Input File
# ================================
# Fill in all values below. Comments show units and valid options.

project:
  name: "MY-PROJECT-001"
  pier_id: "P1"
  designer: "Design Engineer"
  checker: "Checker"
  date: "2025-01-01"

superstructure:
  deck_width: 12.0              # m - Total width of deck
  type: "U Girder"              # Options: U Girder | Box Girder | Slab | I Girder
  left_span_cc_bearing: 58.35   # m - Left span c/c of bearings
  right_span_cc_bearing: 58.35  # m - Right span c/c of bearings
  num_girders: 2                # Number of girders per span
  girder_spacing: 6.0           # m - C/C spacing of girders
  cg_below_deck_top: 1.023      # m - CG of superstructure below deck top
  depth_incl_slab: 3.23         # m - Depth of superstructure including deck slab
  deck_slab_thickness: 0.25     # m
  girder_flange_width: 0.3      # m - Width of top flange
  girder_inertia: 3.245         # m4 - Moment of inertia about centroid
  continuity: "Deck Continuity" # Options: Deck Continuity | Simply Supported
  num_continuous_spans: 2
  camber_superelevation: 2.5    # %
  bearing_pedestal_depth: 0.5   # m
  skew_angle: 0                 # degrees
  radius_of_curvature: 2000     # m (use 0 for straight)
  wearing_coat_thickness: 45    # mm
  crash_barrier_width: 0.5      # m
  median_barrier_width: 0       # m (0 if not applicable)
  crash_barrier_height: 1.1     # m
  noise_barrier_height: 1.5     # m
  dl_reaction_lhs: 0            # kN (0 = auto-calculate)
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
  coordinates_lhs:
    - [-0.825, -3.7]
    - [-0.825, -2.3]
    - [-0.825,  2.3]
    - [-0.825,  3.7]
  coordinates_rhs:
    - [ 0.825, -3.7]
    - [ 0.825, -2.3]
    - [ 0.825,  2.3]
    - [ 0.825,  3.7]
  pedestal_length: 0.8          # m
  pedestal_width: 0.8           # m
  pedestal_beyond: 0.2          # m

levels:
  frl: 315.685                  # m - Finished Road Level
  gl: 303.236                   # m - Ground Level
  gwt: 303.236                  # m - Ground Water Table Level
  hfl: 303.236                  # m - High Flood Level (= GL for land bridges)

materials:
  pier:
    fck: 50                     # MPa
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
  relative_humidity: 54.5       # %
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
  pile_capacity_top: 550        # tonnes - Safe bearing capacity at pile top
  pile_uplift_capacity: -187    # tonnes - Uplift capacity (negative)
  soil_density: 20              # kN/m3
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


@main.command()
def template() -> None:
    """Print a sample input YAML to stdout."""
    click.echo(_SAMPLE_YAML)


# ---------------------------------------------------------------------------
# validate
# ---------------------------------------------------------------------------

_REQUIRED_SECTIONS = [
    "project",
    "superstructure",
    "pier",
    "pier_cap",
    "bearings",
    "levels",
    "materials",
    "creep",
    "foundation",
    "seismic",
    "wind",
]

_REQUIRED_PROJECT_KEYS = ["name", "pier_id"]
_REQUIRED_LEVEL_KEYS = ["frl", "gl"]


def _validate_data(data: dict) -> list[str]:
    """Return a list of validation error strings (empty means valid)."""
    errors: list[str] = []

    if not isinstance(data, dict):
        return ["Input file does not contain a valid YAML mapping."]

    # Check top-level sections
    for section in _REQUIRED_SECTIONS:
        if section not in data:
            errors.append(f"Missing required section: '{section}'")

    # Project
    project = data.get("project", {})
    if isinstance(project, dict):
        for key in _REQUIRED_PROJECT_KEYS:
            if key not in project:
                errors.append(f"Missing 'project.{key}'")

    # Levels
    levels = data.get("levels", {})
    if isinstance(levels, dict):
        for key in _REQUIRED_LEVEL_KEYS:
            if key not in levels:
                errors.append(f"Missing 'levels.{key}'")
        frl = levels.get("frl", 0)
        gl = levels.get("gl", 0)
        if isinstance(frl, (int, float)) and isinstance(gl, (int, float)):
            if frl <= gl:
                errors.append(
                    f"FRL ({frl}) must be above GL ({gl})"
                )

    # Pier
    pier = data.get("pier", {})
    if isinstance(pier, dict):
        ptype = pier.get("type")
        if ptype not in ("Circular", "Rectangular"):
            errors.append(
                f"pier.type must be 'Circular' or 'Rectangular', got '{ptype}'"
            )

    # Materials - concrete grades
    mat = data.get("materials", {})
    if isinstance(mat, dict):
        valid_grades = {20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 90}
        for component in ("pier", "pier_cap", "pedestal", "pile_pilecap"):
            comp_data = mat.get(component, {})
            if isinstance(comp_data, dict):
                fck = comp_data.get("fck")
                if fck is not None and fck not in valid_grades:
                    errors.append(
                        f"materials.{component}.fck = {fck} is not a standard "
                        f"IRC 112 grade. Valid: {sorted(valid_grades)}"
                    )

        exposure = mat.get("exposure")
        valid_exposures = {"Mild", "Moderate", "Severe", "Very Severe", "Extreme"}
        if exposure not in valid_exposures:
            errors.append(
                f"materials.exposure = '{exposure}' is invalid. "
                f"Valid: {sorted(valid_exposures)}"
            )

    # Seismic
    seismic = data.get("seismic", {})
    if isinstance(seismic, dict):
        zone = seismic.get("zone")
        if zone not in ("II", "III", "IV", "V"):
            errors.append(
                f"seismic.zone = '{zone}' is invalid. Valid: II, III, IV, V"
            )

    # Foundation
    fnd = data.get("foundation", {})
    if isinstance(fnd, dict):
        pd = fnd.get("pile_diameter", 0)
        if isinstance(pd, (int, float)) and pd <= 0:
            errors.append("foundation.pile_diameter must be > 0")

    return errors


@main.command()
@click.argument("input_file", type=click.Path(exists=True))
def validate(input_file: str) -> None:
    """Validate an input YAML file without running the design."""
    input_path = Path(input_file)
    click.echo(f"Validating: {input_path}")

    # --- Basic YAML parse ---
    try:
        with open(input_path, encoding="utf-8") as fh:
            data = yaml.safe_load(fh)
    except yaml.YAMLError as exc:
        click.secho(f"YAML syntax error:\n  {exc}", fg="red", err=True)
        raise SystemExit(1) from exc

    # --- Try module-level parser first (may add its own checks) ---
    try:
        from substructure.input_parser import parse_input
        data = parse_input(input_path)
        click.echo("  input_parser.parse_input() succeeded.")
    except ImportError:
        click.echo("  [input_parser module not found -- using built-in checks]")
    except Exception as exc:
        click.secho(f"  input_parser error: {exc}", fg="red", err=True)
        raise SystemExit(1) from exc

    # --- Built-in validation ---
    errors = _validate_data(data)
    if errors:
        click.secho(f"\nFound {len(errors)} issue(s):\n", fg="yellow")
        for i, err in enumerate(errors, 1):
            click.echo(f"  {i}. {err}")
        raise SystemExit(1)
    else:
        click.secho("\nInput file is valid.", fg="green")


# ---------------------------------------------------------------------------
# Allow ``python -m substructure.cli``
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    main()

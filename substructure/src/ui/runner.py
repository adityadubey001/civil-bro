"""Design engine wrapper for the Streamlit app.

Replicates the cli.py pipeline (steps 2-12) but returns structured
results in memory instead of printing to console or writing files.
"""

from __future__ import annotations

from dataclasses import asdict
from typing import Any


def run_design(config: dict) -> dict:
    """Run the full substructure design pipeline.

    Parameters
    ----------
    config : dict
        YAML-compatible configuration dict (same structure as sample_input.yaml).

    Returns
    -------
    dict with keys:
        geometry, materials, creep, loads, combinations,
        pier_cap, pier, pile_capacity, pile_design, pilecap, boq,
        errors (list[str]), warnings (list[str]),
        steps_completed (list[str])
    """
    results: dict[str, Any] = {
        "geometry": None,
        "materials": {},
        "creep": None,
        "loads": None,
        "combinations": None,
        "pier_cap": None,
        "pier": None,
        "pile_capacity": None,
        "pile_design": None,
        "pilecap": None,
        "boq": None,
        "errors": [],
        "warnings": [],
        "steps_completed": [],
    }

    materials_cfg = config.get("materials", {})
    creep_cfg = config.get("creep", {})

    # ------------------------------------------------------------------
    # 1. Geometry
    # ------------------------------------------------------------------
    try:
        from substructure.geometry import calculate_geometry
        geometry = calculate_geometry(config)
        results["geometry"] = geometry
        results["steps_completed"].append("geometry")
    except Exception as exc:
        results["errors"].append(f"Geometry: {exc}")
        return results  # Can't continue without geometry

    # ------------------------------------------------------------------
    # 2. Materials
    # ------------------------------------------------------------------
    try:
        from substructure.materials import get_concrete_properties, get_steel_properties

        material_props = {}
        for comp in ("pier", "pier_cap", "pile_pilecap"):
            comp_cfg = materials_cfg.get(comp, {})
            if isinstance(comp_cfg, dict) and "fck" in comp_cfg:
                cp = get_concrete_properties(
                    comp_cfg["fck"],
                    aggregate=comp_cfg.get("aggregate", "Quartzite"),
                    exposure=materials_cfg.get("exposure"),
                )
                material_props[comp] = cp

        steel_cfg = materials_cfg.get("steel", {})
        fyk = steel_cfg.get("fyk", 550)
        sp = get_steel_properties(float(fyk))
        material_props["steel"] = sp

        results["materials"] = material_props
        results["steps_completed"].append("materials")
    except Exception as exc:
        results["errors"].append(f"Materials: {exc}")

    # ------------------------------------------------------------------
    # 3. Creep
    # ------------------------------------------------------------------
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

        creep_result = calculate_creep(
            fcm=pier_cp.fcm,
            Ecm=pier_cp.Ecm,
            Es=200.0,
            h0=h0,
            RH=float(creep_cfg.get("relative_humidity", 70)),
            t0=float(creep_cfg.get("age_of_loading", 7)),
            t=float(creep_cfg.get("design_life", 36500)),
        )
        results["creep"] = creep_result
        results["steps_completed"].append("creep")
    except Exception as exc:
        results["warnings"].append(f"Creep: {exc}")
        creep_result = None

    # ------------------------------------------------------------------
    # 4. Loads
    # ------------------------------------------------------------------
    try:
        from substructure.loads import calculate_loads
        load_results = calculate_loads(config, geometry)
        results["loads"] = load_results
        results["steps_completed"].append("loads")
    except Exception as exc:
        results["errors"].append(f"Loads: {exc}")
        return results  # Can't continue without loads

    # ------------------------------------------------------------------
    # 5. Load Combinations
    # ------------------------------------------------------------------
    try:
        from substructure.load_combinations import calculate_combinations
        combo_results = calculate_combinations(config, geometry, load_results)
        results["combinations"] = combo_results
        results["steps_completed"].append("combinations")
    except Exception as exc:
        results["errors"].append(f"Load combinations: {exc}")
        return results  # Can't continue without combinations

    # ------------------------------------------------------------------
    # 6. Pier Cap Design
    # ------------------------------------------------------------------
    try:
        from substructure.pier_cap_design import design_pier_cap
        piercap_result = design_pier_cap(config, geometry, combo_results, creep_result)
        results["pier_cap"] = piercap_result
        results["steps_completed"].append("pier_cap")
    except Exception as exc:
        results["errors"].append(f"Pier cap design: {exc}")
        piercap_result = None

    # ------------------------------------------------------------------
    # 7. Pier Design
    # ------------------------------------------------------------------
    try:
        from substructure.pier_design import design_pier
        pier_result = design_pier(config, geometry, combo_results)
        results["pier"] = pier_result
        results["steps_completed"].append("pier")
    except Exception as exc:
        results["errors"].append(f"Pier design: {exc}")
        pier_result = None

    # ------------------------------------------------------------------
    # 8. Pile Capacity
    # ------------------------------------------------------------------
    try:
        from substructure.pile_capacity import calculate_pile_capacity
        pile_cap_result = calculate_pile_capacity(config, geometry, combo_results)
        results["pile_capacity"] = pile_cap_result
        results["steps_completed"].append("pile_capacity")
    except Exception as exc:
        results["errors"].append(f"Pile capacity: {exc}")
        pile_cap_result = None

    # ------------------------------------------------------------------
    # 9. Pile Design
    # ------------------------------------------------------------------
    if pile_cap_result is not None:
        try:
            from substructure.pile_design import design_pile
            pile_design_result = design_pile(
                config, geometry, combo_results, pile_cap_result
            )
            results["pile_design"] = pile_design_result
            results["steps_completed"].append("pile_design")
        except Exception as exc:
            results["errors"].append(f"Pile design: {exc}")
            pile_design_result = None
    else:
        pile_design_result = None

    # ------------------------------------------------------------------
    # 10. Pile Cap Design
    # ------------------------------------------------------------------
    if pile_cap_result is not None:
        try:
            from substructure.pilecap_design import design_pilecap
            pilecap_design_result = design_pilecap(
                config, geometry, combo_results, pile_cap_result, creep_result
            )
            results["pilecap"] = pilecap_design_result
            results["steps_completed"].append("pilecap")
        except Exception as exc:
            results["errors"].append(f"Pile cap design: {exc}")
            pilecap_design_result = None
    else:
        pilecap_design_result = None

    # ------------------------------------------------------------------
    # 11. Bill of Quantities
    # ------------------------------------------------------------------
    if all(r is not None for r in [piercap_result, pier_result,
                                    pile_design_result, pilecap_design_result]):
        try:
            from substructure.boq import calculate_boq
            boq_result = calculate_boq(
                config, geometry, piercap_result, pier_result,
                pile_design_result, pilecap_design_result,
            )
            results["boq"] = boq_result
            results["steps_completed"].append("boq")
        except Exception as exc:
            results["warnings"].append(f"BOQ: {exc}")

    return results


def safe_attr(obj: Any, attr: str, default: Any = None) -> Any:
    """Safely get an attribute from a result object."""
    if obj is None:
        return default
    return getattr(obj, attr, default)


def get_status(obj: Any) -> str | None:
    """Get status string from a result object."""
    if obj is None:
        return None
    return getattr(obj, "status", None)


def all_checks_ok(results: dict) -> bool:
    """Return True if all design checks passed."""
    checks = ["pier_cap", "pier", "pile_capacity", "pile_design", "pilecap"]
    for key in checks:
        status = get_status(results.get(key))
        if status is None or status != "OK":
            return False
    return True

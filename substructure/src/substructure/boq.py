"""
Bill of Quantities (BOQ) calculation module.

Calculates concrete volumes and steel weights for all structural elements,
then computes cost estimates.
"""

from dataclasses import dataclass
from typing import Any
import math

from .geometry import GeometryResults


@dataclass
class BOQItem:
    """Represents a single BOQ line item for a structural element."""
    element: str            # e.g. "Pier", "Pier Cap", "Pile Cap", "Piles"
    concrete_volume: float  # m3
    steel_weight: float     # kg
    concrete_cost: float    # INR
    steel_cost: float       # INR
    total_cost: float       # INR


@dataclass
class BOQResult:
    """Complete BOQ with all items and totals."""
    items: list[BOQItem]
    total_concrete: float   # m3
    total_steel: float      # kg
    total_cost: float       # INR


def calculate_boq(
    config: dict,
    geometry: GeometryResults,
    pier_cap_result: Any,
    pier_result: Any,
    pile_design_result: Any,
    pilecap_result: Any,
) -> BOQResult:
    """
    Calculate Bill of Quantities for all structural elements.

    Args:
        config: Configuration dictionary
        geometry: Geometry calculation results
        pier_cap_result: PierCapDesignResult
        pier_result: PierDesignResult
        pile_design_result: PileDesignResult
        pilecap_result: PilecapDesignResult

    Returns:
        BOQResult with itemized and total quantities
    """
    boq_cfg = config.get("boq", {})
    concrete_rate = boq_cfg.get("concrete_rate", 8000)    # INR/m3
    steel_rate = boq_cfg.get("steel_rate", 70000)         # INR/tonne

    items: list[BOQItem] = []

    # --- Pier ---
    pier_concrete = _pier_concrete(config, geometry)
    pier_steel = _pier_steel(geometry, pier_result)
    items.append(_make_item("Pier", pier_concrete, pier_steel, concrete_rate, steel_rate))

    # --- Pier Cap ---
    pc_concrete = _pier_cap_concrete(config)
    pc_steel = _pier_cap_steel(config, pier_cap_result)
    items.append(_make_item("Pier Cap", pc_concrete, pc_steel, concrete_rate, steel_rate))

    # --- Pile Cap ---
    pcap_concrete = _pilecap_concrete(geometry)
    pcap_steel = _pilecap_steel(geometry, pilecap_result)
    items.append(_make_item("Pile Cap", pcap_concrete, pcap_steel, concrete_rate, steel_rate))

    # --- Piles ---
    piles_concrete = _piles_concrete(geometry)
    piles_steel = _piles_steel(geometry, pile_design_result)
    items.append(_make_item("Piles", piles_concrete, piles_steel, concrete_rate, steel_rate))

    total_concrete = sum(i.concrete_volume for i in items)
    total_steel = sum(i.steel_weight for i in items)
    total_cost = sum(i.total_cost for i in items)

    return BOQResult(
        items=items,
        total_concrete=total_concrete,
        total_steel=total_steel,
        total_cost=total_cost,
    )


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

_STEEL_DENSITY = 7850.0  # kg/m3


def _make_item(name: str, concrete_m3: float, steel_kg: float,
               concrete_rate: float, steel_rate: float) -> BOQItem:
    cc = concrete_m3 * concrete_rate
    sc = (steel_kg / 1000.0) * steel_rate
    return BOQItem(
        element=name,
        concrete_volume=concrete_m3,
        steel_weight=steel_kg,
        concrete_cost=cc,
        steel_cost=sc,
        total_cost=cc + sc,
    )


# --- Concrete volumes ---

def _pier_concrete(config: dict, geometry: GeometryResults) -> float:
    """Pier concrete volume in m3.  Geometry values are in SI (metres)."""
    ps = geometry.pier_section
    h = geometry.pier_height  # m
    if ps.pier_type == "Circular":
        d = ps.width_long  # m (diameter)
        return (math.pi / 4) * d ** 2 * h
    else:
        return ps.width_long * ps.width_trans * h


def _pier_cap_concrete(config: dict) -> float:
    """Pier cap concrete volume in m3 (trapezoidal / haunched)."""
    pc = config.get("pier_cap", {})
    L = pc.get("length_trans", 0)    # m
    W = pc.get("width_long", 0)      # m
    d_max = pc.get("depth_max", 0)   # m
    d_min = pc.get("depth_min", 0)   # m
    return L * W * (d_max + d_min) / 2.0


def _pilecap_concrete(geometry: GeometryResults) -> float:
    """Pile cap concrete volume in m3.  Geometry stores metres."""
    return (geometry.pilecap_length_long
            * geometry.pilecap_width_trans
            * geometry.pilecap_thickness)


def _piles_concrete(geometry: GeometryResults) -> float:
    """Total concrete volume for all piles in m3."""
    d = geometry.pile_props.diameter        # m
    L = geometry.pile_props.pile_length     # m
    n = geometry.num_piles
    return (math.pi / 4) * d ** 2 * L * n


# --- Steel weights ---

def _pier_steel(geometry: GeometryResults, pier_result: Any) -> float:
    """Pier steel weight in kg.  pier_result fields: Ast_provided (mm2), height (m)."""
    Ast_mm2 = pier_result.Ast_provided          # mm2
    h_mm = pier_result.height * 1000.0           # m -> mm
    main_wt = Ast_mm2 * h_mm / 1e9 * _STEEL_DENSITY   # kg
    return main_wt * 1.15   # +15 % for spirals / ties / laps


def _pier_cap_steel(config: dict, pier_cap_result: Any) -> float:
    """Pier cap steel weight in kg."""
    # Main flexural steel at pier face
    Ast = pier_cap_result.Ast_provided_pier   # mm2
    W = pier_cap_result.width                 # mm (longitudinal width)
    main_wt = Ast * W / 1e9 * _STEEL_DENSITY

    # Shear reinforcement
    Asw_s = pier_cap_result.Asw_s_provided    # mm2/mm
    L_trans = config.get("pier_cap", {}).get("length_trans", 0) * 1000  # m -> mm
    shear_wt = Asw_s * L_trans / 1e9 * _STEEL_DENSITY

    return (main_wt + shear_wt) * 1.20  # +20 % secondary / laps


def _pilecap_steel(geometry: GeometryResults, pilecap_result: Any) -> float:
    """Pile cap steel weight in kg."""
    # Bottom mat — longitudinal bars
    Ast_L = pilecap_result.Ast_prov_long   # mm2
    W_t = pilecap_result.width_trans       # mm
    vol_L = Ast_L * W_t / 1e9

    # Bottom mat — transverse bars
    Ast_T = pilecap_result.Ast_prov_trans  # mm2
    L_l = pilecap_result.length_long       # mm
    vol_T = Ast_T * L_l / 1e9

    bottom_wt = (vol_L + vol_T) * _STEEL_DENSITY
    return bottom_wt * 1.15  # +15 % top mat / side faces / laps


def _piles_steel(geometry: GeometryResults, pile_result: Any) -> float:
    """Steel weight for all piles in kg."""
    Ast = pile_result.Ast_provided              # mm2
    L_mm = geometry.pile_props.pile_length * 1000  # m -> mm
    n = geometry.num_piles
    per_pile = Ast * L_mm / 1e9 * _STEEL_DENSITY
    return per_pile * n * 1.15  # +15 % spirals / laps


# ---------------------------------------------------------------------------
# summary
# ---------------------------------------------------------------------------

def summarise_boq(result: BOQResult) -> str:
    """Return a formatted text summary of the BOQ."""
    w = 100
    lines = [
        "=" * w,
        "BILL OF QUANTITIES (BOQ)",
        "=" * w,
        f"{'Element':<15} {'Concrete (m3)':>14} {'Steel (kg)':>14} "
        f"{'Conc. Cost':>14} {'Steel Cost':>14} {'Total Cost':>14}",
        "-" * w,
    ]
    for item in result.items:
        lines.append(
            f"{item.element:<15} "
            f"{item.concrete_volume:>14.2f} "
            f"{item.steel_weight:>14.1f} "
            f"{item.concrete_cost:>14,.0f} "
            f"{item.steel_cost:>14,.0f} "
            f"{item.total_cost:>14,.0f}"
        )
    lines.append("-" * w)
    lines.append(
        f"{'TOTAL':<15} "
        f"{result.total_concrete:>14.2f} "
        f"{result.total_steel:>14.1f} "
        f"{'':>14} "
        f"{'':>14} "
        f"{result.total_cost:>14,.0f}"
    )
    lines.append("=" * w)
    lines.append("")
    lines.append(f"Total Concrete Volume : {result.total_concrete:.2f} m3")
    lines.append(f"Total Steel Weight    : {result.total_steel:.1f} kg "
                 f"({result.total_steel / 1000:.2f} tonnes)")
    lines.append(f"Total Project Cost    : INR {result.total_cost:,.0f}")
    lines.append("=" * w)
    return "\n".join(lines)

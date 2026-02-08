"""
Pile capacity analysis module.

This module distributes forces from the pile cap to individual piles and checks pile capacities
according to IRC 78 guidelines.
"""

from dataclasses import dataclass
from typing import Any

from .loads import ForceVector
from .geometry import GeometryResults, PileCoordinate, PileProperties
from .load_combinations import CombinationResult, CombinationResults


@dataclass
class PileLoad:
    """Load on a single pile."""
    pile_index: int
    pile_coord: PileCoordinate
    P: float           # kN -- vertical load
    HL: float          # kN -- horizontal longitudinal per pile
    HT: float          # kN -- horizontal transverse per pile


@dataclass
class PileCombinationResult:
    """Pile loads for a single load combination."""
    combo_name: str
    category: str
    pile_loads: list[PileLoad]
    max_pile_P: float   # kN -- max compression
    min_pile_P: float   # kN -- min (most negative = uplift)


@dataclass
class PileCapacityResult:
    """Complete pile capacity analysis results."""
    num_piles: int
    pile_diameter: float              # m
    pile_capacity_compression: float  # kN
    pile_capacity_uplift: float       # kN (negative)
    sum_x2: float                     # m2
    sum_y2: float                     # m2
    all_pile_results: list[PileCombinationResult]
    max_compression_pile: PileLoad
    max_compression_combo: str
    max_uplift_pile: PileLoad
    max_uplift_combo: str
    compression_util: float
    uplift_util: float
    status: str


def calculate_pile_capacity(
    config: dict[str, Any],
    geometry: GeometryResults,
    combinations: CombinationResults,
) -> PileCapacityResult:
    """
    Calculate pile capacity and distribution of forces to individual piles.

    Uses the pile group formula:
    P_pile_i = P/n ± ML * y_i / sum(y_j^2) ± MT * x_i / sum(x_j^2)

    Args:
        config: Configuration dictionary with foundation parameters
        geometry: Geometry results containing pile coordinates
        combinations: Load combination results

    Returns:
        PileCapacityResult with detailed pile load distribution and capacity checks
    """
    # Extract pile coordinates and properties
    pile_coords = geometry.pile_coordinates
    num_piles = geometry.num_piles
    pile_diameter = geometry.pile_props.diameter

    # Compute sum of squared coordinates (in meters)
    sum_x2 = sum(p.x ** 2 for p in pile_coords)
    sum_y2 = sum(p.y ** 2 for p in pile_coords)

    # Extract pile capacities from config and convert tonnes to kN (using 9.81)
    foundation_config = config.get('foundation', {})
    pile_capacity_compression = foundation_config.get('pile_capacity_top', 550) * 9.81  # kN
    pile_capacity_uplift = foundation_config.get('pile_uplift_capacity', -187) * 9.81   # kN (negative)

    # Process each load combination
    # Per IRC 78, pile capacity is checked using SERVICE loads against
    # safe bearing capacity (which already includes FOS per IS:2911).
    # ULS combinations are still processed for structural pile design
    # (bending, shear) but NOT used for capacity utilization check.
    all_pile_results = []
    global_max_compression = float('-inf')
    global_min_compression = float('inf')  # Most negative = max uplift
    max_compression_pile = None
    max_compression_combo = ""
    max_uplift_pile = None
    max_uplift_combo = ""
    # Track service-load maxima separately for capacity check
    sls_max_compression = float('-inf')
    sls_min_compression = float('inf')
    sls_max_compression_pile = None
    sls_max_compression_combo = ""
    sls_max_uplift_pile = None
    sls_max_uplift_combo = ""

    for combo in combinations.all_combinations:
        fv = combo.forces_pilecap_top  # ForceVector at top of pile cap

        # Calculate loads on each pile
        pile_loads = []
        combo_max_P = float('-inf')
        combo_min_P = float('inf')

        for i, pile_coord in enumerate(pile_coords):
            # Vertical load using pile group formula
            # P_i = P/n + ML * y_i / sum_y2 + MT * x_i / sum_x2
            P_avg = fv.P / num_piles
            P_from_ML = (fv.ML * pile_coord.y / sum_y2) if sum_y2 != 0 else 0.0
            P_from_MT = (fv.MT * pile_coord.x / sum_x2) if sum_x2 != 0 else 0.0
            P_i = P_avg + P_from_ML + P_from_MT

            # Horizontal loads distributed equally
            HL_i = fv.FL / num_piles
            HT_i = fv.FT / num_piles

            pile_load = PileLoad(
                pile_index=i,
                pile_coord=pile_coord,
                P=P_i,
                HL=HL_i,
                HT=HT_i
            )
            pile_loads.append(pile_load)

            # Track max/min for this combination
            if P_i > combo_max_P:
                combo_max_P = P_i
            if P_i < combo_min_P:
                combo_min_P = P_i

            # Track global max compression and max uplift
            if P_i > global_max_compression:
                global_max_compression = P_i
                max_compression_pile = pile_load
                max_compression_combo = combo.name

            if P_i < global_min_compression:
                global_min_compression = P_i
                max_uplift_pile = pile_load
                max_uplift_combo = combo.name

            # Track SLS maxima for capacity check (service loads)
            is_sls = combo.category.startswith("sls_")
            is_uls_seismic = combo.category == "uls_seismic"
            if is_sls or is_uls_seismic:
                if P_i > sls_max_compression:
                    sls_max_compression = P_i
                    sls_max_compression_pile = pile_load
                    sls_max_compression_combo = combo.name
                if P_i < sls_min_compression:
                    sls_min_compression = P_i
                    sls_max_uplift_pile = pile_load
                    sls_max_uplift_combo = combo.name

        # Store results for this combination
        pile_combo_result = PileCombinationResult(
            combo_name=combo.name,
            category=combo.category,
            pile_loads=pile_loads,
            max_pile_P=combo_max_P,
            min_pile_P=combo_min_P
        )
        all_pile_results.append(pile_combo_result)

    # Determine allowable capacities based on IRC 78 Table 5 / Cl. 709.1.6
    # Seismic and wind combinations allow 25% increase in safe capacity.
    # Normal combinations use 100% of capacity.
    # Per IRC 78, capacity check uses SERVICE loads against safe bearing
    # capacity (which already includes FOS per IS:2911).

    compression_util = 0.0
    uplift_util = 0.0

    # Per IRC 78: check each SLS/seismic combination individually and find
    # the one with the worst utilization (highest P/capacity ratio).
    for pile_result in all_pile_results:
        is_sls = pile_result.category.startswith("sls_")
        is_uls_seismic = pile_result.category == "uls_seismic"
        if not (is_sls or is_uls_seismic):
            continue

        # Transient load combos (seismic, wind-leading) get 25% bonus
        is_transient = (
            is_uls_seismic
            or "Wind Leading" in pile_result.combo_name
        )
        bonus = 1.25 if is_transient else 1.0
        allow_comp = pile_capacity_compression * bonus
        allow_upl = pile_capacity_uplift * bonus

        # Check compression
        util_c = pile_result.max_pile_P / allow_comp if allow_comp > 0 else 0.0
        if util_c > compression_util:
            compression_util = util_c
            sls_max_compression = pile_result.max_pile_P
            sls_max_compression_combo = pile_result.combo_name
            # Find the pile with max P in this combo
            for pl in pile_result.pile_loads:
                if abs(pl.P - pile_result.max_pile_P) < 0.1:
                    sls_max_compression_pile = pl
                    break

        # Check uplift
        if allow_upl != 0:
            util_u = pile_result.min_pile_P / allow_upl
        else:
            util_u = 0.0
        if util_u > uplift_util:
            uplift_util = util_u
            sls_min_compression = pile_result.min_pile_P
            sls_max_uplift_combo = pile_result.combo_name
            for pl in pile_result.pile_loads:
                if abs(pl.P - pile_result.min_pile_P) < 0.1:
                    sls_max_uplift_pile = pl
                    break

    # Fall back to global values if no SLS/seismic combos found
    if sls_max_compression_pile is None:
        sls_max_compression_pile = max_compression_pile
        sls_max_compression_combo = max_compression_combo
        sls_max_compression = global_max_compression
        compression_util = global_max_compression / pile_capacity_compression
    if sls_max_uplift_pile is None:
        sls_max_uplift_pile = max_uplift_pile
        sls_max_uplift_combo = max_uplift_combo
        sls_min_compression = global_min_compression
        if pile_capacity_uplift != 0:
            uplift_util = global_min_compression / pile_capacity_uplift

    # Determine status
    if compression_util > 1.0 or uplift_util > 1.0:
        status = "NOT OK"
    else:
        status = "OK"

    return PileCapacityResult(
        num_piles=num_piles,
        pile_diameter=pile_diameter,
        pile_capacity_compression=pile_capacity_compression,
        pile_capacity_uplift=pile_capacity_uplift,
        sum_x2=sum_x2,
        sum_y2=sum_y2,
        all_pile_results=all_pile_results,
        max_compression_pile=max_compression_pile,
        max_compression_combo=max_compression_combo,
        max_uplift_pile=max_uplift_pile,
        max_uplift_combo=max_uplift_combo,
        compression_util=compression_util,
        uplift_util=uplift_util,
        status=status
    )


def summarise_pile_capacity(result: PileCapacityResult) -> str:
    """
    Generate a formatted summary of pile capacity analysis results.

    Args:
        result: PileCapacityResult object

    Returns:
        Formatted string summary
    """
    lines = []
    lines.append("=" * 80)
    lines.append("PILE CAPACITY ANALYSIS")
    lines.append("=" * 80)
    lines.append("")

    # Pile configuration
    lines.append("Pile Configuration:")
    lines.append(f"  Number of piles: {result.num_piles}")
    lines.append(f"  Pile diameter: {result.pile_diameter:.3f} m")
    lines.append(f"  Sum of x^2: {result.sum_x2:.3f} m^2")
    lines.append(f"  Sum of y^2: {result.sum_y2:.3f} m^2")
    lines.append("")

    # Pile capacities
    lines.append("Pile Capacities:")
    lines.append(f"  Compression capacity: {result.pile_capacity_compression:.1f} kN "
                 f"({result.pile_capacity_compression / 9.81:.1f} tonnes)")
    lines.append(f"  Uplift capacity: {result.pile_capacity_uplift:.1f} kN "
                 f"({result.pile_capacity_uplift / 9.81:.1f} tonnes)")
    lines.append("")

    # Maximum compression
    lines.append("Maximum Compression:")
    lines.append(f"  Load combination: {result.max_compression_combo}")
    lines.append(f"  Pile index: {result.max_compression_pile.pile_index}")
    lines.append(f"  Location: x={result.max_compression_pile.pile_coord.x:.3f} m, "
                 f"y={result.max_compression_pile.pile_coord.y:.3f} m")
    lines.append(f"  Vertical load (P): {result.max_compression_pile.P:.1f} kN "
                 f"({result.max_compression_pile.P / 9.81:.1f} tonnes)")
    lines.append(f"  Horizontal load (HL): {result.max_compression_pile.HL:.1f} kN")
    lines.append(f"  Horizontal load (HT): {result.max_compression_pile.HT:.1f} kN")
    lines.append(f"  Utilization: {result.compression_util:.2%}")
    lines.append("")

    # Maximum uplift
    lines.append("Maximum Uplift:")
    lines.append(f"  Load combination: {result.max_uplift_combo}")
    lines.append(f"  Pile index: {result.max_uplift_pile.pile_index}")
    lines.append(f"  Location: x={result.max_uplift_pile.pile_coord.x:.3f} m, "
                 f"y={result.max_uplift_pile.pile_coord.y:.3f} m")
    lines.append(f"  Vertical load (P): {result.max_uplift_pile.P:.1f} kN "
                 f"({result.max_uplift_pile.P / 9.81:.1f} tonnes)")
    lines.append(f"  Horizontal load (HL): {result.max_uplift_pile.HL:.1f} kN")
    lines.append(f"  Horizontal load (HT): {result.max_uplift_pile.HT:.1f} kN")
    lines.append(f"  Utilization: {result.uplift_util:.2%}")
    lines.append("")

    # Overall status
    lines.append("=" * 80)
    lines.append(f"PILE CAPACITY CHECK: {result.status}")
    lines.append("=" * 80)
    lines.append("")

    # Summary table of all combinations
    lines.append("Summary of All Load Combinations:")
    lines.append("-" * 80)
    lines.append(f"{'Combination':<30} {'Category':<12} {'Max P (kN)':<15} {'Min P (kN)':<15}")
    lines.append("-" * 80)

    for pile_result in result.all_pile_results:
        lines.append(f"{pile_result.combo_name:<30} "
                    f"{pile_result.category:<12} "
                    f"{pile_result.max_pile_P:>14.1f} "
                    f"{pile_result.min_pile_P:>14.1f}")

    lines.append("-" * 80)
    lines.append("")

    return "\n".join(lines)

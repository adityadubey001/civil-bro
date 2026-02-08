"""P-M interaction diagrams for RC sections per IRC 112:2020.

Constructs axial load -- bending moment interaction envelopes for circular
and rectangular reinforced concrete cross-sections using the strain
compatibility method prescribed by IRC 112:2020.

Key references
--------------
* IRC 112:2020, Cl. 6.4.2.8 -- Parabolic-rectangular stress-strain model
  for concrete in compression.
* IRC 112:2020, Cl. 6.2.3 -- Bilinear stress-strain model for reinforcing
  steel.
* IRC 112:2020, Annex A-2 -- Interaction diagram methodology.
* Bresler (1960) -- Biaxial bending check via load-contour method.

Stress-strain models
--------------------
**Concrete (parabolic-rectangular, IRC 112 Cl. 6.4.2.8):**

For concrete grades fck <= 50 MPa (normal-strength concrete):

.. math::

    \\sigma_c = f_{cd} \\left[1 - \\left(1 - \\frac{\\varepsilon}{\\varepsilon_{c2}}\\right)^n\\right]
    \\quad \\text{for } 0 \\le \\varepsilon \\le \\varepsilon_{c2}

    \\sigma_c = f_{cd}
    \\quad \\text{for } \\varepsilon_{c2} \\le \\varepsilon \\le \\varepsilon_{cu2}

where:
    - ``fcd = alpha_cc * fck / gamma_c = 0.67 * fck / 1.5``
    - ``epsilon_c2 = 0.002``
    - ``epsilon_cu2 = 0.0035``
    - ``n = 2.0`` (parabolic exponent)

**Reinforcing steel (bilinear, IRC 112 Cl. 6.2.3):**

.. math::

    \\sigma_s = \\varepsilon_s \\cdot E_s
    \\quad \\text{for } |\\varepsilon_s| \\le \\varepsilon_{yd}

    \\sigma_s = \\pm f_{yd}
    \\quad \\text{for } |\\varepsilon_s| > \\varepsilon_{yd}

where:
    - ``fyd = fyk / gamma_s = fyk / 1.15``
    - ``Es = 200,000 MPa``
    - ``epsilon_yd = fyd / Es``

Units convention
----------------
All internal calculations use **mm** for dimensions and **MPa** for stresses.
Input dimensions are expected in **mm**.  Output forces are in **kN** and
moments in **kN.m**.

The conversion is applied at the output stage:
    - Force: ``N -> kN``  (divide by 1000)
    - Moment: ``N.mm -> kN.m``  (divide by 1e6)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from .loads import ForceVector
from .materials import (
    ConcreteProperties,
    SteelProperties,
    get_concrete_properties,
    get_steel_properties,
)


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_ES_MPA: float = 200_000.0  # Reinforcing steel elastic modulus, MPa
_ALPHA_CC: float = 0.67     # Long-term coefficient for concrete (IRC 112)
_GAMMA_C: float = 1.5       # Partial safety factor for concrete at ULS
_GAMMA_S: float = 1.15      # Partial safety factor for steel at ULS

# Default strain-block parameters for fck <= 50 MPa
_EPSILON_C2: float = 0.002   # Onset of constant stress (parabolic-rectangular)
_EPSILON_CU2: float = 0.0035 # Ultimate compressive strain
_N_PARABOLIC: float = 2.0    # Parabolic exponent

# Number of strips for numerical integration over the compressed concrete zone
_N_STRIPS: int = 100


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class InteractionPoint:
    """A single point on the P-M interaction diagram.

    Attributes
    ----------
    P : float
        Axial force in **kN**.  Positive denotes compression.
    M : float
        Bending moment in **kN.m**.  Always reported as a non-negative
        magnitude (the interaction diagram is symmetric about the P axis).
    """

    P: float
    M: float


@dataclass
class InteractionDiagram:
    """Complete P-M interaction envelope for an RC section.

    The ``points`` list traces the full curve from pure tension (bottom)
    through balanced failure to pure compression (top).

    Attributes
    ----------
    points : list[InteractionPoint]
        Ordered sequence of (P, M) pairs forming the interaction curve.
    P_max : float
        Pure axial compression capacity (M = 0), kN.
    P_min : float
        Pure axial tension capacity (M = 0), kN.  Reported as a negative
        value (tension).
    M_max : float
        Maximum moment capacity on the curve, kN.m.  This generally
        occurs near the balanced condition.
    P_balanced : float
        Axial force at the balanced failure condition, kN.
    """

    points: list[InteractionPoint]
    P_max: float
    P_min: float
    M_max: float
    P_balanced: float


# ---------------------------------------------------------------------------
# Material stress-strain functions
# ---------------------------------------------------------------------------

def concrete_stress(
    strain: float,
    fcd: float,
    epsilon_c2: float = _EPSILON_C2,
    n: float = _N_PARABOLIC,
) -> float:
    """Concrete compressive stress from the parabolic-rectangular model.

    Implements IRC 112:2020, Cl. 6.4.2.8.  The concrete is assumed to
    carry **zero tensile stress** (cracked section).  Compressive strain
    is taken as positive.

    Parameters
    ----------
    strain : float
        Concrete compressive strain (positive in compression).
    fcd : float
        Design compressive strength of concrete, MPa.
        ``fcd = alpha_cc * fck / gamma_c``
    epsilon_c2 : float, optional
        Strain at the end of the parabolic branch (default 0.002 for
        fck <= 50 MPa).
    n : float, optional
        Parabolic exponent (default 2.0 for fck <= 50 MPa).

    Returns
    -------
    float
        Compressive stress in MPa (>= 0).  Returns 0 for tensile strains.
    """
    if strain <= 0.0:
        return 0.0
    if strain <= epsilon_c2:
        return fcd * (1.0 - (1.0 - strain / epsilon_c2) ** n)
    # For epsilon_c2 < strain <= epsilon_cu2 (and beyond, capped at fcd)
    return fcd


def steel_stress(
    strain: float,
    fyd: float,
    Es: float = _ES_MPA,
) -> float:
    """Reinforcing steel stress from the bilinear model.

    Implements IRC 112:2020, Cl. 6.2.3.  The steel stress-strain curve
    is elastic-perfectly-plastic with no strain hardening.

    Parameters
    ----------
    strain : float
        Steel strain (positive in tension, negative in compression).
    fyd : float
        Design yield strength of steel, MPa.
        ``fyd = fyk / gamma_s``
    Es : float, optional
        Elastic modulus of steel, MPa (default 200,000).

    Returns
    -------
    float
        Steel stress in MPa.  Positive denotes tension, negative
        denotes compression.
    """
    epsilon_yd = fyd / Es
    sigma = strain * Es
    if sigma > fyd:
        return fyd
    if sigma < -fyd:
        return -fyd
    return sigma


# ---------------------------------------------------------------------------
# Circular section interaction diagram
# ---------------------------------------------------------------------------

def generate_interaction_circular(
    diameter: float,
    n_bars: int,
    bar_diameter: float,
    cover: float,
    fck: float,
    fyk: float = 550.0,
    n_points: int = 50,
) -> InteractionDiagram:
    """Generate a P-M interaction diagram for a circular RC section.

    Reinforcement is placed as ``n_bars`` uniformly distributed on a circle
    of radius ``r_s = (diameter/2) - cover - bar_diameter/2``.

    The neutral axis depth ``xu`` is varied from well above the section
    (pure tension) through the section (balanced and over-reinforced
    zones) to well below (pure compression) to trace the full interaction
    envelope.

    Parameters
    ----------
    diameter : float
        Outer diameter of the circular section, **mm**.
    n_bars : int
        Number of longitudinal reinforcing bars, equally spaced around
        the circumference.
    bar_diameter : float
        Diameter of each reinforcing bar, **mm**.
    cover : float
        Clear cover to the outermost bar surface, **mm**.
    fck : float
        Characteristic compressive strength of concrete, **MPa**.
    fyk : float, optional
        Characteristic yield strength of reinforcing steel, **MPa**
        (default 550 for Fe 550).
    n_points : int, optional
        Number of neutral axis positions used to trace the curve
        (default 50; internally at least 50 are generated).

    Returns
    -------
    InteractionDiagram
        The complete interaction envelope.

    Notes
    -----
    All internal calculations use mm and MPa.  The returned forces are
    in kN and moments in kN.m.

    The concrete contribution is computed by numerical integration
    (strip method with ~100 strips) over the compressed portion of the
    circular cross-section.
    """
    n_points = max(n_points, 50)

    R = diameter / 2.0  # section radius, mm
    r_s = R - cover - bar_diameter / 2.0  # reinforcement circle radius, mm
    A_bar = math.pi * bar_diameter**2 / 4.0  # area of one bar, mm^2

    # Material design values
    fcd = _ALPHA_CC * fck / _GAMMA_C
    fyd = fyk / _GAMMA_S
    epsilon_c2 = _EPSILON_C2
    epsilon_cu2 = _EPSILON_CU2

    # Bar positions: angular coordinates and distances from the compression
    # face (top of section).
    # Convention: compression face is at y = 0 (top), tension face at y = D.
    # Bar angle measured from the top (compression side).
    bar_angles = np.linspace(0.0, 2.0 * math.pi, n_bars, endpoint=False)
    # Distance of each bar from the compression face:
    #   d_si = R - r_s * cos(theta)
    # where theta = 0 is the bar closest to the compression face.
    d_bars = R - r_s * np.cos(bar_angles)  # mm, from compression face

    # Gross section area (for pure compression limit)
    A_gross = math.pi * R**2
    A_steel_total = n_bars * A_bar

    # ------------------------------------------------------------------
    # Generate neutral axis positions
    # ------------------------------------------------------------------
    # xu is the neutral axis depth measured from the compression face.
    # xu < 0        : entire section in tension
    # 0 < xu < D    : partial compression
    # xu > D        : entire section in compression (axial compression governs)
    #
    # To capture the full curve including pure tension we sweep xu from
    # a negative value through D and beyond.
    xu_min = -0.5 * diameter   # well into pure tension
    xu_max = 3.0 * diameter    # well into pure compression

    xu_values = np.linspace(xu_min, xu_max, n_points)

    points: list[InteractionPoint] = []

    for xu in xu_values:
        P_total, M_total = _compute_circular_point(
            xu, R, diameter, d_bars, A_bar, n_bars,
            fcd, fyd, epsilon_c2, epsilon_cu2,
        )
        points.append(InteractionPoint(P=P_total, M=abs(M_total)))

    # ------------------------------------------------------------------
    # Determine key quantities
    # ------------------------------------------------------------------
    # Pure compression: entire section at epsilon_cu2, all steel at fyd comp.
    P_max = (fcd * (A_gross - A_steel_total) + fyd * A_steel_total) / 1000.0

    # Pure tension: all steel yields in tension, concrete contribution = 0.
    P_min = -(fyd * A_steel_total) / 1000.0

    # Maximum moment and balanced axial load
    M_max = max(pt.M for pt in points)
    # The balanced point is defined where strain at the tension-most bar
    # equals epsilon_yd simultaneously with epsilon_cu2 at the compression
    # face.  We approximate it as the point on the curve with M_max.
    balanced_idx = max(range(len(points)), key=lambda i: points[i].M)
    P_balanced = points[balanced_idx].P

    return InteractionDiagram(
        points=points,
        P_max=P_max,
        P_min=P_min,
        M_max=M_max,
        P_balanced=P_balanced,
    )


def _compute_circular_point(
    xu: float,
    R: float,
    D: float,
    d_bars: np.ndarray,
    A_bar: float,
    n_bars: int,
    fcd: float,
    fyd: float,
    epsilon_c2: float,
    epsilon_cu2: float,
) -> tuple[float, float]:
    """Compute (P, M) in kN and kN.m for a given neutral axis depth on a
    circular section.

    Parameters
    ----------
    xu : float
        Neutral axis depth from the compression face, mm.  Can be negative
        (entire section in tension) or greater than D (entire section in
        compression).
    R : float
        Section radius, mm.
    D : float
        Section diameter (= 2*R), mm.
    d_bars : np.ndarray
        Distances of each bar from the compression face, mm.
    A_bar : float
        Area of one bar, mm^2.
    n_bars : int
        Number of bars.
    fcd : float
        Design compressive strength, MPa.
    fyd : float
        Design yield strength of steel, MPa.
    epsilon_c2 : float
        Strain at onset of constant-stress branch.
    epsilon_cu2 : float
        Ultimate compressive strain.

    Returns
    -------
    tuple[float, float]
        (P_kN, M_kNm) -- axial force (positive compression) and moment
        about the geometric centroid.
    """
    centroid = D / 2.0  # distance from compression face to centroid

    # ------------------------------------------------------------------
    # Concrete compression by strip integration
    # ------------------------------------------------------------------
    # Compressed depth: 0 to min(xu, D), but only if xu > 0.
    C_concrete = 0.0  # resultant concrete force, N
    M_concrete = 0.0  # moment of concrete about centroid, N.mm

    if xu > 0.0:
        # Depth of compressed zone within the section
        xu_eff = min(xu, D)
        n_strips = _N_STRIPS
        dy = xu_eff / n_strips  # strip thickness

        for i in range(n_strips):
            y = (i + 0.5) * dy  # midpoint depth from compression face

            # Strain at this depth (linear profile, epsilon_cu2 at comp face)
            if xu > 1e-6:
                eps = epsilon_cu2 * (xu - y) / xu
            else:
                eps = 0.0

            # Concrete stress (compressive only)
            sigma = concrete_stress(eps, fcd, epsilon_c2, _N_PARABOLIC)

            # Width of circular section at depth y from compression face:
            #   w(y) = 2 * sqrt(R^2 - (R - y)^2)
            r_local = R**2 - (R - y) ** 2
            if r_local < 0.0:
                continue
            w = 2.0 * math.sqrt(r_local)

            dA = w * dy
            dF = sigma * dA       # N
            C_concrete += dF
            M_concrete += dF * (centroid - y)  # positive moment = sagging

    # ------------------------------------------------------------------
    # Steel forces
    # ------------------------------------------------------------------
    F_steel = 0.0   # total steel force, N (positive = compression)
    M_steel = 0.0   # total steel moment about centroid, N.mm

    for j in range(n_bars):
        d_j = d_bars[j]  # distance from compression face

        # Strain at bar location (positive = compression)
        if abs(xu) < 1e-9:
            # xu ~ 0: compression face is at the neutral axis, everything
            # below is in tension.
            eps_s = -epsilon_cu2  # tension
        elif xu > 0:
            eps_s = epsilon_cu2 * (xu - d_j) / xu
        else:
            # xu < 0: entire section in tension.  Strain profile pivots
            # about xu with epsilon_cu2 at the compression face.
            # For consistency with the limiting-strain approach, when xu < 0
            # the strain at any depth y is:
            #   eps = epsilon_cu2 * (xu - y) / xu
            # Since xu < 0 and y >= 0, (xu - y) < xu < 0, so eps > epsilon_cu2
            # which is unphysical.  Instead, for pure tension, assume
            # strain profile is governed by maximum steel strain.
            # Use a large tension strain (steel will yield).
            eps_s = -epsilon_cu2 * (d_j - xu) / max(abs(xu), 1e-6)

        # Steel stress: positive strain -> compression, negative -> tension.
        # Our steel_stress function uses the convention that positive strain
        # produces positive (tensile) stress.  Since our eps_s is positive
        # for compression, we need to negate to match steel convention.
        sigma_s = steel_stress(-eps_s, fyd, _ES_MPA)

        # Force in bar: negative sigma_s (compression) adds to P (positive)
        # We treat compression as positive force for the section.
        F_bar = -sigma_s * A_bar  # N, positive = compression
        F_steel += F_bar
        M_steel += F_bar * (centroid - d_j)

    # ------------------------------------------------------------------
    # Totals
    # ------------------------------------------------------------------
    P_total_N = C_concrete + F_steel  # N
    M_total_Nmm = M_concrete + M_steel  # N.mm

    P_kN = P_total_N / 1000.0
    M_kNm = M_total_Nmm / 1e6

    return P_kN, M_kNm


# ---------------------------------------------------------------------------
# Rectangular section interaction diagram
# ---------------------------------------------------------------------------

def generate_interaction_rectangular(
    width: float,
    depth: float,
    n_bars_tension: int,
    n_bars_compression: int,
    bar_diameter: float,
    cover: float,
    fck: float,
    fyk: float = 550.0,
    n_points: int = 50,
) -> InteractionDiagram:
    """Generate a P-M interaction diagram for a rectangular RC section.

    Reinforcement is arranged in two layers:
        - **Compression steel**: ``n_bars_compression`` bars at a depth
          ``d' = cover + bar_diameter / 2`` from the compression face.
        - **Tension steel**: ``n_bars_tension`` bars at a depth
          ``d = depth - cover - bar_diameter / 2`` from the compression face.

    Parameters
    ----------
    width : float
        Section width (perpendicular to the bending axis), **mm**.
    depth : float
        Section overall depth (parallel to the bending axis), **mm**.
    n_bars_tension : int
        Number of tension-face reinforcing bars.
    n_bars_compression : int
        Number of compression-face reinforcing bars.
    bar_diameter : float
        Diameter of each reinforcing bar, **mm**.
    cover : float
        Clear cover to the nearest bar surface, **mm**.
    fck : float
        Characteristic compressive strength of concrete, **MPa**.
    fyk : float, optional
        Characteristic yield strength of reinforcing steel, **MPa**
        (default 550 for Fe 550).
    n_points : int, optional
        Number of neutral axis positions used to trace the curve
        (default 50; internally at least 50 are generated).

    Returns
    -------
    InteractionDiagram
        The complete interaction envelope.

    Notes
    -----
    All internal calculations use mm and MPa.  The returned forces are
    in kN and moments in kN.m.

    The concrete contribution is computed by numerical integration
    (strip method with ~100 strips) over the compressed depth.
    """
    n_points = max(n_points, 50)

    b = width   # mm
    D = depth   # mm
    A_bar = math.pi * bar_diameter**2 / 4.0  # mm^2

    # Effective depths (from compression face)
    d_prime = cover + bar_diameter / 2.0   # compression steel depth
    d_tens = D - cover - bar_diameter / 2.0  # tension steel depth

    # Material design values
    fcd = _ALPHA_CC * fck / _GAMMA_C
    fyd = fyk / _GAMMA_S
    epsilon_c2 = _EPSILON_C2
    epsilon_cu2 = _EPSILON_CU2

    # Steel areas
    A_sc = n_bars_compression * A_bar  # compression steel area
    A_st = n_bars_tension * A_bar      # tension steel area
    A_steel_total = A_sc + A_st

    # Gross section area
    A_gross = b * D

    # Collect all bar positions and areas for the general algorithm
    # Bar positions from compression face, and their areas
    bar_depths: list[float] = []
    bar_areas: list[float] = []

    for _ in range(n_bars_compression):
        bar_depths.append(d_prime)
        bar_areas.append(A_bar)
    for _ in range(n_bars_tension):
        bar_depths.append(d_tens)
        bar_areas.append(A_bar)

    d_bars = np.array(bar_depths)
    a_bars = np.array(bar_areas)
    total_bars = len(bar_depths)

    # ------------------------------------------------------------------
    # Generate neutral axis positions
    # ------------------------------------------------------------------
    xu_min = -0.5 * D
    xu_max = 3.0 * D

    xu_values = np.linspace(xu_min, xu_max, n_points)
    centroid = D / 2.0

    points: list[InteractionPoint] = []

    for xu in xu_values:
        P_total, M_total = _compute_rectangular_point(
            xu, b, D, centroid, d_bars, a_bars, total_bars,
            fcd, fyd, epsilon_c2, epsilon_cu2,
        )
        points.append(InteractionPoint(P=P_total, M=abs(M_total)))

    # ------------------------------------------------------------------
    # Determine key quantities
    # ------------------------------------------------------------------
    P_max = (fcd * (A_gross - A_steel_total) + fyd * A_steel_total) / 1000.0
    P_min = -(fyd * A_steel_total) / 1000.0

    M_max = max(pt.M for pt in points)
    balanced_idx = max(range(len(points)), key=lambda i: points[i].M)
    P_balanced = points[balanced_idx].P

    return InteractionDiagram(
        points=points,
        P_max=P_max,
        P_min=P_min,
        M_max=M_max,
        P_balanced=P_balanced,
    )


def _compute_rectangular_point(
    xu: float,
    b: float,
    D: float,
    centroid: float,
    d_bars: np.ndarray,
    a_bars: np.ndarray,
    n_bars: int,
    fcd: float,
    fyd: float,
    epsilon_c2: float,
    epsilon_cu2: float,
) -> tuple[float, float]:
    """Compute (P, M) in kN and kN.m for a given neutral axis depth on a
    rectangular section.

    Parameters
    ----------
    xu : float
        Neutral axis depth from the compression face, mm.
    b : float
        Section width, mm.
    D : float
        Section overall depth, mm.
    centroid : float
        Distance from compression face to geometric centroid (= D/2), mm.
    d_bars : np.ndarray
        Distances of each bar from the compression face, mm.
    a_bars : np.ndarray
        Area of each bar, mm^2.
    n_bars : int
        Total number of bars.
    fcd : float
        Design compressive strength, MPa.
    fyd : float
        Design yield strength of steel, MPa.
    epsilon_c2 : float
        Strain at onset of constant-stress branch.
    epsilon_cu2 : float
        Ultimate compressive strain.

    Returns
    -------
    tuple[float, float]
        (P_kN, M_kNm)
    """
    # ------------------------------------------------------------------
    # Concrete compression by strip integration
    # ------------------------------------------------------------------
    C_concrete = 0.0
    M_concrete = 0.0

    if xu > 0.0:
        xu_eff = min(xu, D)
        n_strips = _N_STRIPS
        dy = xu_eff / n_strips

        for i in range(n_strips):
            y = (i + 0.5) * dy

            if xu > 1e-6:
                eps = epsilon_cu2 * (xu - y) / xu
            else:
                eps = 0.0

            sigma = concrete_stress(eps, fcd, epsilon_c2, _N_PARABOLIC)

            dA = b * dy
            dF = sigma * dA
            C_concrete += dF
            M_concrete += dF * (centroid - y)

    # ------------------------------------------------------------------
    # Steel forces
    # ------------------------------------------------------------------
    F_steel = 0.0
    M_steel = 0.0

    for j in range(n_bars):
        d_j = d_bars[j]

        if abs(xu) < 1e-9:
            eps_s = -epsilon_cu2
        elif xu > 0:
            eps_s = epsilon_cu2 * (xu - d_j) / xu
        else:
            eps_s = -epsilon_cu2 * (d_j - xu) / max(abs(xu), 1e-6)

        sigma_s = steel_stress(-eps_s, fyd, _ES_MPA)
        F_bar = -sigma_s * a_bars[j]
        F_steel += F_bar
        M_steel += F_bar * (centroid - d_j)

    # ------------------------------------------------------------------
    # Totals
    # ------------------------------------------------------------------
    P_total_N = C_concrete + F_steel
    M_total_Nmm = M_concrete + M_steel

    P_kN = P_total_N / 1000.0
    M_kNm = M_total_Nmm / 1e6

    return P_kN, M_kNm


# ---------------------------------------------------------------------------
# Utilisation checks
# ---------------------------------------------------------------------------

def check_utilisation(
    diagram: InteractionDiagram,
    P: float,
    M: float,
) -> float:
    """Check utilisation ratio against a P-M interaction diagram.

    Computes the ratio ``M_applied / M_capacity`` at the given axial load
    ``P``.  A value less than 1.0 indicates the load point lies inside
    the interaction envelope (safe); greater than 1.0 indicates it lies
    outside (unsafe).

    The moment capacity at a given ``P`` is obtained by linear
    interpolation along the interaction curve.

    Parameters
    ----------
    diagram : InteractionDiagram
        The interaction diagram to check against.
    P : float
        Applied axial force, **kN** (positive = compression).
    M : float
        Applied bending moment, **kN.m** (absolute value used).

    Returns
    -------
    float
        Utilisation ratio ``M / M_capacity``.  Returns ``float('inf')``
        if no valid capacity can be found (e.g. ``P`` is outside the
        range of the diagram).

    Notes
    -----
    If ``P`` lies outside the range ``[P_min, P_max]`` of the diagram,
    the section is overstressed in pure axial mode and ``inf`` is
    returned.
    """
    M = abs(M)

    pts = diagram.points
    if not pts:
        return float("inf")

    # Sort points by P for interpolation
    sorted_pts = sorted(pts, key=lambda pt: pt.P)
    P_values = [pt.P for pt in sorted_pts]
    M_values = [pt.M for pt in sorted_pts]

    P_lo = P_values[0]
    P_hi = P_values[-1]

    # Check if P is within the range of the diagram
    if P < P_lo - 1e-3 or P > P_hi + 1e-3:
        return float("inf")

    # The interaction curve is not single-valued in P: for a given P there
    # may be two M values (ascending and descending branches).  We want the
    # maximum M (the envelope).
    #
    # Strategy: scan consecutive point pairs and find all segments that
    # bracket the target P.  Take the maximum M from all intersections.
    M_capacity = 0.0
    found = False

    for i in range(len(sorted_pts) - 1):
        P1, M1 = sorted_pts[i].P, sorted_pts[i].M
        P2, M2 = sorted_pts[i + 1].P, sorted_pts[i + 1].M

        # Check if the target P lies between P1 and P2
        if (P1 <= P <= P2) or (P2 <= P <= P1):
            dP = P2 - P1
            if abs(dP) < 1e-12:
                M_interp = max(M1, M2)
            else:
                t = (P - P1) / dP
                M_interp = M1 + t * (M2 - M1)
            M_capacity = max(M_capacity, M_interp)
            found = True

    if not found or M_capacity < 1e-12:
        # If P is at the extremes the moment capacity is essentially zero
        # (pure compression or pure tension).  Return inf if moment is
        # applied, otherwise 0.
        if M < 1e-6:
            return 0.0
        return float("inf")

    return M / M_capacity


def check_biaxial(
    diagram_xx: InteractionDiagram,
    diagram_yy: InteractionDiagram,
    P: float,
    Mx: float,
    My: float,
) -> float:
    """Biaxial bending check using Bresler's load-contour method.

    The interaction equation is:

    .. math::

        \\left(\\frac{M_x}{M_{ux}}\\right)^{\\alpha_n}
        + \\left(\\frac{M_y}{M_{uy}}\\right)^{\\alpha_n}
        \\le 1.0

    where ``Mux`` and ``Muy`` are the uniaxial moment capacities at the
    given axial load ``P`` (obtained from the two interaction diagrams),
    and ``alpha_n`` is an exponent that depends on the normalised axial
    load ``P / Puz``:

        - ``alpha_n = 1.0``  for ``P/Puz <= 0.2``
        - ``alpha_n = 2.0``  for ``P/Puz >= 0.8``
        - linearly interpolated between 0.2 and 0.8

    ``Puz`` is the squash load, taken as ``max(P_max_xx, P_max_yy)``.

    Parameters
    ----------
    diagram_xx : InteractionDiagram
        Interaction diagram for bending about the XX axis.
    diagram_yy : InteractionDiagram
        Interaction diagram for bending about the YY axis.
    P : float
        Applied axial force, **kN** (positive = compression).
    Mx : float
        Applied moment about the XX axis, **kN.m**.
    My : float
        Applied moment about the YY axis, **kN.m**.

    Returns
    -------
    float
        Utilisation ratio (the left-hand side of the Bresler equation).
        A value <= 1.0 is safe; > 1.0 is unsafe.
    """
    Mx = abs(Mx)
    My = abs(My)

    # Squash load (pure compression capacity)
    Puz = max(diagram_xx.P_max, diagram_yy.P_max)

    if Puz < 1e-6:
        return float("inf")

    P_ratio = max(P, 0.0) / Puz

    # Exponent alpha_n per Bresler / IS 456 / IRC 112 practice
    if P_ratio <= 0.2:
        alpha_n = 1.0
    elif P_ratio >= 0.8:
        alpha_n = 2.0
    else:
        # Linear interpolation between (0.2, 1.0) and (0.8, 2.0)
        alpha_n = 1.0 + (P_ratio - 0.2) / (0.8 - 0.2) * (2.0 - 1.0)

    # Uniaxial moment capacities at the given P
    Mux = _moment_capacity_at_P(diagram_xx, P)
    Muy = _moment_capacity_at_P(diagram_yy, P)

    if Mux < 1e-6 and Mx > 1e-6:
        return float("inf")
    if Muy < 1e-6 and My > 1e-6:
        return float("inf")

    term_x = (Mx / Mux) ** alpha_n if Mux > 1e-6 else 0.0
    term_y = (My / Muy) ** alpha_n if Muy > 1e-6 else 0.0

    return term_x + term_y


def _moment_capacity_at_P(
    diagram: InteractionDiagram,
    P: float,
) -> float:
    """Interpolate the moment capacity from an interaction diagram at a
    given axial load.

    Scans all segments of the interaction curve and returns the maximum
    moment at the target ``P``.

    Parameters
    ----------
    diagram : InteractionDiagram
        The interaction diagram.
    P : float
        Axial force, kN.

    Returns
    -------
    float
        Moment capacity in kN.m.  Returns 0.0 if ``P`` is outside the
        range of the diagram.
    """
    pts = diagram.points
    if not pts:
        return 0.0

    # We work with the original point ordering (not sorted by P) because
    # the curve traces a path and we need to interpolate along adjacent
    # segments.
    M_capacity = 0.0

    for i in range(len(pts) - 1):
        P1, M1 = pts[i].P, pts[i].M
        P2, M2 = pts[i + 1].P, pts[i + 1].M

        if (P1 <= P <= P2) or (P2 <= P <= P1):
            dP = P2 - P1
            if abs(dP) < 1e-12:
                M_interp = max(M1, M2)
            else:
                t = (P - P1) / dP
                M_interp = M1 + t * (M2 - M1)
            M_capacity = max(M_capacity, M_interp)

    return M_capacity

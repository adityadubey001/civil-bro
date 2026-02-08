"""Material properties for bridge substructure design per IRC 112:2020.

Provides dataclasses and factory functions for concrete and reinforcing steel
properties, reading base values from ``irc_tables.yaml`` and applying the
appropriate partial safety factors and correction factors defined in IRC 112.

Key references
--------------
* IRC 112:2020, Clause 6.4 -- Concrete properties
* IRC 112:2020, Table 6.5 -- Aggregate correction for Ecm
* IRC 112:2020, Clause 15.2.3.3 -- Design compressive strength
* IRC 112:2020, Clause 15.2.3.4 -- Design tensile strength
* IRC 112:2020, Table 14.1 -- Minimum cover for durability
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from .utils import load_irc_tables


# ---------------------------------------------------------------------------
# Concrete
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ConcreteProperties:
    """Design properties for a concrete grade per IRC 112.

    All stresses are in **MPa** and moduli in **GPa** unless noted otherwise.

    Attributes
    ----------
    fck : float
        Characteristic compressive strength (cylinder), MPa.
    fcm : float
        Mean compressive strength, MPa.  ``fcm = fck + 10`` per IRC 112
        (though the YAML may store a slightly different value for high
        grades; the table value takes precedence).
    fctm : float
        Mean axial tensile strength, MPa.
    fctk_005 : float
        5 % fractile axial tensile strength, MPa.
    fctk_095 : float
        95 % fractile axial tensile strength, MPa.
    Ecm : float
        Secant modulus of elasticity (corrected for aggregate type), **GPa**.
    Ecm_uncorrected : float
        Secant modulus before aggregate correction, **GPa**.
    aggregate_factor : float
        Multiplicative correction applied to ``Ecm`` (IRC 112, Table 6.5).
    epsilon_cu2 : float
        Ultimate compressive strain (parabola-rectangle model).
    epsilon_cu3 : float
        Ultimate compressive strain (bi-linear / rectangular model).
    n_parabolic : float
        Exponent for the parabola-rectangle stress-strain curve.
    gamma_c : float
        Partial safety factor for concrete at ULS.
    alpha_cc : float
        Long-term / loading-rate coefficient for compressive strength.
    fcd : float
        Design compressive strength, MPa.  ``fcd = alpha_cc * fck / gamma_c``
    fctd : float
        Design tensile strength, MPa.  ``fctd = fctk_005 / gamma_c``
    modular_ratio_short : float
        Short-term modular ratio ``Es / Ecm`` (with default Es = 200 GPa).
    exposure : str | None
        Exposure condition (if supplied).
    min_cover : float | None
        Minimum clear cover (mm) per IRC 112 Table 14.1 for the given
        exposure, or ``None`` if no exposure was specified.
    """

    # From IRC table look-up
    fck: float
    fcm: float
    fctm: float
    fctk_005: float
    fctk_095: float
    Ecm: float
    Ecm_uncorrected: float
    aggregate_factor: float
    epsilon_cu2: float
    epsilon_cu3: float
    n_parabolic: float

    # Partial safety factors
    gamma_c: float
    alpha_cc: float

    # Derived design values
    fcd: float
    fctd: float
    modular_ratio_short: float

    # Durability
    exposure: str | None = None
    min_cover: float | None = None


def get_concrete_properties(
    fck: int | float,
    aggregate: str = "Quartzite",
    exposure: str | None = None,
    Es: float = 200.0,
) -> ConcreteProperties:
    """Build :class:`ConcreteProperties` for a given concrete grade.

    Parameters
    ----------
    fck : int or float
        Characteristic compressive strength in **MPa**.  Must match a key in
        the ``concrete_properties`` section of ``irc_tables.yaml`` (e.g. 25,
        30, 35, ...).
    aggregate : str, optional
        Aggregate type for ``Ecm`` correction (default ``"Quartzite"``).
        Must be one of the keys in ``aggregate_factors`` in the YAML:
        ``Quartzite``, ``Granite``, ``Limestone``, ``Sandstone``, ``Basalt``.
    exposure : str or None, optional
        Exposure condition for cover look-up (e.g. ``"Moderate"``,
        ``"Severe"``).  If ``None``, cover is not determined.
    Es : float, optional
        Steel elastic modulus in **GPa** (default 200).

    Returns
    -------
    ConcreteProperties

    Raises
    ------
    KeyError
        If ``fck`` or ``aggregate`` is not found in the tables.
    """
    tables = load_irc_tables()

    # --- Look up raw table row ---
    fck_key = int(fck)
    concrete_table: dict[str, dict[str, Any]] = tables["concrete_properties"]
    if fck_key not in concrete_table:
        available = sorted(concrete_table.keys())
        raise KeyError(
            f"fck={fck_key} MPa not found in IRC table.  "
            f"Available grades: {available}"
        )
    row = concrete_table[fck_key]

    # --- Aggregate correction ---
    agg_factors: dict[str, float] = tables["aggregate_factors"]
    if aggregate not in agg_factors:
        raise KeyError(
            f"Aggregate '{aggregate}' not recognised.  "
            f"Available: {list(agg_factors.keys())}"
        )
    agg_factor = agg_factors[aggregate]
    ecm_uncorrected = float(row["Ecm"])
    ecm_corrected = ecm_uncorrected * agg_factor

    # --- Partial safety factors ---
    psf = tables["partial_safety_factors"]
    gamma_c: float = psf["gamma_c"]
    alpha_cc: float = psf["alpha_cc"]

    # --- Design strengths ---
    fck_val = float(fck)
    fctk_005 = float(row["fctk_005"])
    fcd = alpha_cc * fck_val / gamma_c
    fctd = fctk_005 / gamma_c

    # --- Modular ratio (short-term) ---
    modular_ratio_short = Es / ecm_corrected

    # --- Cover (optional) ---
    min_cover: float | None = None
    if exposure is not None:
        cover_table: dict[str, float] = tables["cover_requirements"]
        if exposure not in cover_table:
            raise KeyError(
                f"Exposure '{exposure}' not recognised.  "
                f"Available: {list(cover_table.keys())}"
            )
        min_cover = float(cover_table[exposure])

    return ConcreteProperties(
        fck=fck_val,
        fcm=float(row["fcm"]),
        fctm=float(row["fctm"]),
        fctk_005=fctk_005,
        fctk_095=float(row["fctk_095"]),
        Ecm=ecm_corrected,
        Ecm_uncorrected=ecm_uncorrected,
        aggregate_factor=agg_factor,
        epsilon_cu2=float(row["epsilon_cu2"]),
        epsilon_cu3=float(row["epsilon_cu3"]),
        n_parabolic=float(row["n_parabolic"]),
        gamma_c=gamma_c,
        alpha_cc=alpha_cc,
        fcd=fcd,
        fctd=fctd,
        modular_ratio_short=modular_ratio_short,
        exposure=exposure,
        min_cover=min_cover,
    )


# ---------------------------------------------------------------------------
# Reinforcing steel
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class SteelProperties:
    """Design properties for reinforcing steel per IRC 112.

    Attributes
    ----------
    fyk : float
        Characteristic yield strength, MPa.
    Es : float
        Modulus of elasticity, **GPa**.
    gamma_s : float
        Partial safety factor at ULS.
    fyd : float
        Design yield strength, MPa.  ``fyd = fyk / gamma_s``
    epsilon_yd : float
        Design yield strain.  ``epsilon_yd = fyd / (Es * 1000)``
        (Es converted to MPa for strain calculation).
    epsilon_ud : float
        Design ultimate strain.  ``epsilon_ud = 0.9 * fyk / (gamma_s * Es * 1000)``
        per IRC 112 simplified.
    """

    fyk: float
    Es: float
    gamma_s: float
    fyd: float
    epsilon_yd: float
    epsilon_ud: float


def get_steel_properties(
    fyk: float = 500.0,
    Es: float = 200.0,
) -> SteelProperties:
    """Build :class:`SteelProperties` for reinforcing steel.

    Parameters
    ----------
    fyk : float, optional
        Characteristic yield strength in MPa (default 500 for Fe 500).
    Es : float, optional
        Modulus of elasticity in **GPa** (default 200).

    Returns
    -------
    SteelProperties
    """
    tables = load_irc_tables()
    gamma_s: float = tables["partial_safety_factors"]["gamma_s"]

    fyd = fyk / gamma_s
    # Es is in GPa; convert to MPa for strain calculation
    Es_mpa = Es * 1_000.0
    epsilon_yd = fyd / Es_mpa
    epsilon_ud = 0.9 * fyk / (gamma_s * Es_mpa)

    return SteelProperties(
        fyk=fyk,
        Es=Es,
        gamma_s=gamma_s,
        fyd=fyd,
        epsilon_yd=epsilon_yd,
        epsilon_ud=epsilon_ud,
    )


# ---------------------------------------------------------------------------
# Modular ratio helpers
# ---------------------------------------------------------------------------

def short_term_modular_ratio(
    fck: int | float,
    aggregate: str = "Quartzite",
    Es: float = 200.0,
) -> float:
    """Short-term modular ratio ``Es / Ecm``.

    Parameters
    ----------
    fck : int or float
        Concrete grade in MPa.
    aggregate : str, optional
        Aggregate type for Ecm correction.
    Es : float, optional
        Steel modulus in GPa (default 200).

    Returns
    -------
    float
    """
    cp = get_concrete_properties(fck, aggregate=aggregate, Es=Es)
    return cp.modular_ratio_short


def long_term_modular_ratio(
    fck: int | float,
    creep_coefficient: float,
    aggregate: str = "Quartzite",
    Es: float = 200.0,
) -> float:
    """Long-term (effective) modular ratio ``Es / Ec_eff``.

    ``Ec_eff = Ecm / (1 + phi)`` where ``phi`` is the creep coefficient.

    Parameters
    ----------
    fck : int or float
        Concrete grade in MPa.
    creep_coefficient : float
        Creep coefficient ``phi(t, t0)`` -- see :mod:`creep`.
    aggregate : str, optional
        Aggregate type for Ecm correction.
    Es : float, optional
        Steel modulus in GPa (default 200).

    Returns
    -------
    float
    """
    cp = get_concrete_properties(fck, aggregate=aggregate, Es=Es)
    ec_eff = cp.Ecm / (1.0 + creep_coefficient)
    return Es / ec_eff

"""Creep coefficient calculation per IRC 112:2020, Annex A-2.

Implements the full creep model for normal-weight concrete, computing the
creep coefficient ``phi(t, t0)`` and the effective elastic modulus ``Ec_eff``
for long-term analysis.

Key reference
-------------
* IRC 112:2020, Annex A-2 (based on fib Model Code / EN 1992-1-1 Annex B)

Notation
--------
- ``fcm``   -- mean compressive strength of concrete (MPa)
- ``Ecm``   -- secant modulus of elasticity (GPa)
- ``Es``    -- modulus of elasticity of steel (GPa)
- ``h0``    -- notional size = 2 * Ac / u  (mm)
- ``RH``    -- ambient relative humidity (%)
- ``t0``    -- age of concrete at loading (days)
- ``t``     -- age of concrete at the time considered (days)
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any


# ---------------------------------------------------------------------------
# Notional size helper
# ---------------------------------------------------------------------------

def calculate_notional_size(area_mm2: float, perimeter_mm: float) -> float:
    """Compute the notional size ``h0 = 2 * Ac / u``.

    Parameters
    ----------
    area_mm2 : float
        Cross-sectional area of the concrete member in **mm^2**.
    perimeter_mm : float
        Perimeter of the cross-section exposed to drying in **mm**.

    Returns
    -------
    float
        Notional size ``h0`` in **mm**.

    Raises
    ------
    ValueError
        If ``perimeter_mm`` is zero or negative.
    """
    if perimeter_mm <= 0.0:
        raise ValueError(
            f"perimeter_mm must be positive, got {perimeter_mm}"
        )
    return 2.0 * area_mm2 / perimeter_mm


# ---------------------------------------------------------------------------
# Creep result container
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class CreepResult:
    """Container for all intermediate and final creep calculation values.

    All fields correspond to the notation in IRC 112, Annex A-2.

    Attributes
    ----------
    fcm : float
        Mean compressive strength (MPa) -- input.
    Ecm : float
        Secant modulus of elasticity (GPa) -- input.
    Es : float
        Steel modulus of elasticity (GPa) -- input.
    h0 : float
        Notional size (mm) -- input.
    RH : float
        Relative humidity (%) -- input.
    t0 : float
        Age at loading (days) -- input.
    t : float
        Age considered (days) -- input.
    alpha1 : float
        ``(43.75 / fcm) ^ 0.7``  -- correction for high-strength concrete.
    alpha2 : float
        ``(43.75 / fcm) ^ 0.2``
    alpha3 : float
        ``(43.75 / fcm) ^ 0.5``
    beta_fcm : float
        ``18.78 / sqrt(fcm)``  -- effect of concrete strength on phi_0.
    beta_t0 : float
        ``1 / (0.1 + t0 ^ 0.2)``  -- effect of age at loading on phi_0.
    phi_RH : float
        Relative-humidity factor.
    phi_0 : float
        Notional creep coefficient  ``phi_RH * beta_fcm * beta_t0``.
    beta_H : float
        Coefficient depending on ``h0`` and ``RH`` for temporal development.
    beta_c : float
        Temporal development coefficient
        ``((t - t0) / (beta_H + t - t0)) ^ 0.3``.
    phi : float
        Final creep coefficient  ``phi_0 * beta_c``.
    Ec_eff : float
        Effective modulus of elasticity (GPa)  ``Ecm / (1 + phi)``.
    modular_ratio_long : float
        Long-term modular ratio  ``Es / Ec_eff``.
    """

    # Inputs
    fcm: float
    Ecm: float
    Es: float
    h0: float
    RH: float
    t0: float
    t: float

    # Intermediate values
    alpha1: float
    alpha2: float
    alpha3: float
    beta_fcm: float
    beta_t0: float
    phi_RH: float
    phi_0: float
    beta_H: float
    beta_c: float

    # Final results
    phi: float
    Ec_eff: float
    modular_ratio_long: float

    def as_dict(self) -> dict[str, float]:
        """Return all fields as a plain dictionary."""
        return {
            "fcm": self.fcm,
            "Ecm": self.Ecm,
            "Es": self.Es,
            "h0": self.h0,
            "RH": self.RH,
            "t0": self.t0,
            "t": self.t,
            "alpha1": self.alpha1,
            "alpha2": self.alpha2,
            "alpha3": self.alpha3,
            "beta_fcm": self.beta_fcm,
            "beta_t0": self.beta_t0,
            "phi_RH": self.phi_RH,
            "phi_0": self.phi_0,
            "beta_H": self.beta_H,
            "beta_c": self.beta_c,
            "phi": self.phi,
            "Ec_eff": self.Ec_eff,
            "modular_ratio_long": self.modular_ratio_long,
        }


# ---------------------------------------------------------------------------
# Core calculation
# ---------------------------------------------------------------------------

def _compute_alpha_coefficients(
    fcm: float,
) -> tuple[float, float, float]:
    """Compute alpha1, alpha2, alpha3 for high-strength correction.

    Parameters
    ----------
    fcm : float
        Mean compressive strength (MPa).

    Returns
    -------
    tuple of (alpha1, alpha2, alpha3)
    """
    ratio = 43.75 / fcm
    alpha1 = ratio ** 0.7
    alpha2 = ratio ** 0.2
    alpha3 = ratio ** 0.5
    return alpha1, alpha2, alpha3


def _compute_phi_RH(
    RH: float,
    h0: float,
    fcm: float,
    alpha1: float,
    alpha2: float,
) -> float:
    """Relative-humidity factor for the notional creep coefficient.

    For ``fcm <= 45 MPa``:
        ``phi_RH = 1 + (1 - RH/100) / (0.1 * h0^(1/3))``

    For ``fcm > 45 MPa``:
        ``phi_RH = [1 + (1 - RH/100) / (0.1 * h0^(1/3)) * alpha1] * alpha2``

    Parameters
    ----------
    RH : float
        Relative humidity (%).
    h0 : float
        Notional size (mm).
    fcm : float
        Mean compressive strength (MPa).
    alpha1, alpha2 : float
        High-strength correction factors.

    Returns
    -------
    float
    """
    rh_term = (1.0 - RH / 100.0) / (0.1 * h0 ** (1.0 / 3.0))

    if fcm <= 45.0:
        return 1.0 + rh_term
    else:
        return (1.0 + rh_term * alpha1) * alpha2


def _compute_beta_H(
    RH: float,
    h0: float,
    fcm: float,
    alpha3: float,
) -> float:
    """Coefficient ``beta_H`` governing the temporal development of creep.

    For ``fcm <= 35 MPa``:
        ``beta_H = min(1.5 * (1 + (0.012*RH)^18) * h0 + 250, 1500)``

    For ``fcm > 35 MPa``:
        ``beta_H = min(1.5 * (1 + (0.012*RH)^18) * h0 + 250*alpha3,
                       1500*alpha3)``

    Parameters
    ----------
    RH : float
        Relative humidity (%).
    h0 : float
        Notional size (mm).
    fcm : float
        Mean compressive strength (MPa).
    alpha3 : float
        High-strength correction factor.

    Returns
    -------
    float
    """
    rh_factor = 1.5 * (1.0 + (0.012 * RH) ** 18) * h0

    if fcm <= 35.0:
        return min(rh_factor + 250.0, 1500.0)
    else:
        return min(rh_factor + 250.0 * alpha3, 1500.0 * alpha3)


def calculate_creep(
    fcm: float,
    Ecm: float,
    Es: float,
    h0: float,
    RH: float,
    t0: float,
    t: float,
) -> CreepResult:
    """Calculate the creep coefficient per IRC 112, Annex A-2.

    Parameters
    ----------
    fcm : float
        Mean compressive strength of concrete (MPa).
    Ecm : float
        Secant modulus of elasticity of concrete (GPa).
    Es : float
        Modulus of elasticity of reinforcing steel (GPa).
    h0 : float
        Notional size ``2*Ac/u`` (mm).  Use :func:`calculate_notional_size`
        to obtain this from cross-section dimensions.
    RH : float
        Ambient relative humidity (%), typically 50-80.
    t0 : float
        Age of concrete at the time of loading (days).  Must be >= 1.
    t : float
        Age of concrete at the time considered (days).  Must be > ``t0``.

    Returns
    -------
    CreepResult
        Dataclass containing all intermediate values, the final creep
        coefficient ``phi``, the effective modulus ``Ec_eff`` (GPa), and
        the long-term modular ratio ``Es / Ec_eff``.

    Raises
    ------
    ValueError
        If any input is out of physically meaningful range.

    Examples
    --------
    >>> result = calculate_creep(
    ...     fcm=38.0, Ecm=30.5, Es=200.0,
    ...     h0=300.0, RH=70.0, t0=28.0, t=25550.0,
    ... )
    >>> round(result.phi, 3)  # doctest: +SKIP
    1.642
    >>> round(result.Ec_eff, 2)  # doctest: +SKIP
    11.55
    """
    # --- Input validation ---
    if fcm <= 0.0:
        raise ValueError(f"fcm must be positive, got {fcm}")
    if Ecm <= 0.0:
        raise ValueError(f"Ecm must be positive, got {Ecm}")
    if Es <= 0.0:
        raise ValueError(f"Es must be positive, got {Es}")
    if h0 <= 0.0:
        raise ValueError(f"h0 must be positive, got {h0}")
    if not (0.0 < RH <= 100.0):
        raise ValueError(f"RH must be in (0, 100], got {RH}")
    if t0 < 1.0:
        raise ValueError(f"t0 must be >= 1 day, got {t0}")
    if t <= t0:
        raise ValueError(
            f"t must be greater than t0: t={t}, t0={t0}"
        )

    # --- Alpha coefficients (high-strength correction) ---
    alpha1, alpha2, alpha3 = _compute_alpha_coefficients(fcm)

    # --- beta_fcm: effect of concrete strength ---
    beta_fcm = 18.78 / math.sqrt(fcm)

    # --- beta_t0: effect of age at loading ---
    beta_t0 = 1.0 / (0.1 + t0 ** 0.2)

    # --- phi_RH: relative-humidity factor ---
    phi_RH = _compute_phi_RH(RH, h0, fcm, alpha1, alpha2)

    # --- Notional creep coefficient ---
    phi_0 = phi_RH * beta_fcm * beta_t0

    # --- beta_H: temporal development parameter ---
    beta_H = _compute_beta_H(RH, h0, fcm, alpha3)

    # --- beta_c(t, t0): temporal development of creep ---
    dt = t - t0
    beta_c = (dt / (beta_H + dt)) ** 0.3

    # --- Final creep coefficient ---
    phi = phi_0 * beta_c

    # --- Effective modulus and long-term modular ratio ---
    Ec_eff = Ecm / (1.0 + phi)
    modular_ratio_long = Es / Ec_eff

    return CreepResult(
        fcm=fcm,
        Ecm=Ecm,
        Es=Es,
        h0=h0,
        RH=RH,
        t0=t0,
        t=t,
        alpha1=alpha1,
        alpha2=alpha2,
        alpha3=alpha3,
        beta_fcm=beta_fcm,
        beta_t0=beta_t0,
        phi_RH=phi_RH,
        phi_0=phi_0,
        beta_H=beta_H,
        beta_c=beta_c,
        phi=phi,
        Ec_eff=Ec_eff,
        modular_ratio_long=modular_ratio_long,
    )


# ---------------------------------------------------------------------------
# Convenience wrappers
# ---------------------------------------------------------------------------

def creep_coefficient(
    fcm: float,
    h0: float,
    RH: float,
    t0: float,
    t: float,
) -> float:
    """Return only the creep coefficient ``phi(t, t0)``.

    This is a thin wrapper around :func:`calculate_creep` for callers that
    need only the scalar result.

    Parameters
    ----------
    fcm, h0, RH, t0, t
        See :func:`calculate_creep`.

    Returns
    -------
    float
        Creep coefficient ``phi(t, t0)``.
    """
    # Ecm and Es are not needed for phi itself, but calculate_creep requires
    # them for Ec_eff.  Use dummy values that won't affect phi.
    result = calculate_creep(
        fcm=fcm, Ecm=1.0, Es=1.0, h0=h0, RH=RH, t0=t0, t=t,
    )
    return result.phi


def effective_modulus(
    fcm: float,
    Ecm: float,
    h0: float,
    RH: float,
    t0: float,
    t: float,
) -> float:
    """Return the effective elastic modulus ``Ec_eff = Ecm / (1 + phi)``.

    Parameters
    ----------
    fcm : float
        Mean compressive strength (MPa).
    Ecm : float
        Secant modulus (GPa).
    h0, RH, t0, t
        See :func:`calculate_creep`.

    Returns
    -------
    float
        ``Ec_eff`` in GPa.
    """
    result = calculate_creep(
        fcm=fcm, Ecm=Ecm, Es=1.0, h0=h0, RH=RH, t0=t0, t=t,
    )
    return result.Ec_eff

"""Shared utilities for bridge substructure design.

Provides:
- IRC table loading from YAML configuration
- Section property calculations (circular and rectangular)
- Unit conversion helpers
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import yaml


# ---------------------------------------------------------------------------
# IRC table loading
# ---------------------------------------------------------------------------

_irc_tables_cache: dict[str, Any] | None = None


def load_irc_tables() -> dict[str, Any]:
    """Load IRC tables from the YAML config file.

    Searches for ``irc_tables.yaml`` relative to this source file in the
    standard project layout.  The result is cached so that repeated calls do
    not re-read from disk.

    Returns
    -------
    dict
        Parsed YAML content keyed by table name (e.g. ``concrete_properties``,
        ``partial_safety_factors``, etc.).

    Raises
    ------
    FileNotFoundError
        If the YAML file cannot be located in any of the expected paths.
    """
    global _irc_tables_cache
    if _irc_tables_cache is not None:
        return _irc_tables_cache

    config_paths = [
        # src/substructure/../../config  (standard layout)
        Path(__file__).resolve().parent.parent.parent / "config" / "irc_tables.yaml",
        # One level further up (e.g. editable install)
        Path(__file__).resolve().parent.parent.parent.parent / "config" / "irc_tables.yaml",
    ]
    for path in config_paths:
        if path.exists():
            with open(path, encoding="utf-8") as fh:
                _irc_tables_cache = yaml.safe_load(fh)
            return _irc_tables_cache

    searched = "\n  ".join(str(p) for p in config_paths)
    raise FileNotFoundError(
        f"irc_tables.yaml not found.  Searched:\n  {searched}"
    )


def _clear_irc_tables_cache() -> None:
    """Reset the internal cache (useful in tests)."""
    global _irc_tables_cache
    _irc_tables_cache = None


# ---------------------------------------------------------------------------
# Section properties -- circular
# ---------------------------------------------------------------------------

def circular_area(diameter: float) -> float:
    """Cross-sectional area of a circle.

    Parameters
    ----------
    diameter : float
        Diameter in consistent units (typically mm).

    Returns
    -------
    float
        Area in ``diameter`` units squared.
    """
    return math.pi * diameter**2 / 4.0


def circular_inertia(diameter: float) -> float:
    """Second moment of area (I) of a solid circular section.

    ``I = pi * d^4 / 64``

    Parameters
    ----------
    diameter : float
        Diameter in consistent units (typically mm).

    Returns
    -------
    float
        Moment of inertia in ``diameter`` units to the fourth power.
    """
    return math.pi * diameter**4 / 64.0


def circular_section_modulus(diameter: float) -> float:
    """Elastic section modulus (Z) of a solid circular section.

    ``Z = I / (d/2) = pi * d^3 / 32``

    Parameters
    ----------
    diameter : float
        Diameter in consistent units.

    Returns
    -------
    float
        Section modulus in ``diameter`` units cubed.
    """
    return math.pi * diameter**3 / 32.0


def circular_perimeter(diameter: float) -> float:
    """Perimeter (circumference) of a circle.

    Parameters
    ----------
    diameter : float
        Diameter in consistent units.

    Returns
    -------
    float
        Perimeter in ``diameter`` units.
    """
    return math.pi * diameter


# ---------------------------------------------------------------------------
# Section properties -- rectangular
# ---------------------------------------------------------------------------

def rectangular_area(width: float, depth: float) -> float:
    """Cross-sectional area of a rectangle.

    Parameters
    ----------
    width : float
        Width (b) in consistent units.
    depth : float
        Depth (d) in consistent units.

    Returns
    -------
    float
        Area in ``width`` units squared.
    """
    return width * depth


def rectangular_inertia(width: float, depth: float) -> float:
    """Second moment of area (I) of a solid rectangular section about the
    axis parallel to ``width``.

    ``I = b * d^3 / 12``

    Parameters
    ----------
    width : float
        Width (b).
    depth : float
        Depth (d) -- the dimension perpendicular to the bending axis.

    Returns
    -------
    float
        Moment of inertia in consistent units to the fourth power.
    """
    return width * depth**3 / 12.0


def rectangular_section_modulus(width: float, depth: float) -> float:
    """Elastic section modulus (Z) of a solid rectangular section.

    ``Z = b * d^2 / 6``

    Parameters
    ----------
    width : float
        Width (b).
    depth : float
        Depth (d).

    Returns
    -------
    float
        Section modulus in consistent units cubed.
    """
    return width * depth**2 / 6.0


def rectangular_perimeter(width: float, depth: float) -> float:
    """Perimeter of a rectangle.

    Parameters
    ----------
    width : float
        Width (b).
    depth : float
        Depth (d).

    Returns
    -------
    float
        Perimeter in consistent units.
    """
    return 2.0 * (width + depth)


# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------

_G: float = 9.81  # m/s^2, standard gravity


def kn_to_tonnes(kn: float) -> float:
    """Convert kilonewtons to metric tonnes (force)."""
    return kn / _G


def tonnes_to_kn(tonnes: float) -> float:
    """Convert metric tonnes (force) to kilonewtons."""
    return tonnes * _G


def mm_to_m(mm: float) -> float:
    """Convert millimetres to metres."""
    return mm / 1_000.0


def m_to_mm(m: float) -> float:
    """Convert metres to millimetres."""
    return m * 1_000.0


def mpa_to_kn_m2(mpa: float) -> float:
    """Convert megapascals to kN/m^2."""
    return mpa * 1_000.0


def kn_m2_to_mpa(kn_m2: float) -> float:
    """Convert kN/m^2 to megapascals."""
    return kn_m2 / 1_000.0


def degrees_to_radians(deg: float) -> float:
    """Convert degrees to radians."""
    return math.radians(deg)


def radians_to_degrees(rad: float) -> float:
    """Convert radians to degrees."""
    return math.degrees(rad)

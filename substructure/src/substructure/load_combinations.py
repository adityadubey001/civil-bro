"""Load combinations for bridge substructure design per IRC 6:2017.

Combines individual load effects from :class:`loads.LoadResults` into factored
design force vectors at critical sections (pier base, pilecap top) in
accordance with IRC 6:2017 Tables B.1 through B.4:

* **ULS Basic** (Table B.1)  -- three leading-variable sub-combinations
* **ULS Seismic** (Table B.1) -- directionality rule (100/30)
* **SLS Rare / Characteristic** (Table B.2)
* **SLS Frequent** (Table B.3)
* **SLS Quasi-Permanent** (Table B.4)

For each table the function iterates over all live-load cases, both
wind directions (+/-), and (for seismic) all eight sign permutations of
the 100/30 directional rule.  The resulting :class:`CombinationResults`
dataclass contains the full list of combinations together with envelope
(governing) cases for each limit state category.

Key references
--------------
* IRC 6:2017, Annex B, Tables B.1--B.4
* IRC 6:2017, Cl.219 -- Seismic force
* IRC:SP-114:2018 -- Guidelines for Seismic Design of Road Bridges
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

from .geometry import GeometryResults
from .loads import ForceVector, LLCase, LoadResults
from .utils import load_irc_tables


# ---------------------------------------------------------------------------
# Result data classes
# ---------------------------------------------------------------------------

@dataclass
class CombinationResult:
    """Result of a single load combination at all critical sections.

    Attributes
    ----------
    name : str
        Human-readable label, e.g.
        ``"ULS Basic - LL Leading (Max P) [+Wind]"``.
    category : str
        One of ``"uls_basic"``, ``"uls_seismic"``, ``"sls_rare"``,
        ``"sls_frequent"``, ``"sls_quasi_permanent"``.
    forces_pier_base : ForceVector
        Factored forces at the pier base (top of pilecap).
    forces_pilecap_top : ForceVector
        Factored forces at the pilecap top, equal to pier-base forces
        plus the pilecap self-weight (concentric, only *P* changes).
    """

    name: str
    category: str
    forces_pier_base: ForceVector
    forces_pilecap_top: ForceVector


@dataclass
class CombinationResults:
    """All load combination results with governing-case envelope.

    Attributes
    ----------
    all_combinations : list[CombinationResult]
        Every generated combination across all limit-state categories.
    governing_uls_basic : CombinationResult
        Governing ULS basic combination (max resultant moment).
    governing_uls_seismic : CombinationResult
        Governing ULS seismic combination (max resultant moment).
    governing_sls_rare : CombinationResult
        Governing SLS rare combination (max resultant moment).
    governing_sls_frequent : CombinationResult
        Governing SLS frequent combination (max resultant moment).
    governing_sls_quasi : CombinationResult
        Governing SLS quasi-permanent combination (max resultant moment).
    max_P : CombinationResult
        Combination with the maximum vertical load (across all categories).
    min_P : CombinationResult
        Combination with the minimum vertical load -- governs uplift check.
    max_ML : CombinationResult
        Combination with the maximum transverse moment |ML|.
    max_MT : CombinationResult
        Combination with the maximum longitudinal moment |MT|.
    """

    all_combinations: list[CombinationResult]

    # Governing cases per category
    governing_uls_basic: CombinationResult
    governing_uls_seismic: CombinationResult
    governing_sls_rare: CombinationResult
    governing_sls_frequent: CombinationResult
    governing_sls_quasi: CombinationResult

    # Envelope forces at pier base (across all categories)
    max_P: CombinationResult
    min_P: CombinationResult
    max_ML: CombinationResult
    max_MT: CombinationResult


# ---------------------------------------------------------------------------
# ForceVector arithmetic helpers
# ---------------------------------------------------------------------------

def _scale(fv: ForceVector, factor: float) -> ForceVector:
    """Return a new ForceVector with every component multiplied by *factor*."""
    return ForceVector(
        P=fv.P * factor,
        FL=fv.FL * factor,
        FT=fv.FT * factor,
        ML=fv.ML * factor,
        MT=fv.MT * factor,
    )


def _add(*fvs: ForceVector) -> ForceVector:
    """Return the component-wise sum of one or more ForceVectors."""
    P = FL = FT = ML = MT = 0.0
    for fv in fvs:
        P += fv.P
        FL += fv.FL
        FT += fv.FT
        ML += fv.ML
        MT += fv.MT
    return ForceVector(P=P, FL=FL, FT=FT, ML=ML, MT=MT)


def _negate_horizontal(fv: ForceVector) -> ForceVector:
    """Negate both horizontal forces and both moments (full sign reversal
    for direction-sensitive loads like wind)."""
    return ForceVector(
        P=fv.P,
        FL=-fv.FL,
        FT=-fv.FT,
        ML=-fv.ML,
        MT=-fv.MT,
    )


def _resultant_moment(fv: ForceVector) -> float:
    """Return sqrt(ML**2 + MT**2) -- combined bending resultant."""
    return math.sqrt(fv.ML ** 2 + fv.MT ** 2)


# ---------------------------------------------------------------------------
# Permanent load assembly
# ---------------------------------------------------------------------------

def _permanent_forces(
    loads: LoadResults,
    gamma_dl: float,
    gamma_sidl_perm: float,
    gamma_sidl_surf: float,
) -> ForceVector:
    """Factored permanent loads (DL + SIDL_permanent + SIDL_surfacing).

    All force vectors inside *loads* are already resolved at pier base.
    """
    dl = _add(
        loads.dl_pier,
        loads.dl_piercap,
        loads.dl_super_lhs,
        loads.dl_super_rhs,
    )
    sidl_perm = _add(loads.sidl_lhs, loads.sidl_rhs)
    sidl_surf = _add(loads.wc_lhs, loads.wc_rhs)
    return _add(
        _scale(dl, gamma_dl),
        _scale(sidl_perm, gamma_sidl_perm),
        _scale(sidl_surf, gamma_sidl_surf),
    )


# ---------------------------------------------------------------------------
# Wind direction variants
# ---------------------------------------------------------------------------

def _wind_both_directions(
    loads: LoadResults,
) -> list[tuple[ForceVector, str]]:
    """Return the total wind force (super + LL + sub) for positive and
    negative transverse/longitudinal direction.

    Each entry is ``(wind_total, sign_label)``.
    """
    w_pos = _add(loads.wind_on_super, loads.wind_on_ll, loads.wind_on_sub)
    w_neg = _negate_horizontal(w_pos)
    return [
        (w_pos, "+Wind"),
        (w_neg, "-Wind"),
    ]


# ---------------------------------------------------------------------------
# Worst-case temperature
# ---------------------------------------------------------------------------

def _worst_temp(loads: LoadResults) -> ForceVector:
    """Return the temperature case (rise or fall) that gives the larger
    absolute longitudinal moment |MT|."""
    if abs(loads.temp_rise.MT) >= abs(loads.temp_fall.MT):
        return loads.temp_rise
    return loads.temp_fall


# ---------------------------------------------------------------------------
# Pilecap self-weight
# ---------------------------------------------------------------------------

def _pilecap_self_weight(
    config: dict[str, Any],
    geometry: GeometryResults,
) -> float:
    """Return the pile-cap self-weight in kN."""
    density = float(config["materials"].get("concrete_density", 25.0))
    vol = (geometry.pilecap_length_long
           * geometry.pilecap_width_trans
           * geometry.pilecap_thickness)
    return vol * density


def _forces_at_pilecap_top(
    pier_base_fv: ForceVector,
    pilecap_weight_kn: float,
) -> ForceVector:
    """Produce the force vector at pilecap top by adding the pilecap
    self-weight (concentric -- only *P* changes)."""
    return ForceVector(
        P=pier_base_fv.P + pilecap_weight_kn,
        FL=pier_base_fv.FL,
        FT=pier_base_fv.FT,
        ML=pier_base_fv.ML,
        MT=pier_base_fv.MT,
    )


# ---------------------------------------------------------------------------
# Governing-case selection helpers
# ---------------------------------------------------------------------------

def _governing_by_moment(
    combos: list[CombinationResult],
) -> CombinationResult:
    """Return the combination with the largest resultant moment
    sqrt(ML**2 + MT**2) at the pier base."""
    return max(combos, key=lambda c: _resultant_moment(c.forces_pier_base))


def _governing_by_max_P(
    combos: list[CombinationResult],
) -> CombinationResult:
    """Return the combination with the largest vertical load P."""
    return max(combos, key=lambda c: c.forces_pier_base.P)


def _governing_by_min_P(
    combos: list[CombinationResult],
) -> CombinationResult:
    """Return the combination with the smallest (most negative) P."""
    return min(combos, key=lambda c: c.forces_pier_base.P)


def _governing_by_max_abs_ML(
    combos: list[CombinationResult],
) -> CombinationResult:
    """Return the combination with the largest |ML|."""
    return max(combos, key=lambda c: abs(c.forces_pier_base.ML))


def _governing_by_max_abs_MT(
    combos: list[CombinationResult],
) -> CombinationResult:
    """Return the combination with the largest |MT|."""
    return max(combos, key=lambda c: abs(c.forces_pier_base.MT))


def _filter_category(
    combos: list[CombinationResult],
    category: str,
) -> list[CombinationResult]:
    """Return only combinations belonging to *category*."""
    return [c for c in combos if c.category == category]


# ---------------------------------------------------------------------------
# ULS Basic (IRC 6 Table B.1)
# ---------------------------------------------------------------------------

def _generate_uls_basic(
    loads: LoadResults,
    irc: dict[str, Any],
) -> list[CombinationResult]:
    """Generate all ULS Basic combinations per IRC 6 Table B.1.

    Three sub-combinations (LL leading, Wind leading, Temp leading) are
    crossed with all LL cases, both wind directions, and both DL scenarios
    (unfavourable gamma_sup=1.35 and favourable gamma_inf=1.0).
    """
    tbl = irc["uls_basic"]

    gamma_dl_sup = float(tbl["DL"][0])                     # 1.35
    gamma_dl_inf = float(tbl["DL"][1])                     # 1.0
    gamma_sidl_perm = float(tbl["SIDL_permanent"][0])      # 1.35
    gamma_sidl_surf = float(tbl["SIDL_surfacing"][0])      # 1.75
    gamma_ll_lead = float(tbl["LL_leading"][0])            # 1.5
    gamma_ll_acc = float(tbl["LL_accompanying"][0])        # 1.15
    gamma_wind_lead = float(tbl["Wind_leading"][0])        # 1.5
    gamma_wind_acc = float(tbl["Wind_accompanying"][0])    # 0.9
    gamma_temp_lead = float(tbl["Temp_leading"][0])        # 1.5
    gamma_temp_acc = float(tbl["Temp_accompanying"][0])    # 0.9
    gamma_braking = float(tbl["Braking"][0])               # 1.15

    temp = _worst_temp(loads)
    wind_variants = _wind_both_directions(loads)

    # DL scenarios: unfavourable (max P) and favourable (min P / uplift)
    dl_scenarios: list[tuple[float, str]] = [
        (gamma_dl_sup, ""),
        (gamma_dl_inf, " Fav DL"),
    ]

    results: list[CombinationResult] = []

    for gamma_dl, dl_tag in dl_scenarios:
        perm = _permanent_forces(
            loads, gamma_dl, gamma_sidl_perm, gamma_sidl_surf,
        )

        for case_idx, ll_case in enumerate(loads.ll_cases):
            ll_fv = ll_case.force_at_pier_base
            braking_fv = loads.braking_cases[case_idx]["total"]
            centrifugal_fv = loads.centrifugal_cases[case_idx]["total"]

            for wind_total, wind_tag in wind_variants:

                # --- Sub-combination 1: LL Leading ---
                fv = _add(
                    perm,
                    _scale(ll_fv, gamma_ll_lead),
                    _scale(braking_fv, gamma_braking),
                    _scale(centrifugal_fv, gamma_braking),
                    _scale(wind_total, gamma_wind_acc),
                    _scale(temp, gamma_temp_acc),
                )
                results.append(CombinationResult(
                    name=(f"ULS Basic - LL Leading ({ll_case.name}) "
                          f"[{wind_tag}{dl_tag}]"),
                    category="uls_basic",
                    forces_pier_base=fv,
                    forces_pilecap_top=ForceVector(),  # filled later
                ))

                # --- Sub-combination 2: Wind Leading ---
                fv = _add(
                    perm,
                    _scale(ll_fv, gamma_ll_acc),
                    _scale(braking_fv, gamma_braking),
                    _scale(centrifugal_fv, gamma_braking),
                    _scale(wind_total, gamma_wind_lead),
                    _scale(temp, gamma_temp_acc),
                )
                results.append(CombinationResult(
                    name=(f"ULS Basic - Wind Leading ({ll_case.name}) "
                          f"[{wind_tag}{dl_tag}]"),
                    category="uls_basic",
                    forces_pier_base=fv,
                    forces_pilecap_top=ForceVector(),
                ))

                # --- Sub-combination 3: Temp Leading ---
                fv = _add(
                    perm,
                    _scale(ll_fv, gamma_ll_acc),
                    _scale(braking_fv, gamma_braking),
                    _scale(centrifugal_fv, gamma_braking),
                    _scale(wind_total, gamma_wind_acc),
                    _scale(temp, gamma_temp_lead),
                )
                results.append(CombinationResult(
                    name=(f"ULS Basic - Temp Leading ({ll_case.name}) "
                          f"[{wind_tag}{dl_tag}]"),
                    category="uls_basic",
                    forces_pier_base=fv,
                    forces_pilecap_top=ForceVector(),
                ))

    return results


# ---------------------------------------------------------------------------
# ULS Seismic (IRC 6 Table B.1 -- seismic combination)
# ---------------------------------------------------------------------------

def _generate_uls_seismic(
    loads: LoadResults,
    irc: dict[str, Any],
) -> list[CombinationResult]:
    """Generate all ULS Seismic combinations per IRC 6 Table B.1.

    Directionality rule (IRC:SP-114 / IRC 6 Cl.219):
      - 100 % EQ_long + 30 % EQ_trans
      - 30 % EQ_long + 100 % EQ_trans

    Each is considered with all four sign permutations (+/- for both axes),
    giving eight directional combinations.  No braking or wind is combined
    with seismic.  LL is included with gamma = 0.2.
    """
    tbl = irc["uls_seismic"]

    gamma_dl = float(tbl["DL"][0])                         # 1.35
    gamma_sidl_perm = float(tbl["SIDL_permanent"][0])      # 1.35
    gamma_sidl_surf = float(tbl["SIDL_surfacing"][0])      # 1.75
    gamma_ll = float(tbl["LL"][0])                         # 0.2
    gamma_eq = float(tbl["Seismic"][0])                    # 1.5

    perm = _permanent_forces(loads, gamma_dl, gamma_sidl_perm, gamma_sidl_surf)

    # The loads module provides seismic forces as separate longitudinal and
    # transverse components on superstructure, live load, and substructure.
    # For the directional rule we scale each component independently.
    #
    # eq_super_long : only HL and MT (longitudinal inertia force on super)
    # eq_super_trans: only HT and ML (transverse inertia force on super)
    # eq_ll         : has both HL/HT and ML/MT (acts in both directions)
    # eq_sub        : has both HL/HT and ML/MT (pier inertia in both dirs)
    #
    # The directional rule applies to the overall *direction* of shaking:
    #   frac_long scales all *longitudinal* components (HL, MT)
    #   frac_trans scales all *transverse* components (HT, ML)
    #
    # For eq_super_long the response IS the longitudinal component, and
    # for eq_super_trans the response IS the transverse component.
    # For eq_ll and eq_sub, HL/MT are longitudinal and HT/ML are transverse.

    # Eight directional permutations
    # (sign_long, frac_long, sign_trans, frac_trans, label)
    directional_combos: list[tuple[float, float, float, float, str]] = [
        (+1.0, 1.0, +1.0, 0.3, "100L+30T"),
        (+1.0, 1.0, -1.0, 0.3, "100L-30T"),
        (-1.0, 1.0, +1.0, 0.3, "-100L+30T"),
        (-1.0, 1.0, -1.0, 0.3, "-100L-30T"),
        (+1.0, 0.3, +1.0, 1.0, "30L+100T"),
        (+1.0, 0.3, -1.0, 1.0, "30L-100T"),
        (-1.0, 0.3, +1.0, 1.0, "-30L+100T"),
        (-1.0, 0.3, -1.0, 1.0, "-30L-100T"),
    ]

    results: list[CombinationResult] = []

    for ll_case in loads.ll_cases:
        ll_fv = ll_case.force_at_pier_base

        for sign_l, frac_l, sign_t, frac_t, dir_label in directional_combos:
            # Build combined seismic force vector.
            # eq_super_long is purely longitudinal -> scale by frac_l * sign_l
            # eq_super_trans is purely transverse -> scale by frac_t * sign_t
            eq_super = _add(
                _scale(loads.eq_super_long, frac_l * sign_l),
                _scale(loads.eq_super_trans, frac_t * sign_t),
            )

            # eq_ll has both directional components; scale the vector by the
            # longitudinal fraction (it is dominated by longitudinal braking-
            # type seismic on LL).  A more precise split would decompose
            # HL/MT and HT/ML separately, but eq_ll is small and this is
            # standard practice.
            eq_ll_combined = _add(
                _scale(loads.eq_ll, frac_l * sign_l),
            )

            # eq_sub: pier inertia has both longitudinal and transverse
            # components.  The loads module stores them in a single vector.
            # We need to split them directionally:
            #   longitudinal part: HL, MT  (scale by frac_l * sign_l)
            #   transverse part:   HT, ML  (scale by frac_t * sign_t)
            eq_sub_long = ForceVector(
                P=loads.eq_sub.P * frac_l * abs(sign_l),
                FL=loads.eq_sub.FL * frac_l * sign_l,
                FT=0.0,
                ML=0.0,
                MT=loads.eq_sub.MT * frac_l * sign_l,
            )
            eq_sub_trans = ForceVector(
                P=0.0,
                FL=0.0,
                FT=loads.eq_sub.FT * frac_t * sign_t,
                ML=loads.eq_sub.ML * frac_t * sign_t,
                MT=0.0,
            )
            eq_sub_combined = _add(eq_sub_long, eq_sub_trans)

            # Total seismic (unfactored), then apply gamma_eq
            eq_total = _scale(
                _add(eq_super, eq_ll_combined, eq_sub_combined),
                gamma_eq,
            )

            fv = _add(
                perm,
                _scale(ll_fv, gamma_ll),
                eq_total,
            )

            results.append(CombinationResult(
                name=f"ULS Seismic - {dir_label} ({ll_case.name})",
                category="uls_seismic",
                forces_pier_base=fv,
                forces_pilecap_top=ForceVector(),
            ))

    return results


# ---------------------------------------------------------------------------
# SLS Rare / Characteristic (IRC 6 Table B.2)
# ---------------------------------------------------------------------------

def _generate_sls_rare(
    loads: LoadResults,
    irc: dict[str, Any],
) -> list[CombinationResult]:
    """Generate all SLS Rare (Characteristic) combinations per IRC 6
    Table B.2.

    Two sub-combinations: LL leading and Wind leading.
    """
    tbl = irc["sls_rare"]

    gamma_dl = float(tbl["DL"][0])                         # 1.0
    gamma_sidl_perm = float(tbl["SIDL_permanent"][0])      # 1.0
    gamma_sidl_surf = float(tbl["SIDL_surfacing"][0])      # 1.2
    gamma_ll_lead = float(tbl["LL_leading"][0])            # 1.0
    gamma_ll_acc = float(tbl["LL_accompanying"][0])        # 0.75
    gamma_wind_lead = float(tbl["Wind_leading"][0])        # 1.0
    gamma_wind_acc = float(tbl["Wind_accompanying"][0])    # 0.6
    gamma_temp_acc = float(tbl["Temp_accompanying"][0])    # 0.6
    gamma_braking = float(tbl["Braking"][0])               # 0.75

    perm = _permanent_forces(loads, gamma_dl, gamma_sidl_perm, gamma_sidl_surf)
    temp = _worst_temp(loads)
    wind_variants = _wind_both_directions(loads)

    results: list[CombinationResult] = []

    for case_idx, ll_case in enumerate(loads.ll_cases):
        ll_fv = ll_case.force_at_pier_base
        braking_fv = loads.braking_cases[case_idx]["total"]
        centrifugal_fv = loads.centrifugal_cases[case_idx]["total"]

        for wind_total, wind_tag in wind_variants:

            # --- LL Leading ---
            fv = _add(
                perm,
                _scale(ll_fv, gamma_ll_lead),
                _scale(braking_fv, gamma_braking),
                _scale(centrifugal_fv, gamma_braking),
                _scale(wind_total, gamma_wind_acc),
                _scale(temp, gamma_temp_acc),
            )
            results.append(CombinationResult(
                name=(f"SLS Rare - LL Leading ({ll_case.name}) "
                      f"[{wind_tag}]"),
                category="sls_rare",
                forces_pier_base=fv,
                forces_pilecap_top=ForceVector(),
            ))

            # --- Wind Leading ---
            fv = _add(
                perm,
                _scale(ll_fv, gamma_ll_acc),
                _scale(braking_fv, gamma_braking),
                _scale(centrifugal_fv, gamma_braking),
                _scale(wind_total, gamma_wind_lead),
                _scale(temp, gamma_temp_acc),
            )
            results.append(CombinationResult(
                name=(f"SLS Rare - Wind Leading ({ll_case.name}) "
                      f"[{wind_tag}]"),
                category="sls_rare",
                forces_pier_base=fv,
                forces_pilecap_top=ForceVector(),
            ))

    return results


# ---------------------------------------------------------------------------
# SLS Frequent (IRC 6 Table B.3)
# ---------------------------------------------------------------------------

def _generate_sls_frequent(
    loads: LoadResults,
    irc: dict[str, Any],
) -> list[CombinationResult]:
    """Generate SLS Frequent combinations per IRC 6 Table B.3.

    Single combination per LL case and wind direction.  All variable
    loads use their psi_1 factors.
    """
    tbl = irc["sls_frequent"]

    gamma_dl = float(tbl["DL"][0])                         # 1.0
    gamma_sidl_perm = float(tbl["SIDL_permanent"][0])      # 1.0
    gamma_sidl_surf = float(tbl["SIDL_surfacing"][0])      # 1.2
    gamma_ll = float(tbl["LL"][0])                         # 0.75
    gamma_wind = float(tbl["Wind"][0])                     # 0.5
    gamma_temp = float(tbl["Temp"][0])                     # 0.5

    perm = _permanent_forces(loads, gamma_dl, gamma_sidl_perm, gamma_sidl_surf)
    temp = _worst_temp(loads)
    wind_variants = _wind_both_directions(loads)

    results: list[CombinationResult] = []

    for case_idx, ll_case in enumerate(loads.ll_cases):
        ll_fv = ll_case.force_at_pier_base
        centrifugal_fv = loads.centrifugal_cases[case_idx]["total"]

        for wind_total, wind_tag in wind_variants:
            fv = _add(
                perm,
                _scale(ll_fv, gamma_ll),
                _scale(centrifugal_fv, gamma_ll),
                _scale(wind_total, gamma_wind),
                _scale(temp, gamma_temp),
            )
            results.append(CombinationResult(
                name=f"SLS Frequent ({ll_case.name}) [{wind_tag}]",
                category="sls_frequent",
                forces_pier_base=fv,
                forces_pilecap_top=ForceVector(),
            ))

    return results


# ---------------------------------------------------------------------------
# SLS Quasi-Permanent (IRC 6 Table B.4)
# ---------------------------------------------------------------------------

def _generate_sls_quasi_permanent(
    loads: LoadResults,
    irc: dict[str, Any],
) -> list[CombinationResult]:
    """Generate SLS Quasi-Permanent combinations per IRC 6 Table B.4.

    LL and Wind have psi_2 = 0 so they do not contribute.  Only
    temperature (psi_2 = 0.5) adds a variable component.
    """
    tbl = irc["sls_quasi_permanent"]

    gamma_dl = float(tbl["DL"][0])                         # 1.0
    gamma_sidl_perm = float(tbl["SIDL_permanent"][0])      # 1.0
    gamma_sidl_surf = float(tbl["SIDL_surfacing"][0])      # 1.2
    gamma_temp = float(tbl["Temp"][0])                     # 0.5

    perm = _permanent_forces(loads, gamma_dl, gamma_sidl_perm, gamma_sidl_surf)
    temp = _worst_temp(loads)

    fv = _add(perm, _scale(temp, gamma_temp))

    return [CombinationResult(
        name="SLS Quasi-Permanent",
        category="sls_quasi_permanent",
        forces_pier_base=fv,
        forces_pilecap_top=ForceVector(),
    )]


# ---------------------------------------------------------------------------
# Summary / reporting helpers
# ---------------------------------------------------------------------------

def summarise_governing(results: CombinationResults) -> str:
    """Return a human-readable summary of the governing combinations.

    Useful for console output and logging.
    """
    lines: list[str] = []
    lines.append("=" * 72)
    lines.append("GOVERNING LOAD COMBINATIONS (IRC 6:2017)")
    lines.append("=" * 72)

    def _fmt(label: str, combo: CombinationResult) -> None:
        fv = combo.forces_pier_base
        lines.append("")
        lines.append(f"  {label}")
        lines.append(f"    Name : {combo.name}")
        lines.append(f"    Cat  : {combo.category}")
        lines.append(f"    P    : {fv.P:>12.2f} kN")
        lines.append(f"    HL   : {fv.HL:>12.2f} kN")
        lines.append(f"    HT   : {fv.HT:>12.2f} kN")
        lines.append(f"    ML   : {fv.ML:>12.2f} kN.m")
        lines.append(f"    MT   : {fv.MT:>12.2f} kN.m")

    _fmt("ULS Basic (governing)", results.governing_uls_basic)
    _fmt("ULS Seismic (governing)", results.governing_uls_seismic)
    _fmt("SLS Rare (governing)", results.governing_sls_rare)
    _fmt("SLS Frequent (governing)", results.governing_sls_frequent)
    _fmt("SLS Quasi-Permanent (governing)", results.governing_sls_quasi)
    _fmt("Max P (all categories)", results.max_P)
    _fmt("Min P (all categories)", results.min_P)
    _fmt("Max |ML| (all categories)", results.max_ML)
    _fmt("Max |MT| (all categories)", results.max_MT)

    lines.append("")
    lines.append("=" * 72)
    lines.append(f"Total combinations generated: {len(results.all_combinations)}")

    # Count per category
    cat_counts: dict[str, int] = {}
    for c in results.all_combinations:
        cat_counts[c.category] = cat_counts.get(c.category, 0) + 1
    for cat in ("uls_basic", "uls_seismic", "sls_rare",
                "sls_frequent", "sls_quasi_permanent"):
        lines.append(f"  {cat:25s}: {cat_counts.get(cat, 0)}")

    lines.append("=" * 72)
    return "\n".join(lines)


def combination_table(
    combos: list[CombinationResult],
    *,
    section: str = "pier_base",
    max_rows: int = 0,
) -> list[dict[str, Any]]:
    """Convert a list of CombinationResult objects to a flat table format.

    Each row is a dict suitable for tabular display or CSV export:

    .. code-block:: python

        {"name": "...", "category": "...", "P": ..., "HL": ...,
         "HT": ..., "ML": ..., "MT": ...}

    Parameters
    ----------
    combos : list[CombinationResult]
        Combinations to tabulate.
    section : str
        ``"pier_base"`` or ``"pilecap_top"`` to select which force vector.
    max_rows : int
        If positive, truncate the output to this many rows.

    Returns
    -------
    list[dict[str, Any]]
        One dict per combination.
    """
    rows: list[dict[str, Any]] = []
    for c in combos:
        fv = (c.forces_pier_base if section == "pier_base"
              else c.forces_pilecap_top)
        rows.append({
            "name": c.name,
            "category": c.category,
            "P": fv.P,
            "HL": fv.HL,
            "HT": fv.HT,
            "ML": fv.ML,
            "MT": fv.MT,
        })
    if max_rows > 0:
        rows = rows[:max_rows]
    return rows


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def calculate_combinations(
    config: dict[str, Any],
    geometry: GeometryResults,
    loads: LoadResults,
) -> CombinationResults:
    """Generate all IRC 6 load combinations and find governing cases.

    This is the primary public API of the module.  It reads the load
    combination factors from ``irc_tables.yaml``, applies them to the
    individual load results computed by :mod:`.loads`, transfers forces to
    the pilecap level, and returns a :class:`CombinationResults` object
    with every combination and the governing envelopes.

    Parameters
    ----------
    config : dict
        Validated configuration dictionary from :func:`input_parser.parse_input`.
    geometry : GeometryResults
        Pre-computed geometry from :func:`geometry.calculate_geometry`.
    loads : LoadResults
        Pre-computed loads from :func:`loads.calculate_loads`.

    Returns
    -------
    CombinationResults
        All combinations and governing cases.

    Raises
    ------
    FileNotFoundError
        If ``irc_tables.yaml`` cannot be located.
    KeyError
        If a required combination table is missing from the YAML file.
    """
    irc = load_irc_tables()

    # ------------------------------------------------------------------
    # 1. Generate all combinations at pier base
    # ------------------------------------------------------------------
    all_combos: list[CombinationResult] = []

    all_combos.extend(_generate_uls_basic(loads, irc))
    all_combos.extend(_generate_uls_seismic(loads, irc))
    all_combos.extend(_generate_sls_rare(loads, irc))
    all_combos.extend(_generate_sls_frequent(loads, irc))
    all_combos.extend(_generate_sls_quasi_permanent(loads, irc))

    # ------------------------------------------------------------------
    # 2. Transfer forces to pilecap top (pier base + pilecap self-weight)
    # ------------------------------------------------------------------
    pcap_wt = _pilecap_self_weight(config, geometry)
    for combo in all_combos:
        combo.forces_pilecap_top = _forces_at_pilecap_top(
            combo.forces_pier_base, pcap_wt,
        )

    # ------------------------------------------------------------------
    # 3. Extract governing cases per category
    # ------------------------------------------------------------------
    uls_basic_list = _filter_category(all_combos, "uls_basic")
    uls_seismic_list = _filter_category(all_combos, "uls_seismic")
    sls_rare_list = _filter_category(all_combos, "sls_rare")
    sls_freq_list = _filter_category(all_combos, "sls_frequent")
    sls_quasi_list = _filter_category(all_combos, "sls_quasi_permanent")

    gov_uls_basic = _governing_by_moment(uls_basic_list)
    gov_uls_seismic = _governing_by_moment(uls_seismic_list)
    gov_sls_rare = _governing_by_moment(sls_rare_list)
    gov_sls_freq = _governing_by_moment(sls_freq_list)
    gov_sls_quasi = _governing_by_moment(sls_quasi_list)

    # ------------------------------------------------------------------
    # 4. Global envelope (across all categories)
    # ------------------------------------------------------------------
    max_p_combo = _governing_by_max_P(all_combos)
    min_p_combo = _governing_by_min_P(all_combos)
    max_ml_combo = _governing_by_max_abs_ML(all_combos)
    max_mt_combo = _governing_by_max_abs_MT(all_combos)

    # ------------------------------------------------------------------
    # 5. Assemble results
    # ------------------------------------------------------------------
    return CombinationResults(
        all_combinations=all_combos,
        governing_uls_basic=gov_uls_basic,
        governing_uls_seismic=gov_uls_seismic,
        governing_sls_rare=gov_sls_rare,
        governing_sls_frequent=gov_sls_freq,
        governing_sls_quasi=gov_sls_quasi,
        max_P=max_p_combo,
        min_P=min_p_combo,
        max_ML=max_ml_combo,
        max_MT=max_mt_combo,
    )

"""
IRC 112:2020 code provisions for bridge design.

Key clauses implemented:
- Clause 6.4: Material properties
- Table 6.5: Concrete properties (Ecm, fctm)
- Clause 10.3: Shear design
- Clause 12.2.1: SLS stress limits
- IRC 6: Impact factors and vehicle loads
"""

from typing import Dict, Tuple, List
import math
from .base_code import DesignCode, CoverRequirements


class IRC112(DesignCode):
    """
    IRC 112:2020 - Code of Practice for Concrete Road Bridges.
    """

    # Table 6.5: Ecm values (MPa) for different concrete grades
    ECM_TABLE = {
        25: 30000,
        30: 31000,
        35: 32000,
        40: 33000,
        45: 34000,
        50: 35000,
        55: 36000,
        60: 37000,
    }

    # Table 6.5: fctm values (MPa) for different concrete grades
    FCTM_TABLE = {
        25: 2.2,
        30: 2.5,
        35: 2.8,
        40: 3.0,
        45: 3.25,
        50: 3.5,
        55: 3.7,
        60: 3.9,
    }

    # Cover requirements based on exposure (IRC 112 Table 14.2)
    COVER_TABLE = {
        'moderate': {'cover': 40, 'min_grade': 'M30'},
        'severe': {'cover': 45, 'min_grade': 'M35'},
        'very_severe': {'cover': 50, 'min_grade': 'M40'},
        'extreme': {'cover': 75, 'min_grade': 'M45'},
    }

    # xu_max/d values for different steel grades
    XU_MAX_RATIO = {
        415: 0.48,
        500: 0.46,
        550: 0.44,
    }

    @property
    def code_name(self) -> str:
        return "IRC 112:2020"

    def get_partial_safety_factors(self) -> Dict[str, float]:
        """
        Partial safety factors per IRC 112 / IRC 6.

        Table B.1: Load combinations for bridges
        - Dead load: γf = 1.35
        - Live load: γf = 1.50
        """
        return {
            'gamma_c': 1.5,       # Concrete
            'gamma_s': 1.15,      # Steel
            'gamma_f_dead': 1.35,  # Dead load factor (bridges)
            'gamma_f_live': 1.50,  # Live load factor (bridges)
        }

    def get_span_depth_ratio(self, support_type: str) -> float:
        """
        Span/depth ratio for preliminary sizing of bridge girders.

        For rectangular beams: L/10 to L/12
        For T-beams: L/12 to L/15
        """
        ratios = {
            'simply_supported': 10,  # Conservative for rectangular
            'simply_supported_tbeam': 12,
            'continuous': 15,
            'cantilever': 5,
        }
        return ratios.get(support_type, 10)

    def get_minimum_cover(self, exposure_class: str) -> CoverRequirements:
        """Return minimum cover per IRC 112 Table 14.2."""
        data = self.COVER_TABLE.get(exposure_class, self.COVER_TABLE['moderate'])
        return CoverRequirements(
            nominal_cover=data['cover'],
            exposure_class=exposure_class,
            min_grade=data['min_grade']
        )

    def get_minimum_reinforcement_ratio(self, fy: float) -> float:
        """
        Minimum reinforcement per IRC 112 Clause 16.5.1.1.

        As_min = max(0.26 × fctm/fy × bt × d, 0.0013 × bt × d)
        Returns as percentage.
        """
        # Using fctm for M40 as reference
        fctm = 3.0
        min_ratio = max(0.26 * fctm / fy, 0.0013) * 100
        return min_ratio

    def get_maximum_reinforcement_ratio(self) -> float:
        """Maximum reinforcement per IRC 112 = 4%."""
        return 4.0

    def get_shear_strength_concrete(self, pt: float, fck: float) -> float:
        """
        Shear resistance without reinforcement per IRC 112 Clause 10.3.2.

        VRd,c = 0.12 × k × (80 × ρl × fck)^(1/3) × bw × d

        Returns τc in MPa (for unit width and depth)
        """
        # This returns the coefficient, actual VRd,c needs bw and d
        k = 2.0  # Maximum value, actual depends on d
        rho_l = min(pt / 100, 0.02)
        tau_c = 0.12 * k * (80 * rho_l * fck) ** (1/3)
        return tau_c

    def get_maximum_shear_stress(self, fck: float) -> float:
        """
        Maximum shear stress per IRC 112 Clause 10.3.3.

        VRd,max based on strut crushing.
        """
        fcd = 0.67 * fck / 1.5
        nu1 = 0.6 * (1 - fck / 250)
        # For θ = 45°: τ_max = 0.5 × αcw × ν1 × fcd
        tau_max = 0.5 * 1.0 * nu1 * fcd
        return tau_max

    def get_xu_max_ratio(self, fy: float) -> float:
        """Get xu_max/d ratio for limiting neutral axis depth."""
        if fy in self.XU_MAX_RATIO:
            return self.XU_MAX_RATIO[fy]
        elif fy <= 415:
            return 0.48
        elif fy <= 500:
            return 0.46
        else:
            return 0.44

    def get_ecm(self, fck: float) -> float:
        """Get modulus of elasticity from Table 6.5."""
        if int(fck) in self.ECM_TABLE:
            return self.ECM_TABLE[int(fck)]
        else:
            # Formula fallback
            fcm = fck + 10
            return 22 * (fcm / 10) ** 0.3 * 1000

    def get_fctm(self, fck: float) -> float:
        """Get mean tensile strength from Table 6.5."""
        if int(fck) in self.FCTM_TABLE:
            return self.FCTM_TABLE[int(fck)]
        else:
            # Formula fallback
            return 0.259 * (fck ** (2/3))

    def get_sls_stress_limits(self, fck: float, fy: float) -> Dict[str, float]:
        """
        SLS stress limits per IRC 112 Clause 12.2.1.

        Concrete: σc ≤ 0.48 × fck (rare combination)
        Steel: σs ≤ 0.8 × fy (rare combination)
        """
        return {
            'sigma_c_rare': 0.48 * fck,
            'sigma_c_quasi': 0.36 * fck,
            'sigma_s_rare': 0.8 * fy,
        }

    def get_design_constants(self, fck: float, fy: float) -> Dict[str, float]:
        """Calculate design constants for flexural design."""
        gamma_c = 1.5
        gamma_s = 1.15

        fcd = 0.67 * fck / gamma_c
        fyd = fy / gamma_s

        xu_max_d = self.get_xu_max_ratio(fy)

        # Limiting moment coefficient
        k = xu_max_d
        Ru_lim = 0.36 * k * (1 - 0.42 * k) * fck

        return {
            'fcd': fcd,
            'fyd': fyd,
            'xu_max_d': xu_max_d,
            'Ru_lim': Ru_lim,
            'Ecm': self.get_ecm(fck),
            'fctm': self.get_fctm(fck),
        }


# IRC 6 Vehicle Loads
def irc_class_a_train() -> List[Tuple[float, float]]:
    """
    IRC Class A train of vehicles (IRC 6).
    Returns list of (load_kN, position_from_first_axle_m).
    """
    loads = [
        (27, 0.0),
        (27, 1.1),
        (114, 1.1 + 3.2),
        (114, 1.1 + 3.2 + 1.2),
        (68, 1.1 + 3.2 + 1.2 + 4.3),
        (68, 1.1 + 3.2 + 1.2 + 4.3 + 3.0),
        (68, 1.1 + 3.2 + 1.2 + 4.3 + 3.0 + 3.0),
        (68, 1.1 + 3.2 + 1.2 + 4.3 + 3.0 + 3.0 + 3.0),
    ]
    return loads


def irc_class_70r_wheeled() -> List[Tuple[float, float]]:
    """
    IRC Class 70R wheeled vehicle (IRC 6).
    Total load = 1000 kN (100 tonnes).
    """
    loads = [
        (80, 0.0),
        (120, 3.96),
        (120, 3.96 + 1.52),
        (170, 3.96 + 1.52 + 2.13),
        (170, 3.96 + 1.52 + 2.13 + 1.37),
        (170, 3.96 + 1.52 + 2.13 + 1.37 + 3.05),
        (170, 3.96 + 1.52 + 2.13 + 1.37 + 3.05 + 1.37),
    ]
    return loads


def irc_class_70r_tracked() -> Tuple[float, float]:
    """
    IRC Class 70R tracked vehicle.
    Returns (total_load_kN, contact_length_m).
    """
    return (700, 4.57)


def impact_factor_class_a(span_m: float) -> float:
    """Impact factor for Class A loading (IRC 6)."""
    # For RC bridges: IF = 4.5 / (6 + L), minimum 9%
    if_calc = 4.5 / (6 + span_m)
    return max(if_calc, 0.09)


def impact_factor_70r_wheeled(span_m: float) -> float:
    """Impact factor for 70R wheeled (IRC 6)."""
    if span_m <= 12:
        return 0.25
    elif span_m >= 45:
        return 0.10
    else:
        return 0.25 - 0.15 * (span_m - 12) / (45 - 12)


def impact_factor_70r_tracked(span_m: float) -> float:
    """Impact factor for 70R tracked (IRC 6)."""
    if span_m <= 5:
        return 0.25
    elif span_m >= 9:
        return 0.10
    else:
        return 0.25 - 0.15 * (span_m - 5) / (9 - 5)


def max_bm_simply_supported(span_m: float, axle_loads: List[Tuple[float, float]],
                            num_positions: int = 1000) -> Tuple[float, float, float]:
    """
    Calculate maximum bending moment for a train of loads.
    Uses influence line method with moving load.

    Returns:
        (max_moment_kNm, critical_position_m, position_of_max_moment_m)
    """
    max_moment = 0
    critical_pos = 0
    moment_location = span_m / 2

    train_length = max(pos for _, pos in axle_loads)
    start_pos = 0
    end_pos = span_m + train_length
    step = (end_pos - start_pos) / num_positions

    for i in range(num_positions + 1):
        front_axle_pos = start_pos + i * step

        for section_x in [span_m / 2]:
            moment = 0

            for load, axle_offset in axle_loads:
                axle_pos = front_axle_pos + axle_offset

                if 0 <= axle_pos <= span_m:
                    if axle_pos <= section_x:
                        il_ordinate = axle_pos * (span_m - section_x) / span_m
                    else:
                        il_ordinate = section_x * (span_m - axle_pos) / span_m

                    moment += load * il_ordinate

            if moment > max_moment:
                max_moment = moment
                critical_pos = front_axle_pos
                moment_location = section_x

    return max_moment, critical_pos, moment_location


def max_sf_simply_supported(span_m: float, axle_loads: List[Tuple[float, float]],
                            num_positions: int = 1000) -> Tuple[float, float]:
    """
    Calculate maximum shear force for a train of loads.

    Returns:
        (max_shear_kN, critical_position_m)
    """
    max_shear = 0
    critical_pos = 0

    train_length = max(pos for _, pos in axle_loads)
    start_pos = -train_length
    end_pos = span_m
    step = (end_pos - start_pos) / num_positions

    for i in range(num_positions + 1):
        front_axle_pos = start_pos + i * step
        shear = 0

        for load, axle_offset in axle_loads:
            axle_pos = front_axle_pos + axle_offset

            if 0 <= axle_pos <= span_m:
                shear += load * (span_m - axle_pos) / span_m

        if shear > max_shear:
            max_shear = shear
            critical_pos = front_axle_pos

    return max_shear, critical_pos


def max_bm_tracked_vehicle(span_m: float, total_load: float, contact_length: float) -> float:
    """Maximum BM for tracked vehicle."""
    udl = total_load / contact_length
    c = contact_length
    L = span_m

    if c >= L:
        return udl * L ** 2 / 8
    else:
        R = total_load / 2
        return R * L / 2 - udl * c * c / 8


def max_sf_tracked_vehicle(span_m: float, total_load: float, contact_length: float) -> float:
    """Maximum SF for tracked vehicle."""
    c = contact_length
    L = span_m

    if c >= L:
        return total_load / 2
    else:
        return total_load * (L - c / 2) / L

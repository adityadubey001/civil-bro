"""
IS 456:2000 code provisions for reinforced concrete design.

Key clauses implemented:
- Clause 23.2: Effective span and span/depth ratios
- Clause 26.5.1: Minimum and maximum reinforcement
- Clause 38: Limit state of collapse - Flexure
- Clause 40: Limit state of collapse - Shear
- Clause 43: Limit state of serviceability - Deflection
- Table 16: Cover requirements
- Table 19: Design shear strength of concrete
- Table 20: Maximum shear stress
"""

from typing import Dict
import numpy as np
from .base_code import DesignCode, CoverRequirements


class IS456(DesignCode):
    """
    IS 456:2000 - Indian Standard for Plain and Reinforced Concrete.
    Code of Practice for Plain and Reinforced Concrete (Fourth Revision).
    """

    # Table 16: Nominal cover (mm) for different exposure conditions
    COVER_TABLE = {
        'mild': {'cover': 20, 'min_grade': 'M20'},
        'moderate': {'cover': 30, 'min_grade': 'M25'},
        'severe': {'cover': 45, 'min_grade': 'M30'},
        'very_severe': {'cover': 50, 'min_grade': 'M35'},
        'extreme': {'cover': 75, 'min_grade': 'M40'},
    }

    # Table 19: Design shear strength of concrete τc (N/mm²)
    # Interpolation table for pt% vs fck
    SHEAR_STRENGTH_PT = [0.15, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00]
    SHEAR_STRENGTH_TABLE = {
        # fck: [τc values for each pt]
        15: [0.28, 0.35, 0.46, 0.54, 0.60, 0.64, 0.68, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71],
        20: [0.28, 0.36, 0.48, 0.56, 0.62, 0.67, 0.72, 0.75, 0.79, 0.81, 0.82, 0.82, 0.82],
        25: [0.29, 0.36, 0.49, 0.57, 0.64, 0.70, 0.74, 0.78, 0.82, 0.85, 0.88, 0.90, 0.92],
        30: [0.29, 0.37, 0.50, 0.59, 0.66, 0.71, 0.76, 0.80, 0.84, 0.88, 0.91, 0.94, 0.96],
        35: [0.29, 0.37, 0.50, 0.59, 0.67, 0.73, 0.78, 0.82, 0.86, 0.90, 0.93, 0.96, 0.99],
        40: [0.30, 0.38, 0.51, 0.60, 0.68, 0.74, 0.79, 0.84, 0.88, 0.92, 0.95, 0.98, 1.01],
    }

    # Table 20: Maximum shear stress τc,max (N/mm²)
    MAX_SHEAR_STRESS = {
        15: 2.5,
        20: 2.8,
        25: 3.1,
        30: 3.5,
        35: 3.7,
        40: 4.0,
    }

    # xu_max/d values for different steel grades (Clause 38.1, Note)
    XU_MAX_RATIO = {
        250: 0.53,
        415: 0.48,
        500: 0.46,
        550: 0.44,
    }

    # Basic span/depth ratios (Clause 23.2.1)
    BASIC_LD_RATIOS = {
        'cantilever': 7,
        'simply_supported': 20,
        'continuous': 26,
    }

    @property
    def code_name(self) -> str:
        return "IS 456:2000"

    def get_partial_safety_factors(self) -> Dict[str, float]:
        """
        Partial safety factors per IS 456.

        Clause 36.4.2: Partial safety factors for materials
        - Concrete: γm = 1.5
        - Steel: γm = 1.15

        Table 18: Partial safety factors for loads
        - Dead load: γf = 1.5
        - Live load: γf = 1.5
        """
        return {
            'gamma_c': 1.5,       # Concrete
            'gamma_s': 1.15,     # Steel
            'gamma_f_dead': 1.5,  # Dead load factor
            'gamma_f_live': 1.5,  # Live load factor
        }

    def get_span_depth_ratio(self, support_type: str) -> float:
        """
        Basic span/effective depth ratio per Clause 23.2.1.

        Args:
            support_type: 'simply_supported', 'continuous', 'cantilever'

        Returns:
            Basic L/d ratio
        """
        return self.BASIC_LD_RATIOS.get(support_type, 20)

    def get_minimum_cover(self, exposure_class: str) -> CoverRequirements:
        """
        Return nominal cover per Table 16.

        Note: Cover should be increased by 5mm if concrete grade is
        one grade lower than minimum specified.
        """
        data = self.COVER_TABLE.get(exposure_class, self.COVER_TABLE['moderate'])
        return CoverRequirements(
            nominal_cover=data['cover'],
            exposure_class=exposure_class,
            min_grade=data['min_grade']
        )

    def get_minimum_reinforcement_ratio(self, fy: float) -> float:
        """
        Minimum reinforcement per Clause 26.5.1.1.

        As_min/bd >= 0.85/fy
        But not less than:
        - 0.12% for Fe415 and Fe500
        - 0.15% for mild steel (Fe250)

        Args:
            fy: Steel yield strength in MPa

        Returns:
            Minimum reinforcement ratio as percentage
        """
        min_ratio_formula = 0.85 / fy * 100  # Convert to percentage

        if fy <= 250:
            min_ratio_abs = 0.15
        else:
            min_ratio_abs = 0.12

        return max(min_ratio_formula, min_ratio_abs)

    def get_maximum_reinforcement_ratio(self) -> float:
        """
        Maximum reinforcement per Clause 26.5.1.1(b).

        Maximum area of tension reinforcement shall not exceed 4% of
        gross cross-sectional area.

        Returns:
            Maximum reinforcement ratio as percentage (4.0%)
        """
        return 4.0

    def get_shear_strength_concrete(self, pt: float, fck: float) -> float:
        """
        Design shear strength of concrete per Table 19.

        Interpolates between table values for given pt and fck.

        Args:
            pt: Percentage of tension reinforcement (100*As/bd)
            fck: Characteristic compressive strength of concrete (MPa)

        Returns:
            τc in MPa (N/mm²)
        """
        # Clamp pt to table range
        pt = max(0.15, min(pt, 3.0))

        # Get fck keys
        fck_values = sorted(self.SHEAR_STRENGTH_TABLE.keys())

        # Find appropriate fck (use lower bound if between grades)
        if fck <= fck_values[0]:
            fck_use = fck_values[0]
        elif fck >= fck_values[-1]:
            fck_use = fck_values[-1]
        else:
            fck_use = max(f for f in fck_values if f <= fck)

        # Interpolate for pt
        pt_values = self.SHEAR_STRENGTH_PT
        tau_values = self.SHEAR_STRENGTH_TABLE[fck_use]

        tau_c = np.interp(pt, pt_values, tau_values)

        return float(tau_c)

    def get_maximum_shear_stress(self, fck: float) -> float:
        """
        Maximum shear stress per Table 20.

        Args:
            fck: Characteristic compressive strength (MPa)

        Returns:
            τc,max in MPa
        """
        fck_values = sorted(self.MAX_SHEAR_STRESS.keys())

        if fck <= fck_values[0]:
            return self.MAX_SHEAR_STRESS[fck_values[0]]
        elif fck >= fck_values[-1]:
            return self.MAX_SHEAR_STRESS[fck_values[-1]]
        else:
            # Linear interpolation
            lower = max(f for f in fck_values if f <= fck)
            upper = min(f for f in fck_values if f > fck)
            tau_lower = self.MAX_SHEAR_STRESS[lower]
            tau_upper = self.MAX_SHEAR_STRESS[upper]
            return tau_lower + (tau_upper - tau_lower) * (fck - lower) / (upper - lower)

    def get_xu_max_ratio(self, fy: float) -> float:
        """
        Get xu_max/d ratio for limiting neutral axis depth.

        Clause 38.1, Note: For ductile behavior (under-reinforced section),
        xu <= xu_max where xu_max/d depends on steel grade.

        Args:
            fy: Steel yield strength in MPa

        Returns:
            xu_max/d ratio
        """
        if fy in self.XU_MAX_RATIO:
            return self.XU_MAX_RATIO[fy]
        else:
            # Interpolate or use closest
            fy_values = sorted(self.XU_MAX_RATIO.keys())
            if fy <= fy_values[0]:
                return self.XU_MAX_RATIO[fy_values[0]]
            elif fy >= fy_values[-1]:
                return self.XU_MAX_RATIO[fy_values[-1]]
            else:
                # Linear interpolation
                lower = max(f for f in fy_values if f <= fy)
                upper = min(f for f in fy_values if f > fy)
                ratio_lower = self.XU_MAX_RATIO[lower]
                ratio_upper = self.XU_MAX_RATIO[upper]
                return ratio_lower + (ratio_upper - ratio_lower) * (fy - lower) / (upper - lower)

    def get_modification_factor_tension(self, pt: float, fs: float) -> float:
        """
        Modification factor for tension reinforcement per Clause 23.2.1(e).

        Based on area of tension reinforcement and steel stress.

        Args:
            pt: Percentage of tension reinforcement
            fs: Steel stress at service load (typically 0.58*fy)

        Returns:
            Modification factor (kt)
        """
        # Simplified formula based on IS 456 Fig. 4
        # For fs = 0.58*fy (typical service stress)
        # MF varies from ~2.0 at low pt to ~1.0 at high pt

        if pt <= 0.15:
            return 2.0
        elif pt >= 3.0:
            return 1.0
        else:
            # Approximate formula for Fe500 (fs ≈ 290 MPa)
            # Based on curve fitting of Fig. 4
            mf = 2.0 - 0.4 * pt
            return max(1.0, min(2.0, mf))

    def get_modification_factor_compression(self, pc: float) -> float:
        """
        Modification factor for compression reinforcement per Clause 23.2.1(f).

        Based on area of compression reinforcement.

        Args:
            pc: Percentage of compression reinforcement

        Returns:
            Modification factor (kc)
        """
        if pc <= 0:
            return 1.0
        elif pc >= 3.0:
            return 1.5
        else:
            # Linear approximation
            return 1.0 + 0.167 * pc

    def get_design_constants(self, fck: float, fy: float) -> Dict[str, float]:
        """
        Calculate design constants for flexural design.

        Args:
            fck: Characteristic compressive strength (MPa)
            fy: Characteristic yield strength (MPa)

        Returns:
            Dictionary with design constants
        """
        gamma_c = 1.5
        gamma_s = 1.15

        # Design strengths
        fcd = 0.67 * fck / gamma_c  # Design compressive strength
        fsd = fy / gamma_s  # Design steel strength

        # For rectangular stress block (Clause 38.1)
        xu_max_d = self.get_xu_max_ratio(fy)

        # Limiting moment of resistance coefficient
        # Mu_lim = 0.36 * fck * b * xu_max * (d - 0.42*xu_max)
        # Mu_lim / (fck * b * d²) = 0.36 * (xu_max/d) * (1 - 0.42*(xu_max/d))
        k = xu_max_d
        Ru_lim = 0.36 * k * (1 - 0.42 * k) * fck  # Mu_lim / (b*d²)

        return {
            'fcd': fcd,
            'fsd': fsd,
            'xu_max_d': xu_max_d,
            'Ru_lim': Ru_lim,  # MPa (N/mm²)
        }

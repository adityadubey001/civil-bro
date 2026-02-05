"""
Engineering constants for RC beam design.
"""

# Standard bar sizes in mm
STANDARD_BAR_SIZES = [8, 10, 12, 16, 20, 25, 32]

# Stirrup bar sizes in mm
STIRRUP_BAR_SIZES = [6, 8, 10, 12]

# Concrete grades with characteristic strength (fck in MPa)
CONCRETE_GRADES = {
    "M15": 15,
    "M20": 20,
    "M25": 25,
    "M30": 30,
    "M35": 35,
    "M40": 40,
    "M45": 45,
    "M50": 50,
}

# Steel grades with yield strength (fy in MPa)
STEEL_GRADES = {
    "Fe250": 250,
    "Fe415": 415,
    "Fe500": 500,
    "Fe550": 550,
}

# Unit weight of concrete (kN/m³)
CONCRETE_UNIT_WEIGHT = 25.0

# Modulus of elasticity of steel (MPa)
ES = 200000

# Standard bar areas (mm²)
BAR_AREAS = {
    6: 28.27,
    8: 50.27,
    10: 78.54,
    12: 113.10,
    16: 201.06,
    20: 314.16,
    25: 490.87,
    32: 804.25,
}

# Exposure conditions
EXPOSURE_CONDITIONS = [
    "mild",
    "moderate",
    "severe",
    "very_severe",
    "extreme",
]

# Data models for IRC 112 bridge girder design
from .inputs import (
    BridgeGirderInput, DeckInput,
    ConcreteGrade, SteelGrade, ExposureCondition,
    ReinforcementPreference, LiveLoadType, GirderType
)
from .outputs import (
    BridgeGirderOutput, FlexuralDesignOutput, ShearDesignOutput,
    SLSStressCheckOutput, CrackWidthCheckOutput, DeflectionCheckOutput,
    LoadAnalysisOutput, LiveLoadResult, DesignStatus
)

# Backward compatibility aliases
BeamDesignInput = BridgeGirderInput
BeamDesignOutput = BridgeGirderOutput

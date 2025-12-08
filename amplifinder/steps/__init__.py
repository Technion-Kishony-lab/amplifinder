"""Pipeline steps."""

from amplifinder.steps.base import Step
from amplifinder.steps.initialize import InitializingStep
from amplifinder.steps.get_reference import GetRefGenomeStep
from amplifinder.steps.locate_tns import LocateTNsUsingGenbankStep, LocateTNsUsingISfinderStep
from amplifinder.steps.run_breseq import BreseqStep
from amplifinder.steps.create_reference_junctions import CreateRefTnJcStep, CreateRefTnEndSeqsStep
from amplifinder.steps.create_tnjc import CreateTnJcStep
from amplifinder.steps.create_tnjc2 import CreateTnJc2Step

# Backward compatibility alias
ISfinderStep = LocateTNsUsingISfinderStep

__all__ = [
    "Step",
    "InitializingStep",
    "GetRefGenomeStep",
    "LocateTNsUsingGenbankStep",
    "LocateTNsUsingISfinderStep",
    "ISfinderStep",  # deprecated alias
    "BreseqStep",
    "CreateRefTnJcStep",
    "CreateRefTnEndSeqsStep",
    "CreateTnJcStep",
    "CreateTnJc2Step",
]

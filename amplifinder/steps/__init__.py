"""Pipeline steps."""

from amplifinder.steps.base import Step
from amplifinder.steps.initialize import InitializingStep
from amplifinder.steps.get_reference import GetReferenceStep
from amplifinder.steps.locate_tns_genbank import LocateTNsUsingGenbankStep
from amplifinder.steps.locate_tns_isfinder import LocateTNsUsingISfinderStep
from amplifinder.steps.run_breseq import BreseqStep
from amplifinder.steps.create_reference_junctions import CreateRefTnJcsStep, CreateRefTnEndSeqsStep
from amplifinder.steps.create_tnjc import CreateTNJCStep
from amplifinder.steps.create_tnjc2 import CreateTNJC2Step

# Backward compatibility alias
ISfinderStep = LocateTNsUsingISfinderStep

__all__ = [
    "Step",
    "InitializingStep",
    "GetReferenceStep",
    "LocateTNsUsingGenbankStep",
    "LocateTNsUsingISfinderStep",
    "ISfinderStep",  # deprecated alias
    "BreseqStep",
    "CreateRefTnJcsStep",
    "CreateRefTnEndSeqsStep",
    "CreateTNJCStep",
    "CreateTNJC2Step",
]

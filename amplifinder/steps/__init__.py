"""Pipeline steps."""

from amplifinder.steps.base import Step
from amplifinder.steps.initialize import InitializingStep
from amplifinder.steps.get_reference import GetReferenceStep
from amplifinder.steps.locate_tns_genbank import LocateTNsUsingGenbank
from amplifinder.steps.locate_tns_isfinder import LocateTNsUsingISfinder
from amplifinder.steps.run_breseq import BreseqStep

# Backward compatibility alias
ISfinderStep = LocateTNsUsingISfinder

__all__ = [
    "Step",
    "InitializingStep",
    "GetReferenceStep",
    "LocateTNsUsingGenbank",
    "LocateTNsUsingISfinder",
    "ISfinderStep",  # deprecated alias
    "BreseqStep",
]

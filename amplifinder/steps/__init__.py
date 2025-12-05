"""Pipeline steps."""

from amplifinder.steps.base import Step
from amplifinder.steps.initialize import InitializingStep
from amplifinder.steps.get_reference import GetReferenceStep
from amplifinder.steps.locate_iss_genbank import LocateISsUsingGenbank
from amplifinder.steps.locate_iss_isfinder import LocateISsUsingISfinder
from amplifinder.steps.run_breseq import BreseqStep

# Backward compatibility alias
ISfinderStep = LocateISsUsingISfinder

__all__ = [
    "Step",
    "InitializingStep",
    "GetReferenceStep",
    "LocateISsUsingGenbank",
    "LocateISsUsingISfinder",
    "ISfinderStep",  # deprecated alias
    "BreseqStep",
]

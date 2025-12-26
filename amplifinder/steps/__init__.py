"""Pipeline steps."""

from amplifinder.steps.base import Step
from amplifinder.steps.initialize import InitializingStep
from amplifinder.steps.get_reference import GetRefGenomeStep
from amplifinder.steps.locate_tns import LocateTNsUsingGenbankStep, LocateTNsUsingISfinderStep
from amplifinder.steps.run_breseq import BreseqStep
from amplifinder.steps.create_reference_junctions import CreateRefTnJcStep, CreateRefTnEndSeqsStep
from amplifinder.steps.create_tnjc import CreateTnJcStep
from amplifinder.steps.create_tnjc2 import PairTnJc2Step
from amplifinder.steps.classify_structure import ClassifyTnJc2StructureStep
from amplifinder.steps.filter_candidates import FilterTnJc2CandidatesStep
from amplifinder.steps.calc_amplicon_coverage import CalcTnJc2AmpliconCoverageStep
from amplifinder.steps.synthetic_junctions import CreateSyntheticJunctionsStep
from amplifinder.steps.align_reads import AlignReadsToJunctionsStep
from amplifinder.steps.analyze_alignments import AnalyzeTnJc2AlignmentsStep
from amplifinder.steps.classify_candidates import ClassifyTnJc2CandidatesStep
from amplifinder.steps.export import ExportTnJc2Step

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
    "PairTnJc2Step",
    "CalcTnJc2AmpliconCoverageStep",
    "ClassifyTnJc2StructureStep",
    "FilterTnJc2CandidatesStep",
    "CreateSyntheticJunctionsStep",
    "AlignReadsToJunctionsStep",
    "AnalyzeTnJc2AlignmentsStep",
    "ClassifyTnJc2CandidatesStep",
    "ExportTnJc2Step",
]

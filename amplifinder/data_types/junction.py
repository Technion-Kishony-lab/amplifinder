"""Junction (JC) data model from breseq output."""

from dataclasses import dataclass
from typing import Optional

from amplifinder.data_types.tabularable import Tabularable


@dataclass
class Junction(Tabularable):
    """A junction record from breseq output.
    
    Represents a novel junction between two genomic positions,
    potentially indicating a structural variant.
    """
    
    # Fields to ignore when converting from dict (breseq internal fields)
    _ignored_fields = {"type", "dot", "zero"}
    
    # Required fields
    num: int = 0
    scaf1: str = ""
    pos1: int = 0
    dir1: int = 0  # -1 or 1
    scaf2: str = ""
    pos2: int = 0
    dir2: int = 0  # -1 or 1
    
    # Common optional fields
    flanking_left: int = 0
    flanking_right: int = 0
    reject: Optional[str] = None
    prediction: Optional[str] = None
    new_junction_coverage: Optional[float] = None
    new_junction_read_count: Optional[int] = None
    frequency: Optional[float] = None
    
    # Coverage fields
    coverage_minus: Optional[float] = None
    coverage_plus: Optional[float] = None
    side_1_coverage: Optional[float] = None
    side_2_coverage: Optional[float] = None
    
    @property
    def coverage(self) -> float:
        """Total coverage (minus + plus)."""
        return (self.coverage_minus or 0) + (self.coverage_plus or 0)
    
    @property
    def is_rejected(self) -> bool:
        """Check if junction was rejected by breseq."""
        return self.reject is not None and self.reject != ""

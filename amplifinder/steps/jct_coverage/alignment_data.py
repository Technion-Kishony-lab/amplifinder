from abc import ABC, abstractmethod
from dataclasses import dataclass

from amplifinder.data_types.enums import ReadType
from amplifinder.steps.jct_coverage.cigar import Cigar


@dataclass(frozen=True)
class AlignmentData(ABC):
    """Base class for alignment data."""
    read_type: ReadType
    
    @abstractmethod
    def get_bam_indices(self) -> tuple[int, ...]:
        """Return BAM indices as tuple."""
        pass

    @property
    @abstractmethod
    def left(self) -> int:
        """Return left BAM index."""
        pass
    
    @property
    @abstractmethod
    def right(self) -> int:
        pass

    def get_plotting_segments(self) -> list[tuple[int, int]]:
        pass


@dataclass(frozen=True)
class SingleAlignment(AlignmentData):
    start: int
    end: int
    bam_index: int
    cigar: Cigar
    
    def get_bam_indices(self) -> tuple[int]:
        return (self.bam_index,)

    @property
    def left(self) -> int:
        return self.start

    @property
    def right(self) -> int:
        return self.end

    def get_plotting_segments(self) -> list[tuple[int, int]]:
        return [(self.start, self.end)]


@dataclass(frozen=True)
class PairedAlignment(AlignmentData):
    start1: int
    end1: int
    bam_index1: int
    cigar1: Cigar
    start2: int
    end2: int
    bam_index2: int
    cigar2: Cigar
    
    def get_bam_indices(self) -> tuple[int, int]:
        return (self.bam_index1, self.bam_index2)

    @property
    def left(self) -> int:
        return self.start1

    @property
    def right(self) -> int:
        return self.end2

    def get_plotting_segments(self) -> list[tuple[int, int]]:
        return [(self.start1, self.end1), (self.start2, self.end2)]

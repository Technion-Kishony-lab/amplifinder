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

    @abstractmethod
    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar], ...]:
        """Return tuple of (start, cigar) pairs."""
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

    @property
    def middle(self) -> int:
        return (self.left + self.right) / 2


@dataclass(frozen=True)
class BaseSingleAlignment(AlignmentData):
    """Single alignment from one BAM hit."""
    start: int
    end: int

    @property
    def left(self) -> int:
        return self.start

    @property
    def right(self) -> int:
        return self.end


@dataclass(frozen=True)
class SingleAlignment(BaseSingleAlignment):
    """Single alignment from one BAM hit."""
    cigar: Cigar
    bam_index: int

    def get_bam_indices(self) -> tuple[int]:
        return (self.bam_index,)

    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar]]:
        return ((self.start, self.cigar),)


@dataclass(frozen=True)
class CombinedSingleAlignment(BaseSingleAlignment):
    """Single alignment combined from multiple hits with same read_id and orientation."""
    bam_indices: tuple[int, ...]
    cigars: tuple[Cigar, ...]
    starts: tuple[int, ...]
    ends: tuple[int, ...]

    def get_bam_indices(self) -> tuple[int, ...]:
        return self.bam_indices

    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar], ...]:
        return tuple(zip(self.starts, self.cigars))


@dataclass(frozen=True)
class PairedAlignment(AlignmentData):
    alignment1: BaseSingleAlignment
    alignment2: BaseSingleAlignment

    def get_bam_indices(self) -> tuple[int, ...]:
        return (*self.alignment1.get_bam_indices(), *self.alignment2.get_bam_indices())

    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar], ...]:
        return (*self.alignment1.get_starts_and_cigars(), *self.alignment2.get_starts_and_cigars())

    @property
    def left(self) -> int:
        return self.alignment1.left

    @property
    def right(self) -> int:
        return self.alignment2.right

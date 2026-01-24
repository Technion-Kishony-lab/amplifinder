from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field

from amplifinder.steps.jct_coverage.cigar import Cigar
from amplifinder.data_types import ReadType


@dataclass
class AlignmentData(ABC):
    """Base class for alignment data."""

    read_type: ReadType | None = field(default=None, kw_only=True)

    @abstractmethod
    def get_bam_indices(self) -> tuple[int, ...]:
        """Return BAM indices as tuple."""
        pass

    @abstractmethod
    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar], ...]:
        """Return tuple of (start, cigar) pairs."""
        pass

    @abstractmethod
    def get_all_single_alignments(self) -> tuple["SingleAlignment", ...]:
        """Return all SingleAlignment objects, recursively from nested structures."""
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

    @property
    def length(self) -> int:
        return self.right - self.left + 1


@dataclass
class BaseSingleAlignment(AlignmentData):
    """Single alignment from one BAM hit."""
    start: int
    end: int
    is_reverse: bool

    @property
    def left(self) -> int:
        return self.start

    @property
    def right(self) -> int:
        return self.end

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}(start={self.start}, end={self.end}, is_reverse={self.is_reverse})"


@dataclass
class SingleAlignment(BaseSingleAlignment):
    """Single alignment from one BAM hit."""
    cigar: Cigar
    bam_index: int
    alignment_score: int
    read_id: str = ""

    def get_bam_indices(self) -> tuple[int]:
        return (self.bam_index,)

    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar]]:
        return ((self.start, self.cigar),)

    def get_all_single_alignments(self) -> tuple[SingleAlignment, ...]:
        return (self,)


@dataclass
class CombinedSingleAlignment(BaseSingleAlignment):
    """Single alignment combined from multiple hits with same read_id and orientation."""
    alignments: tuple[SingleAlignment, ...]

    @classmethod
    def from_alignments(cls, alignments: list[SingleAlignment]) -> "CombinedSingleAlignment":
        # Sort by reference start position
        sorted_alignments = tuple(sorted(alignments, key=lambda a: a.start))

        # Overall span
        start = min(a.start for a in sorted_alignments)
        end = max(a.end for a in sorted_alignments)
        assert all(a.is_reverse == sorted_alignments[0].is_reverse for a in sorted_alignments)
        is_reverse = sorted_alignments[0].is_reverse

        return cls(start=start, end=end, is_reverse=is_reverse, alignments=sorted_alignments)

    def get_bam_indices(self) -> tuple[int, ...]:
        return tuple(a.bam_index for a in self.alignments)

    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar], ...]:
        return tuple((a.start, a.cigar) for a in self.alignments)

    def get_all_single_alignments(self) -> tuple[SingleAlignment, ...]:
        return self.alignments


@dataclass
class PairedAlignment(AlignmentData):
    forward_alignment: BaseSingleAlignment
    reverse_alignment: BaseSingleAlignment
    is_swapped: bool = False

    def get_bam_indices(self) -> tuple[int, ...]:
        return (*self.forward_alignment.get_bam_indices(), *self.reverse_alignment.get_bam_indices())

    def get_starts_and_cigars(self) -> tuple[tuple[int, Cigar], ...]:
        return (
            *self.forward_alignment.get_starts_and_cigars(),
            *self.reverse_alignment.get_starts_and_cigars()
        )

    def get_all_single_alignments(self) -> tuple[SingleAlignment, ...]:
        return (
            *self.forward_alignment.get_all_single_alignments(),
            *self.reverse_alignment.get_all_single_alignments()
        )

    @property
    def left(self) -> int:
        return self.forward_alignment.left

    @property
    def right(self) -> int:
        return self.reverse_alignment.right

    @property
    def overlapping_length(self) -> int:
        # >0 if overlapping, 0 if not overlapping, <0 if there is distance between the two alignments
        return self.forward_alignment.right - self.reverse_alignment.left

    def __repr__(self) -> str:
        class_name = self.__class__.__name__
        return f"{class_name}(\n\tfwd={self.forward_alignment},\n\trev={self.reverse_alignment},\n\tis_swapped={self.is_swapped})"

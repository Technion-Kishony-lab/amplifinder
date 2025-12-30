"""
Scaffold coordinate system:

Scaffold is a circular or linear sequence of DNA, with a length and a circularity flag.

Its coordinates are 1-based inclusive.

`start` is the position where the forward strand starts, and 
`end` is the position where the forward strand ends.

Therefore start/end are NOT left/right.

- start <= end: normal segment, span_origin is False

<-------------------length------------------>
~~~~~~~~~~~~|======================|~~~~~~~~~
            start                   end
            |                      |
            ----------------------->           inclusive
                 segment_length    

- start > end and circular scaffold => span_origin is True

<-------------------length------------------>
============|~~~~~~~~~~~~~~~~~~~~~~|=========
            end                   start
            |                      |
            |                      --------->  inclusive
------------>                      |
                                   --------->------------>
                                       segment_length    


- if start > end and scaffold is linear, the segment is not valid => span_origin is None
"""



import numpy as np


from dataclasses import dataclass, field
from typing import TypeVar

from Bio.Seq import reverse_complement

from amplifinder.data_types.enums import Orientation

T = TypeVar('T', str, np.ndarray)


@dataclass
class Scaffold:
    """Coordinate helper with circularity awareness."""

    is_circular: bool
    length: int

    def normalize_range(self, start: int, end: int) -> tuple[int, int, bool | None]:
        """
        Normalize coordinates for range extraction (1-based inclusive).
        Returns:
            start_norm: Normalized start coordinate
            end_norm: Normalized end coordinate
            span_origin: True if segment spans circular origin (start > end and is_circular),
                         False if segment is not circular (start < end)
                         None if segment is not circular but starts after end (start > end and not is_circular)
        """
        if not self.is_circular:
            if start > end:
                return start, end, None
            return start, end, False
        start_norm = (start - 1) % self.length + 1
        end_norm = (end - 1) % self.length + 1
        spans_origin = start_norm > end_norm
        return start_norm, end_norm, spans_origin

    def get_segment_length(self, start: int, end: int) -> int:
        """Length of segment, respecting circular origin rules."""
        start, end, span_origin = self.normalize_range(start, end)
        if span_origin is None:
            return -1
        diff = end - start
        if span_origin:
            return self.length - diff + 1
        return diff + 1

    def __len__(self) -> int:
        return self.length

    @staticmethod
    def _concatenate_sequence(seq1: T, seq2: T) -> T:
        """Concatenate two sequences."""
        if isinstance(seq1, np.ndarray):
            return np.concatenate([seq1, seq2])
        return seq1 + seq2

    def slice(self, start: int, end: int, seq: T, direction: Orientation = Orientation.FORWARD) -> T:
        """Slice arbitrary scaffold-aligned sequence (string or ndarray).

        Args:
            start: Start position (1-based inclusive)
            end: End position (1-based inclusive)
            seq: Sequence (string or ndarray)
            direction: Orientation (reverse only valid for string)
        """
        assert len(seq) == self.length
        start, end, spans_origin = self.normalize_range(start, end)
        if spans_origin is None:
            segment = seq[0:0]  # invalid segment
        if spans_origin:
            segment = self._concatenate_sequence(seq[start - 1:], seq[:end])
        else:
            segment = seq[start - 1:end]
        if direction == Orientation.FORWARD:
            return segment
        if not isinstance(segment, str):
            raise TypeError("reverse slicing requires string sequence")
        return reverse_complement(segment)


@dataclass
class SegmentMixin:
    """Mixin adding start/end helpers."""

    start: int
    end: int

    def slice(self, start: int | None = None, end: int | None = None, seq: T = None, direction: Orientation = Orientation.FORWARD) -> T:
        if seq is None:
            raise ValueError("seq is required")
        start = self.start if start is None else start
        end = self.end if end is None else end
        return super().slice(start, end, seq, direction)

    def normalize_range(self, start: int | None = None, end: int | None = None) -> tuple[int, int, bool | None]:
        start = self.start if start is None else start
        end = self.end if end is None else end
        return super().normalize_range(start, end)

    def get_segment_length(self, start: int | None = None, end: int | None = None) -> int:
        start = self.start if start is None else start
        end = self.end if end is None else end
        return super().get_segment_length(start, end)

    @property
    def segment_length(self) -> int:
        return self.get_segment_length()

    @property
    def span_origin(self) -> bool | None:
        return self.normalize_range()[2]

    @property
    def left(self) -> int:
        return self.normalize_range()[0]

    @property
    def right(self) -> int:
        return self.normalize_range()[1]


@dataclass
class SegmentScaffold(SegmentMixin, Scaffold):
    """Scaffold with persisted start/end."""
    pass


@dataclass
class SeqScaffold(Scaffold):
    """Sequence-backed scaffold."""

    seq: T
    length: int = field(init=False)

    def __post_init__(self):
        self.length = len(self.seq)

    def slice(self, start: int, end: int, seq: T | None = None, direction: Orientation = Orientation.FORWARD) -> T:
        """Get range (1-based inclusive) from sequence."""
        seq = self.seq if seq is None else seq
        return super().slice(start, end, seq, direction)


@dataclass
class SeqSegmentScaffold(SegmentMixin, SeqScaffold):
    """Sequence scaffold with persisted start/end."""
    def slice(self, start: int | None = None, end: int | None = None, seq: T | None = None, direction: Orientation = Orientation.FORWARD) -> T:
        seq = self.seq if seq is None else seq
        return super().slice(start, end, seq, direction)

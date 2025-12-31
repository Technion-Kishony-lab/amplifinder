"""
Scaffold coordinate system:

Scaffold is a circular or linear sequence of DNA, with a length and a circularity flag.

Its coordinates are 1-based inclusive.

`start` is the position where the sequence starts, and 
`end` is the position where the sequence ends.
`orientation` is the orientation of the sequence, specifying the FORWARD or REVERSE strand.

Therefore start/end are NOT left/right.

For FORWARD orientation:
- start <= end: normal segment, span_origin is False

<-------------------length------------------>
~~~~~~~~~~~~|======================|~~~~~~~~~
            start                   end
            |                      |
            ----------------------->           inclusive
                 segment_length    

- start > end and circular scaffold => span_origin is True

============|~~~~~~~~~~~~~~~~~~~~~~|=========
            end                   start
            |                      |
            |                      --------->  inclusive
------------>                      |
                                   --------->------------>
                                       segment_length    


- if start > end and scaffold is linear, the segment is not valid => span_origin is None

For REVERSE orientation:
- start > end: normal reverse-complemented segment, span_origin is False

~~~~~~~~~~~~|======================|~~~~~~~~~
            end                   start
            |                      |
            <-----------------------           inclusive
                 segment_length    

- start < end and scaffold is circular => span_origin is True

============|~~~~~~~~~~~~~~~~~~~~~~|=========
          start                   end
            |                      |
<------------                      |
                                   <---------
"""

from __future__ import annotations

import numpy as np


from dataclasses import dataclass, field
from typing import NamedTuple, TypeVar

from Bio.Seq import reverse_complement

from amplifinder.data_types.enums import Orientation

T = TypeVar('T', str, np.ndarray)


class JcArm(NamedTuple):
    """Junction arm coordinates and orientation."""
    scaf: str
    start: int
    dir: Orientation
    flank: int
    
    @property
    def end(self) -> int:
        """Compute end position based on start position, flank length, and direction."""
        return self.start + self.flank * self.dir
    
    def slice_scaffold(self, scaf_obj: Scaffold, seq: T) -> T:
        """Extract sequence from a Scaffold.
        
        Args:
            scaf_obj: Scaffold with .slice() method
            
        Returns:
            Sequence string with proper orientation
        """
        return scaf_obj.slice(self.start, self.end, self.dir, seq)


@dataclass
class Scaffold:
    """Coordinate helper with circularity awareness."""

    is_circular: bool
    length: int

    def normalize_range(self, start: int, end: int, orientation: Orientation = Orientation.FORWARD) -> tuple[int, int, bool | None]:
        """
        Normalize coordinates for range extraction (1-based inclusive).
        Args:
            start: Start coordinate
            end: End coordinate
            orientation: Orientation.FORWARD or Orientation.REVERSE
        Returns:
            start_norm: Normalized start coordinate
            end_norm: Normalized end coordinate
            span_origin: True if segment spans circular origin (start > end and is_circular),
                         False if segment is not circular (start < end)
                         None if segment is not circular but starts after end (start > end and not is_circular)
        """
        def is_span_origin(start: int, end: int, orientation: Orientation) -> bool:
            if orientation == Orientation.FORWARD:
                return start > end
            if orientation == Orientation.REVERSE:
                return start < end
            assert False, "Invalid orientation"

        if not self.is_circular:
            return start, end, None if is_span_origin(start, end, orientation) else False
        start_norm = (start - 1) % self.length + 1
        end_norm = (end - 1) % self.length + 1
        return start_norm, end_norm, is_span_origin(start_norm, end_norm, orientation)

    def get_segment_length(self, start: int, end: int, orientation: Orientation = Orientation.FORWARD) -> int:
        """Length of segment, respecting circular origin rules."""
        start, end, span_origin = self.normalize_range(start, end, orientation)
        if span_origin is None:
            return -1
        diff = end - start if orientation == Orientation.FORWARD else start - end
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

    def slice(self, start: int, end: int, orientation: Orientation, seq: T) -> T:
        """Slice arbitrary scaffold-aligned sequence (string or ndarray).

        Args:
            start: Start position (1-based inclusive)
            end: End position (1-based inclusive)
            orientation: Orientation (FORWARD or REVERSE, reverse only valid for string)
            seq: Sequence (string or ndarray)
        """
        assert len(seq) == self.length
        start, end, spans_origin = self.normalize_range(start, end, orientation)
        if spans_origin is None:
            segment = seq[0:0]  # invalid segment
        elif spans_origin:
            segment = self._concatenate_sequence(seq[start - 1:], seq[:end])
        else:
            segment = seq[start - 1:end]
        if orientation == Orientation.FORWARD:
            return segment
        if not isinstance(segment, str):
            raise TypeError("reverse slicing requires string sequence")
        return reverse_complement(segment)


@dataclass
class SegmentMixin:
    """Mixin adding start/end helpers."""

    start: int
    end: int
    orientation: Orientation = Orientation.FORWARD

    def _resolve_params(self, start: int | None, end: int | None, 
                        orientation: Orientation | None) -> tuple[int, int, Orientation]:
        """Resolve parameters, using instance values as defaults."""
        return (
            self.start if start is None else start,
            self.end if end is None else end,
            self.orientation if orientation is None else orientation
        )

    def slice(self, start: int | None = None, end: int | None = None, 
              orientation: Orientation = Orientation.FORWARD, seq: T = None) -> T:
        if seq is None:
            raise ValueError("seq is required")
        start, end, orientation = self._resolve_params(start, end, orientation)
        return super().slice(start, end, orientation, seq)

    def normalize_range(self, start: int | None = None, end: int | None = None, 
                        orientation: Orientation = Orientation.FORWARD) -> tuple[int, int, bool | None]:
        start, end, orientation = self._resolve_params(start, end, orientation)
        return super().normalize_range(start, end, orientation)

    def get_segment_length(self, start: int | None = None, end: int | None = None, 
                           orientation: Orientation = Orientation.FORWARD) -> int:
        start, end, orientation = self._resolve_params(start, end, orientation)
        return super().get_segment_length(start, end, orientation)

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

    def slice(self, start: int, end: int, orientation: Orientation = Orientation.FORWARD, seq: T | None = None) -> T:
        """Get range (1-based inclusive) from sequence."""
        seq = self.seq if seq is None else seq
        return super().slice(start, end, orientation, seq)


@dataclass
class SeqSegmentScaffold(SegmentMixin, SeqScaffold):
    """Sequence scaffold with persisted start/end."""
    def slice(self, start: int | None = None, end: int | None = None, 
              orientation: Orientation = Orientation.FORWARD, seq: T | None = None) -> T:
        seq = self.seq if seq is None else seq
        return super().slice(start, end, orientation, seq)

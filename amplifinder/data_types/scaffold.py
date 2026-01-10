"""
Scaffold coordinate system:

Scaffold is a circular or linear sequence of DNA, with a length and a circularity flag.

Its coordinates are 1-based inclusive.

`start` is the position where the sequence starts, and
`end` is the position where the sequence ends.
`orientation` is the orientation of the sequence, specifying the FORWARD or REVERSE strand.

Therefore start/end are NOT left/right.

For FORWARD orientation:
* start <= end: normal segment, span_origin is False

<-------------------length------------------>
~~~~~~~~~~~~|======================|~~~~~~~~~
            start                   end
            |                      |
            ----------------------->           inclusive
                 segment_length

* start > end and circular scaffold => span_origin is True

============|~~~~~~~~~~~~~~~~~~~~~~|=========
            end                   start
            |                      |
            |                      --------->  inclusive
------------>                      |
                                   --------->------------>
                                       segment_length


* if start > end and scaffold is linear, the segment is not valid => span_origin is None

For REVERSE orientation:
- start >= end: normal reverse-complemented segment, span_origin is False

~~~~~~~~~~~~|======================|~~~~~~~~~
            end                   start
            |                      |
            <-----------------------           inclusive
                 segment_length

* start < end and scaffold is circular => span_origin is True

============|~~~~~~~~~~~~~~~~~~~~~~|=========
          start                   end
            |                      |
<------------                      |
                                   <---------
"""

from __future__ import annotations

import numpy as np

from typing import ClassVar, Optional, TypeVar
from pydantic import field_validator, model_validator

from Bio.Seq import reverse_complement

from amplifinder.data_types.records import Record
from amplifinder.data_types.enums import Orientation
from amplifinder.utils.sequence_utils import concatenate_sequence

T = TypeVar('T', str, np.ndarray)


class JcArm(Record):
    """Junction arm coordinates and orientation."""
    NAME: ClassVar[str] = "Junction arms"
    scaf: str
    start: int
    dir: Orientation
    flank: int

    @field_validator('dir')
    @classmethod
    def validate_dir(cls, v: Orientation) -> Orientation:
        if v not in [Orientation.FORWARD, Orientation.REVERSE]:
            raise ValueError("Invalid orientation")
        return v

    @property
    def end(self) -> int:
        """Compute end position based on start position, flank length, and direction."""
        return self.start + (self.flank - 1) * self.dir

    def shift_by_offset(self, offset: int) -> 'JcArm':
        """Shift arm by offset in the direction of the arm.

        Args:
            offset: Distance to shift (positive = in arm direction, negative = opposite)

        Returns:
            New JcArm with shifted start position
        """
        return JcArm(
            scaf=self.scaf,
            start=self.start + offset * self.dir,
            dir=self.dir,
            flank=self.flank
        )

    def mirror(self) -> 'JcArm':
        """Create a mirror arm at the adjacent position.
        ~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~
                           |------>           self
                    <------|                  mirror
        """
        return JcArm(
            scaf=self.scaf,
            start=self.start - self.dir,
            dir=self.dir.opposite(),
            flank=self.flank
        )


class Scaffold(Record):
    """Coordinate helper with circularity awareness."""
    NAME: ClassVar[str] = "Scaffolds"
    is_circular: bool
    length: int
    scaf: Optional[str] = None  # Scaffold name/identifier

    def __len__(self) -> int:
        return self.length

    def circular_normalize_range(
        self, start: int, end: int,
        orientation: Orientation = Orientation.FORWARD
    ) -> tuple[int, int, bool | None]:
        """Normalize coordinates for range extraction (1-based inclusive).
        Returns:
            start_norm: Circularly normalized start coordinate (modulo)
            end_norm: Circularly normalized end coordinate (modulo)
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
        start, end, span_origin = self.circular_normalize_range(start, end, orientation)
        if span_origin is None:
            return -1  # invalid segment (span origin, but scaffold is not circular)
        diff = end - start if orientation == Orientation.FORWARD else start - end
        # diff is negative if we span the circular origin:
        assert diff < 0 and span_origin or diff > 0 and not span_origin
        if span_origin:
            return self.length + diff + 1
        return diff + 1

    def slice(self, start: int, end: int, orientation: Orientation, seq: T) -> T:
        """Slice arbitrary scaffold-aligned sequence (string or ndarray)."""
        assert len(seq) == self.length
        start, end, spans_origin = self.circular_normalize_range(start, end, orientation)
        if spans_origin is None:
            segment = seq[0:0]  # invalid segment
        if orientation == Orientation.FORWARD:
            if spans_origin:
                segment = concatenate_sequence(seq[start - 1:], seq[:end])
            else:
                segment = seq[start - 1:end]
        else:
            if spans_origin:
                segment = concatenate_sequence(seq[end - 1:], seq[:start])
            else:
                segment = seq[end - 1:start]
            segment = reverse_complement(segment)
        return segment


class SeqScaffold(Scaffold):
    """Sequence-backed scaffold."""
    model_config = {"arbitrary_types_allowed": True}

    seq: T
    length: int = 0  # Will be computed from seq

    @model_validator(mode='after')
    def compute_length(self):
        """Compute length from sequence."""
        object.__setattr__(self, 'length', len(self.seq))
        return self

    def slice(self, start: int, end: int, orientation: Orientation = Orientation.FORWARD, seq: T | None = None) -> T:
        """Get range (1-based inclusive) from sequence."""
        seq = self.seq if seq is None else seq
        return super().slice(start, end, orientation, seq)


class SegmentMixin:
    """Mixin adding start/end helpers."""

    start: int
    end: int
    orientation: Orientation = Orientation.FORWARD

    @classmethod
    def from_scaffold_and_jc_arm(cls, scaffold: Scaffold | SeqScaffold, jc_arm: JcArm):
        """Create SegmentScaffold from JcArm.

        Args:
            scaffold: Scaffold (or SeqScaffold) providing circularity, length, and optionally seq
            jc_arm: Junction arm with start, direction, and flank

        Returns:
            SegmentScaffold or SeqSegmentScaffold with JcArm coordinates
        """
        assert scaffold.scaf == jc_arm.scaf
        return cls.from_other(
            scaffold,
            start=jc_arm.start,
            end=jc_arm.end,
            orientation=jc_arm.dir
        )

    @classmethod
    def from_scaffold_left_right_orientation(
        cls, scaffold: Scaffold | SeqScaffold, left: int, right: int,
        orientation: Orientation, **kwargs
    ) -> SegmentScaffold | SeqSegmentScaffold:
        """Create SegmentScaffold from scaffold, left, right, and orientation.

        Additional kwargs are passed to from_other for subclass-specific fields.
        """
        if orientation == Orientation.FORWARD:
            start, end = left, right
        else:
            start, end = right, left
        return cls.from_other(
            scaffold,
            start=start,
            end=end,
            orientation=orientation,
            **kwargs
        )

    def _resolve_params(self, start: int | None, end: int | None,
                        orientation: Orientation | None) -> tuple[int, int, Orientation]:
        """Resolve parameters, using instance values as defaults."""
        return (
            self.start if start is None else start,
            self.end if end is None else end,
            self.orientation if orientation is None else orientation
        )

    def slice(self, start: int | None = None, end: int | None = None,
              orientation: Orientation | None = None, seq: T | None = None) -> T:
        start, end, orientation = self._resolve_params(start, end, orientation)
        return super().slice(start, end, orientation, seq)

    def circular_normalize_range(
        self, start: int | None = None, end: int | None = None,
        orientation: Orientation | None = None
    ) -> tuple[int, int, bool | None]:
        start, end, orientation = self._resolve_params(start, end, orientation)
        return super().circular_normalize_range(start, end, orientation)

    def get_segment_length(self, start: int | None = None, end: int | None = None,
                           orientation: Orientation | None = None) -> int:
        start, end, orientation = self._resolve_params(start, end, orientation)
        return super().get_segment_length(start, end, orientation)

    @property
    def segment_length(self) -> int:
        return self.get_segment_length()

    @property
    def span_origin(self) -> bool | None:
        return self.circular_normalize_range()[2]

    @property
    def left(self) -> int:
        """Leftmost position on the forward strand."""
        return self.start if self.orientation == Orientation.FORWARD else self.end

    @property
    def right(self) -> int:
        """Rightmost position on the forward strand."""
        return self.end if self.orientation == Orientation.FORWARD else self.start

    def get_inward_arms(self, flanks: int | tuple[int, int]) -> tuple[JcArm, JcArm]:
        """Creates inward pointing arms at the segment start and end positions.

        Args:
            flanks:  Flank length(s) for arms. Single int applies to both arms,
                     tuple of (start_flank, end_flank) for different lengths.

        Returns:
            Tuple of (arm_start, arm_end)
        """
        if isinstance(flanks, int):
            flanks = (flanks, flanks)
        start_flank, end_flank = flanks
        arm_S = JcArm(scaf=self.scaf, start=self.start, dir=self.orientation, flank=start_flank)
        arm_E = JcArm(scaf=self.scaf, start=self.end, dir=self.orientation.opposite(), flank=end_flank)
        return arm_S, arm_E

    def get_outward_arms(self, flanks: int | tuple[int, int]) -> tuple[JcArm, JcArm]:
        """Creates outward pointing arms at the segment start and end positions.
        Args, as in get_inward_arms
        """
        arm_S, arm_E = self.get_inward_arms(flanks=flanks)
        return arm_S.mirror(), arm_E.mirror()

    def get_junctions(self, out_flanks: int | tuple[int, int],
                      in_flanks: int | tuple[int, int],
                      junction_class=None):
        """Creates junctions at the segment start and end positions,
        In each junction arm1 is inward and arm2 is outward.
        """
        if junction_class is None:
            from amplifinder.data_types.record_types import Junction
            junction_class = Junction

        inward_arms = self.get_inward_arms(flanks=in_flanks)
        outward_arms = self.get_outward_arms(flanks=out_flanks)

        return (
            junction_class.from_jc_arms(inward_arms[0], outward_arms[0]),
            junction_class.from_jc_arms(inward_arms[1], outward_arms[1]),
        )


class SegmentScaffold(SegmentMixin, Scaffold):
    """Scaffold with persisted start/end."""
    pass


class SeqSegmentScaffold(SegmentMixin, SeqScaffold):
    """Sequence scaffold with persisted start/end."""
    pass

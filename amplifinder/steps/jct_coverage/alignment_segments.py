"""Convert CIGAR alignments to plottable coordinate segments."""
from __future__ import annotations
from typing import Iterator, TypeVar, Generic, Callable
from dataclasses import dataclass, fields

from numpy import ndarray, nan, full, column_stack

from amplifinder.steps.jct_coverage.alignment_data import AlignmentData, PairedAlignment, BaseSingleAlignment
from amplifinder.steps.jct_coverage.cigar import Cigar


# Plot hits coordinate system:
# start, end are 0-based end-exclusive coordinates.
#
#  |---|---|---|---|---|---|---|---|
# -4  -3  -2  -1   0   1   2   3   4  x-axis
#    0   1   2   3   4   5   6   7    seq index
#
# junction_length = 8, arm_len = 4
# read with start=3, end=7 will be plotted as (len = 4)
#              |---|---|---|---|
#             -1               3      x-axis
# x_start = start - arm_len = 3 - 4 = -1
# x_end = end - arm_len = 7 - 4 = 3


# CIGAR operation codes (pysam):
# 0=M (match/mismatch), 1=I (insertion), 2=D (deletion),
# 7== (sequence match), 8=X (sequence mismatch)
MATCH_OPS = {0, 7}
MISMATCH_OPS = {8}
INSERT_OPS = {1}
DELETE_OPS = {2}

SNP_VLINE_HALF_HEIGHT = 0.18


class Coords:
    def __init__(self, x: list = None, y: list = None):
        self.x: list = x or []
        self.y: list = y or []

    def append(self, x, y) -> None:
        self.x.append(x)
        self.y.append(y)

    def extend(self, other: "Coords") -> None:
        self.x.extend(other.x)
        self.y.extend(other.y)


def convert_coords_to_nan_separated_arrays(coords: "Coords") -> tuple[ndarray, ndarray]:
    x = column_stack([coords.x, full(len(coords.x), nan)]).ravel()
    y = column_stack([coords.y, full(len(coords.y), nan)]).ravel()
    return x, y


T = TypeVar('T')


@dataclass
class AlignmentElements(Generic[T]):
    match: T
    mismatch: T
    deletion: T
    insertion: T
    snp: T

    def items(self) -> Iterator[tuple[str, T]]:
        """Iterate over (name, value) pairs."""
        for field in fields(self):
            yield field.name, getattr(self, field.name)

    def __iter__(self) -> Iterator[T]:
        """Iterate over field values."""
        for _, value in self.items():
            yield value

    def __getitem__(self, key: str) -> T:
        return getattr(self, key)

    def apply_in_place(self, func: Callable[[T], None]) -> None:
        for value in self:
            func(value)

    def apply(self, func: Callable[[T], T]) -> "AlignmentElements[T]":
        return type(self)(*[func(value) for value in self])


class CoordsAlignmentElements(AlignmentElements[Coords]):
    @classmethod
    def create(cls) -> "CoordsAlignmentElements":
        return cls(*[Coords() for _ in fields(cls)])

    def extend(self, other: "CoordsAlignmentElements") -> None:
        """Extend segments by adding segments from other."""
        for self_coord, other_coord in zip(self, other):
            self_coord.extend(other_coord)


def get_alignment_segments_from_cigar(
    cigar: Cigar,
    start: int,
    x0: int = 0,
    y0: int = 0,
) -> CoordsAlignmentElements:
    """Extract segments from a single CIGAR with detailed events."""
    segments = CoordsAlignmentElements.create()

    for op, length, x_start in cigar.iter_with_ref_pos(start):
        x_end = x_start + length

        if op in INSERT_OPS:
            segments.insertion.append(x_start + x0, y0)
            continue

        segment = ([x_start + x0, x_end + x0], [y0, y0])

        if op in MATCH_OPS:
            segments.match.append(*segment)
        elif op in DELETE_OPS:
            segments.deletion.append(*segment)
        elif op in MISMATCH_OPS:
            segments.mismatch.append(*segment)
            for i in range(length):
                x_pos = x_start + i + 0.5
                segments.snp.append([x_pos + x0, x_pos + x0], [y0 - SNP_VLINE_HALF_HEIGHT, y0 + SNP_VLINE_HALF_HEIGHT])
        else:
            raise ValueError(f"Unsupported CIGAR operation: {op}")

    return segments


def get_alignment_segments_from_single_alignment(
    alignment: BaseSingleAlignment,
    show_events: bool,
    x0: int = 0,
    y0: int = 0,
) -> CoordsAlignmentElements:
    """Extract segments from a single alignment."""
    segments = CoordsAlignmentElements.create()

    if not show_events:
        segments.match.append([alignment.left + x0, alignment.right + x0], [y0, y0])
        return segments

    # Get all cigars with their start positions
    for start, cigar in alignment.get_starts_and_cigars():
        segments.extend(get_alignment_segments_from_cigar(cigar, start, x0, y0))

    return segments


def get_alignment_segments(
    alignment: AlignmentData,
    show_events: bool,
    x0: int = 0,
    y0: int = 0,
) -> CoordsAlignmentElements:
    """Get segments for an AlignmentData with SNP/indel annotations."""
    if isinstance(alignment, PairedAlignment):
        segments1 = get_alignment_segments_from_single_alignment(alignment.forward_alignment, show_events, x0, y0)
        segments2 = get_alignment_segments_from_single_alignment(alignment.reverse_alignment, show_events, x0, y0)
        segments1.extend(segments2)
        return segments1
    assert isinstance(alignment, BaseSingleAlignment)
    return get_alignment_segments_from_single_alignment(alignment, show_events, x0, y0)

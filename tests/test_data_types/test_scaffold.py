import numpy as np
import pytest

from amplifinder.data_types.enums import Orientation
from amplifinder.data_types.scaffold import Scaffold, SeqScaffold, SegmentScaffold, SeqSegmentScaffold


def test_seqscaffold_slice_forward_linear():
    scaf = SeqScaffold(seq="abcdef", is_circular=False)
    assert scaf.slice(2, 4) == "bcd"


def test_seqscaffold_slice_reverse_circular():
    scaf = SeqScaffold(seq="ATGCAT", is_circular=True)
    assert scaf.slice(5, 2) == "ATAT"


def test_scaffold_slice_numpy_wrap():
    scaf = Scaffold(is_circular=True, length=5)
    arr = np.arange(1, 6)
    np.testing.assert_array_equal(scaf.slice(4, 2, seq=arr), [4, 5, 1, 2])


def test_segment_scaffold_defaults_start_end():
    seg_scaf = SegmentScaffold(is_circular=False, length=6, start=2, end=4)
    arr = np.arange(6)
    np.testing.assert_array_equal(seg_scaf.slice(seq=arr), np.array([1, 2, 3]))


def test_seqsegment_scaffold_defaults_seq_and_coords():
    seg_scaf = SeqSegmentScaffold(seq="abcdef", is_circular=False, start=2, end=4)
    assert seg_scaf.slice() == "bcd"


if __name__ == "__main__":
    test_seqsegment_scaffold_defaults_seq_and_coords()
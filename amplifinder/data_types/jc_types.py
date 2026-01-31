"""Junction types and element enums for AmpliFinder."""
from __future__ import annotations
import numpy as np

from enum import Enum
from typing import TypeAlias
from functools import wraps

from amplifinder.data_types.basic_enums import Side


class Element(str, Enum):
    """Elements (chromosome, cassette, IS)."""
    CHR = "chr"  # chromosome
    AMP = "amp"  # cassette
    TN = "tn"    # IS


class JunctionType(Enum):
    """The 7 synthetic junction types.

    Amplicon structure:

    ~~~~~~~~~~> |---> |========> |---> |========> |---> |~~~~~~~~~~
    chromosome   IS    amplicon   IS    amplicon   IS    chromosome

    """
    CHR_AMP = (1, '~~> |==', Element.CHR, Element.AMP, -1)  # noqa: E221, E203
    CHR_TN  = (2, '~~> |--', Element.CHR, Element.TN , -2)  # noqa: E221, E203
    AMP_TN  = (3, '==> |--', Element.AMP, Element.TN , +3)  # noqa: E221, E203
    AMP_AMP = (4, '==> |==', Element.AMP, Element.AMP,  0)  # noqa: E221, E203
    TN_AMP  = (5, '--> |==', Element.TN,  Element.AMP, -3)  # noqa: E221, E203
    TN_CHR  = (6, '--> |~~', Element.TN,  Element.CHR, +2)  # noqa: E221, E203
    AMP_CHR = (7, '==> |~~', Element.AMP, Element.CHR, +1)  # noqa: E221, E203

    @classmethod
    def from_elements(cls, element_left: Element, element_right: Element):
        for jc in cls:
            if jc.element_pair == (element_left, element_right):
                return jc

    @classmethod
    def get_element_pairs_to_jcs(cls, side: Side) -> dict[tuple[Element, Element], JunctionType]:
        pairs_to_jcs: dict[tuple[Element, Element], JunctionType] = {}
        for jc in JunctionType:
            pair = jc.element_pair
            left, right = pair
            if side == Side.LEFT:
                pairs_to_jcs[pair] = cls.get_element_pairs_to_jcs(left, right)
            elif side == Side.RIGHT:
                pairs_to_jcs[pair] = cls.get_element_pairs_to_jcs(right, left)
        return pairs_to_jcs

    @property
    def num(self) -> int:
        """Return the number of the junction. (1-7, matlab legacy)"""
        return self.value[0]

    @property
    def symbol(self) -> str:
        """Return the symbol of the junction."""
        return self.value[1]

    @property
    def element_pair(self) -> tuple[Element, Element]:
        """Return the elements of the junction."""
        return self.value[2], self.value[3]

    @property
    def order(self) -> int:
        """Return the number of the junction."""
        return self.value[4]

    @property
    def side(self) -> Side:
        """Return the side of the junction."""
        return Side(np.sign(self.order))


class JcCall(Enum):

    POS = True
    NEG = False
    AMBIGIOUS = None

    """Junction call: True if the junction is called, False if not, None if not enough evidence."""
    value: bool | None

    def __bool__(self):
        return bool(self.value)

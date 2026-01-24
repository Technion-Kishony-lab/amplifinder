"""Junction types and element enums for AmpliFinder."""
from __future__ import annotations
import numpy as np

from enum import Enum
from typing import TypeAlias

from amplifinder.data_types.basic_enums import Side


class Element(str, Enum):
    """Elements (chromosome, cassette, IS)."""
    CHR = "chr"  # chromosome
    AMP = "amp"  # cassette
    TN = "tn"    # IS


class JunctionType(Enum):
    """The 7 synthetic junction types.

    Amplicon structure:

    ~~~~~~~~~>>>======>>>======>>>~~~~~~~~~

    legend:
    ~~~    chromosome
    >>>    IS
    ====== cassette (amplicon)

    """
    CHR_TO_AMP_LEFT       = (1, '~~-==', Element.CHR, Element.AMP, -1)  # noqa: E221, E203
    CHR_TO_TN_LEFT        = (2, '~~->>', Element.CHR, Element.TN , -2)  # noqa: E221, E203
    AMP_RIGHT_TO_TN_LEFT  = (3, '==->>', Element.AMP, Element.TN , +3)  # noqa: E221, E203
    AMP_RIGHT_TO_AMP_LEFT = (4, '==-==', Element.AMP, Element.AMP,  0)  # noqa: E221, E203
    TN_RIGHT_TO_AMP_LEFT  = (5, '>>-==', Element.TN,  Element.AMP, -3)  # noqa: E221, E203
    TN_RIGHT_TO_CHR       = (6, '>>-~~', Element.TN,  Element.CHR, +2)  # noqa: E221, E203
    AMP_RIGHT_TO_CHR      = (7, '==-~~', Element.AMP, Element.CHR, +1)  # noqa: E221, E203

    @classmethod
    def sorted(cls) -> list[JunctionType]:
        """Return the JunctionType enum members sorted by value."""
        return sorted(cls, key=lambda x: x.num)

    @property
    def num(self) -> int:
        """Return the number of the junction. (1-7, matlab legacy)"""
        return self.value[0]

    @property
    def symbol(self) -> str:
        """Return the symbol of the junction."""
        return self.value[1]

    @property
    def elements(self) -> tuple[Element, Element]:
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


JcCall: TypeAlias = bool | None
"""Junction call: True if the junction is called, False if not, None if not enough evidence."""

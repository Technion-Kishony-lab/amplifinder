"""Junction types for AmpliFinder."""
from __future__ import annotations

from typing import Any, ClassVar, Dict, List, Optional, TypeVar
from pydantic import ConfigDict, field_validator

from amplifinder.records.base_records import Record
from amplifinder.data_types.basic_enums import Orientation


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

    def shift_by_offset(self, offset: int) -> JcArm:
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
    
    def get_distance_to(self, pos: int) -> int:
        """Get distance to a position.
        ~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~
                           |----arm----> 
                           |        | 
                         start     pos
                           |--dist->|
        return: 
             0 if pos is the arm's start position
            >0 if pos is further in the arm's direction
            <0 if pos is further in the opposite direction
        """
        return pos - self.start if self.dir == Orientation.FORWARD else self.start - pos

    def mirror(self) -> JcArm:
        """Create a mirror arm at the adjacent position.
        ~~~~~~~~~~~~~~~~~~~|~~~~~~~~~~~~~
                           |------> 
                            self
                    <------|                  
                     mirror
        """
        return JcArm(
            scaf=self.scaf,
            start=self.start - self.dir,
            dir=self.dir.opposite(),
            flank=self.flank
        )


JunctionT = TypeVar("JunctionT", bound="Junction")


class Junction(Record):
    """Base junction record with coordinate fields only."""
    NAME: ClassVar[str] = "Junctions"

    # Arm 1 fields
    scaf1: str
    pos1: int
    dir1: Orientation
    flanking1: int   # Length of sequence flanking arm 1 (used for sequence extraction)

    # Arm 2 fields
    scaf2: str
    pos2: int
    dir2: Orientation
    flanking2: int  # Length of sequence flanking arm 2 (used for sequence extraction)

    @classmethod
    def _get_extra_fields(cls) -> List[str]:
        """Get extra fields not part of the arm coordinates."""
        arm_fields = Junction.model_fields.keys()
        return [f for f in cls.model_fields if f not in arm_fields]

    def _get_extra_kwargs(self) -> Dict[str, Any]:
        """Get extra kwargs not part of the arm coordinates."""
        return {field: getattr(self, field) for field in self._get_extra_fields()}

    def swap_sides(self: JunctionT, **kwargs) -> JunctionT:
        """Return new junction with arm 1 and arm 2 swapped."""
        jc_arms = self.get_jc_arms()
        extra_kwargs = self._get_extra_kwargs()
        extra_kwargs.update(kwargs)
        return self.from_jc_arms(jc_arms[1], jc_arms[0], **extra_kwargs)

    def get_jc_arm(self, arm: int) -> JcArm:
        """Get scaffold, position, direction, and flanking length for an arm."""
        return self.get_jc_arms()[arm - 1]

    def get_jc_arms(self) -> tuple[JcArm, JcArm]:
        """Get junction arms."""
        return (
            JcArm(scaf=self.scaf1, start=self.pos1, dir=self.dir1, flank=self.flanking1),
            JcArm(scaf=self.scaf2, start=self.pos2, dir=self.dir2, flank=self.flanking2),
        )

    @classmethod
    def from_jc_arms(cls, arm1: JcArm, arm2: JcArm, **kwargs) -> Junction:
        """Create a Junction from junction arm coordinates."""
        return cls(scaf1=arm1.scaf, pos1=arm1.start, dir1=arm1.dir, flanking1=arm1.flank,
                   scaf2=arm2.scaf, pos2=arm2.start, dir2=arm2.dir, flanking2=arm2.flank,
                   **kwargs)


class NumJunction(Junction):
    """Junction with identifier."""
    NAME: ClassVar[str] = "Numbered Junctions"
    # Junction identifier: breseq junction number (positive), or negative for reference junctions
    num: Optional[int] = None


class BreseqJunction(NumJunction):
    """Breseq junction."""
    NAME: ClassVar[str] = "Breseq junctions"
    model_config = ConfigDict(extra='allow')
    ALLOW_EXTRA: ClassVar[bool] = True

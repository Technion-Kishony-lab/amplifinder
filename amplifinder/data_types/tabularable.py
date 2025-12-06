"""Base class for dataclasses that can be converted to/from DataFrames."""

from __future__ import annotations

from dataclasses import dataclass, field, fields
from typing import Any, ClassVar, Dict, List, Set, Type, TypeVar

import pandas as pd

T = TypeVar("T", bound="Tabularable")


@dataclass
class Tabularable:
    """Base dataclass with an 'extra' field for unknown columns.

    Provides to_dict/from_dict and DataFrame conversion utilities.
    """

    # Fields to ignore when converting from dict (override in subclass)
    _ignored_fields: ClassVar[Set[str]] = set()

    # All other fields stored here
    extra: Dict[str, Any] = field(default_factory=dict)

    @classmethod
    def _field_names(cls) -> List[str]:
        """Get names of dataclass fields (excluding 'extra')."""
        return [f.name for f in fields(cls) if f.name != "extra"]

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        d = {name: getattr(self, name) for name in self._field_names()}
        d.update(self.extra)
        return d

    @classmethod
    def from_dict(cls: Type[T], d: Dict[str, Any]) -> T:
        """Create instance from dictionary."""
        known = set(cls._field_names())
        kwargs = {k: v for k, v in d.items() if k in known}
        extra = {k: v for k, v in d.items() if k not in known and k not in cls._ignored_fields}
        return cls(**kwargs, extra=extra)

    @classmethod
    def list_to_dataframe(cls: Type[T], items: List[T]) -> pd.DataFrame:
        """Convert list of Tabularable objects to DataFrame."""
        if not items:
            return pd.DataFrame(columns=cls._field_names())
        return pd.DataFrame([item.to_dict() for item in items])

    @classmethod
    def dataframe_to_list(cls: Type[T], df: pd.DataFrame) -> List[T]:
        """Convert DataFrame to list of Tabularable objects."""
        if df.empty:
            return []
        return [cls.from_dict(row.to_dict()) for _, row in df.iterrows()]

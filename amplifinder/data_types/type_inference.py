"""Type inference utilities for Record properties and fields."""
from typing import Any, Type, get_origin, get_args, Union, get_type_hints


def infer_property_type(cls: Type, property_name: str) -> tuple[Type, bool]:
    """Infer type and optionality from a property's return annotation.

    Args:
        cls: The class containing the property
        property_name: Name of the property

    Returns:
        Tuple of (dtype, is_optional) where:
        - dtype: The inferred type (or Any if inference fails)
        - is_optional: True if type is Optional[T]
    """
    try:
        prop = getattr(cls, property_name)
        if isinstance(prop, property) and prop.fget is not None:
            hints = get_type_hints(prop.fget)
            dtype = hints.get('return', Any)
            is_optional = (get_origin(dtype) is Union and type(None) in get_args(dtype))
            return dtype, is_optional
        else:
            return Any, True
    except (AttributeError, TypeError):
        # Fallback if property doesn't exist or can't get hints
        return Any, True

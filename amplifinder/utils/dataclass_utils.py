"""Utility functions for dataclass field type conversions."""

from dataclasses import fields
from typing import Any, Dict, Type, TypeVar, Optional, get_args, get_origin

T = TypeVar('T')


TRUE_VALUES = ("1", "true", "yes", "y", "t")
FALSE_VALUES = ("0", "false", "no", "n", "f")
NONE_VALUES = ("", "na", "none", "null")


def str_to_bool(value: str) -> Optional[bool]:
    """Convert string to boolean.
    
    Args:
        value: String value to convert
        
    Returns:
        True if value is in TRUE_VALUES, False if value is in FALSE_VALUES, 
        None if value is in NONE_VALUES
        
    Raises:
        ValueError: If value is not a recognized boolean string
    """
    value_lower = value.lower()
    if value_lower in TRUE_VALUES:
        return True
    elif value_lower in FALSE_VALUES:
        return False
    elif value_lower in NONE_VALUES:
        return None
    raise ValueError(f"Invalid boolean value: {value}")


def get_base_type(field_type: Type) -> Optional[Type]:
    """Extract the base type from a type annotation, unwrapping Optional.
    
    Args:
        field_type: Type annotation (e.g., int, Optional[int], str)
        
    Returns:
        Base type (e.g., int, str, bool) or None if cannot determine
    """
    # Handle Optional[T] which is Union[T, None]
    origin = get_origin(field_type)
    if origin is not None:
        args = get_args(field_type)
        # For Optional[T], args is (T, type(None))
        if len(args) == 2 and type(None) in args:
            return next(arg for arg in args if arg is not type(None))
    
    # Direct type like int, str, bool
    return field_type if field_type in (int, float, bool, str) else None


def convert_value(value: str, target_type: Type) -> Any:
    """Convert a string value to the target type.
    
    Args:
        value: String value to convert
        target_type: Target type (int, float, bool, str, etc.)
        
    Returns:
        Converted value
        
    Raises:
        ValueError: If conversion fails
    """
    # Special handling for bool (custom conversion logic)
    if target_type is bool:
        return str_to_bool(value)
    # For all other types, use the type constructor
    return target_type(value)


def get_field_types(cls: Type[T]) -> Dict[str, Type]:
    """Get a mapping of field names to their base types for a dataclass.
    
    Args:
        cls: Dataclass type to inspect
        
    Returns:
        Dictionary mapping field names to base types (unwrapped from Optional)
    """
    field_types = {}
    for f in fields(cls):
        base_type = get_base_type(f.type)
        if base_type is not None:
            field_types[f.name] = base_type
    return field_types


def convert_csv_row_types(row_args: Dict[str, Any], cls: Type[T]) -> None:
    """Convert CSV string values to appropriate types for a dataclass.
    
    Automatically converts string values to int, float, bool, etc. based on
    the dataclass field types. Modifies row_args in place.
    
    Args:
        row_args: Dictionary of field names to values (from CSV row)
        cls: Dataclass type to use for type information
    """
    field_types = get_field_types(cls)
    
    # Convert string values to appropriate types
    for key, value in list(row_args.items()):
        if key in field_types and isinstance(value, str):
            try:
                row_args[key] = convert_value(value, field_types[key])
            except (ValueError, TypeError) as e:
                raise ValueError(f"Error converting value for field {key}: {e}") from e

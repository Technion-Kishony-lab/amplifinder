"""YAML/JSON loading and annotated-YAML writing utilities."""

import dataclasses
import json
from dataclasses import fields as dc_fields
from pathlib import Path
from typing import Any, Optional, Sequence, get_args

import yaml


def to_commented_map(d: dict):
    """Recursively convert dict to ruamel.yaml CommentedMap with flow-style lists."""
    from ruamel.yaml.comments import CommentedMap, CommentedSeq
    cm = CommentedMap()
    for k, v in d.items():
        if isinstance(v, dict):
            cm[k] = to_commented_map(v)
        elif isinstance(v, list):
            seq = CommentedSeq(v)
            seq.fa.set_flow_style()
            cm[k] = seq
        else:
            cm[k] = v
    return cm


def _get_nested_dataclass(field_type) -> Optional[type]:
    """Return the dataclass type if *field_type* is (Optional[]) a dataclass, else None."""
    if dataclasses.is_dataclass(field_type):
        return field_type
    for arg in get_args(field_type):
        if dataclasses.is_dataclass(arg):
            return arg
    return None


def dataclass_to_commented_yaml(data: dict[str, Any],
                                dc_field_list: Sequence,
                                header: Optional[list[str]] = None):
    """Build a ruamel.yaml CommentedMap from a data dict, annotated via dataclass field metadata.

    Fields may carry ``metadata={"section": "Title", "comment": "inline note"}``.
    A new ``# --- Title ---`` header is emitted whenever *section* changes.
    Nested dataclass fields are recursively annotated.

    Args:
        data: Flat YAML-serializable dict (e.g. from ``to_yaml_dict()``).
        dc_field_list: Sequence of dataclass fields (``dataclasses.fields(instance)``).
        header: Optional header comment lines (without ``#`` prefix).
    """
    from ruamel.yaml.comments import CommentedMap, CommentedSeq

    cm = CommentedMap()
    is_first_section = True
    current_section = None

    for f in dc_field_list:
        if f.name not in data:
            continue

        value = data[f.name]
        if isinstance(value, dict):
            nested_dc = _get_nested_dataclass(f.type)
            if nested_dc is not None:
                cm[f.name] = dataclass_to_commented_yaml(value, dc_fields(nested_dc))
            else:
                cm[f.name] = to_commented_map(value)
        elif isinstance(value, list):
            seq = CommentedSeq(value)
            seq.fa.set_flow_style()
            cm[f.name] = seq
        else:
            cm[f.name] = value

        meta = f.metadata
        section = meta.get("section") if meta else None
        comment = meta.get("comment") if meta else None

        if section and section != current_section:
            parts = []
            if is_first_section and header:
                parts.extend(header)
                parts.append('')
            parts.append(f'--- {section} ---')
            before_text = '\n'.join(parts)
            if not is_first_section:
                before_text = '\n' + before_text
            cm.yaml_set_comment_before_after_key(f.name, before=before_text, indent=0)
            current_section = section
            is_first_section = False
        elif isinstance(value, dict):
            cm.yaml_set_comment_before_after_key(f.name, before='\n', indent=0)

        if comment:
            cm.yaml_add_eol_comment(comment, f.name)

    return cm


def save_annotated_yaml(data: dict[str, Any],
                        dc_field_list: Sequence,
                        path: Path,
                        header: Optional[list[str]] = None) -> None:
    """Build an annotated CommentedMap from dataclass metadata and write it to *path*."""
    cm = dataclass_to_commented_yaml(data, dc_field_list, header)
    dump_yaml(cm, path)


def dump_yaml(data, path: Path) -> None:
    """Dump a CommentedMap (or plain dict) to a YAML file via ruamel.yaml."""
    from ruamel.yaml import YAML as RuamelYAML

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    yw = RuamelYAML()
    yw.default_flow_style = False
    yw.representer.add_representer(
        type(None),
        lambda dumper, _: dumper.represent_scalar('tag:yaml.org,2002:null', 'null'),
    )
    with open(path, 'w') as f:
        yw.dump(data, f)


def load_config(config_path: Path) -> dict[str, Any]:
    """Load configuration from YAML or JSON file.

    Args:
        config_path: Path to config file (.yaml, .yml, or .json)

    Returns:
        Dictionary of configuration values
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    suffix = config_path.suffix.lower()

    with open(config_path) as f:
        if suffix in (".yaml", ".yml"):
            return yaml.safe_load(f) or {}
        elif suffix == ".json":
            return json.load(f)
        else:
            raise ValueError(
                f"Unsupported config format: {suffix}. "
                "Use .yaml, .yml, or .json"
            )

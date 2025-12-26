"""Helpers for default filenames / paths for Record-based tables."""

from pathlib import Path
from typing import Dict, Type

from amplifinder.data_types import (
    RefTnJunction, TnJunction, RawTnJc2, CoveredTnJc2, ClassifiedTnJc2,
    FilteredTnJc2, AnalyzedTnJc2, ExportedTnJc2,
)
from amplifinder.data_types.records import Record


DEFAULT_FILENAMES: Dict[Type[Record], str] = {
    RefTnJunction: "ref_tn_jc.csv",
    TnJunction: "tn_jc.csv",
    RawTnJc2: "tnjc2_raw.csv",
    CoveredTnJc2: "tnjc2_covered.csv",
    ClassifiedTnJc2: "tnjc2_classified.csv",
    FilteredTnJc2: "tnjc2_filtered.csv",
    AnalyzedTnJc2: "tnjc2_analyzed.csv",
    ExportedTnJc2: "tnjc2_exported.csv",
}


def default_filename(record_type: Type[Record]) -> str:
    return DEFAULT_FILENAMES[record_type]


def default_path(output_dir: Path, record_type: Type[Record]) -> Path:
    return output_dir / default_filename(record_type)

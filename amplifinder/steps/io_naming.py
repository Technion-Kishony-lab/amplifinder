"""Helpers for default filenames / paths for Record-based tables."""

from pathlib import Path
from typing import Dict, Type

from amplifinder.data_types import (
    RefTnJunction, TnJunction, RawTnJc2, CoveredTnJc2, SingleLocusLinkedTnJc2,
    SynJctsTnJc2, AnalyzedTnJc2, ExportedTnJc2,
)
from amplifinder.data_types.records import Record


DEFAULT_FILENAMES: Dict[Type[Record], str] = {
    RefTnJunction: "ref_tnjc.csv",
    TnJunction: "tnjc.csv",
    RawTnJc2: "tnjc2_raw.csv",
    CoveredTnJc2: "tnjc2_covered.csv",
    SingleLocusLinkedTnJc2: "tnjc2_classified.csv",
    SynJctsTnJc2: "tnjc2_syn_jcts.csv",
    AnalyzedTnJc2: "tnjc2_analyzed.csv",
    ExportedTnJc2: "tnjc2_exported.csv",
}


def default_filename(record_type: Type[Record]) -> str:
    return DEFAULT_FILENAMES[record_type]


def default_path(run_dir: Path, record_type: Type[Record]) -> Path:
    """Return the default path for a record type within a run directory."""
    return run_dir / default_filename(record_type)

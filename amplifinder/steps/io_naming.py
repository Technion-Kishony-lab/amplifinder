"""Helpers for default filenames / paths for Record-based tables."""

from pathlib import Path
from typing import Dict, Type

from amplifinder.data_types import (
    RefTnJunction, TnJunction, RawTnJc2, CoveredTnJc2, SingleLocusLinkedTnJc2,
    SynJctsTnJc2, AnalyzedTnJc2, ClassifiedTnJc2, ExportedTnJc2,
)
from amplifinder.records.base_records import Record


DEFAULT_FILENAMES: Dict[Type[Record], str] = {
    RefTnJunction: "ref_tnjc.csv",
    TnJunction: "tnjc.csv",
    RawTnJc2: "tnjc2_A_raw.csv",
    CoveredTnJc2: "tnjc2_B_covered.csv",
    SingleLocusLinkedTnJc2: "tnjc2_C_linked.csv",
    SynJctsTnJc2: "tnjc2_D_syn_jcts.csv",
    AnalyzedTnJc2: "tnjc2_E_analyzed.csv",
    ClassifiedTnJc2: "tnjc2_F_classified.csv",
    ExportedTnJc2: "tnjc2_G_final.csv",
}


def default_filename(record_type: Type[Record]) -> str:
    return DEFAULT_FILENAMES[record_type]


def default_path(run_dir: Path, record_type: Type[Record]) -> Path:
    """Return the default path for a record type within a run directory."""
    return run_dir / default_filename(record_type)

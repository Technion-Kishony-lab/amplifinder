"""Helpers for default filenames / paths for Record-based tables."""

from pathlib import Path
from typing import Dict, Type

from amplifinder.data_types import (
    RefTnJunction, TnJunction, TnJc2, CoveredTnJc2, ClassifiedTnJc2,
    CandidateTnJc2, AnalyzedTnJc2,
)
from amplifinder.data_types.records import Record


DEFAULT_FILENAMES: Dict[Type[Record], str] = {
    RefTnJunction: "ref_tn_jc.csv",
    TnJunction: "tn_jc.csv",
    TnJc2: "tn_jc2.csv",
    CoveredTnJc2: "tn_jc2_covered.csv",
    ClassifiedTnJc2: "tn_jc2_classified.csv",
    CandidateTnJc2: "tn_jc2_candidates.csv",
    AnalyzedTnJc2: "tn_jc2_analyzed.csv",
}


def default_filename(record_type: Type[Record]) -> str:
    return DEFAULT_FILENAMES[record_type]


def default_path(output_dir: Path, record_type: Type[Record]) -> Path:
    return output_dir / default_filename(record_type)

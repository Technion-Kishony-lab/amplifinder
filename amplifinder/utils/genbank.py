"""GenBank parsing utilities using BioPython."""

import re
from pathlib import Path
from typing import List, Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature

from amplifinder.data_types.record_types import RefTnLoc


def find_tn_elements(genbank_path: Path, ref_name: str) -> List[RefTnLoc]:
    """Find TN elements from GenBank annotations.

    Looks for mobile_element features with 'insertion sequence' in the type.
    """
    record = SeqIO.read(genbank_path, "genbank")

    records = []
    id_counter = 0

    for feature in record.features:
        tn_name = _extract_tn_name_from_feature(feature)
        if tn_name is None:
            # Not a TN element
            continue

        location = feature.location
        id_counter += 1
        records.append(RefTnLoc(
            tn_id=id_counter,
            tn_name=tn_name,
            tn_scaf=ref_name,
            loc_left=int(location.start) + 1,  # BioPython is 0-based
            loc_right=int(location.end),
            complement=location.strand == -1,
            join=hasattr(location, 'parts') and len(location.parts) > 1,
        ))

    return records


def _extract_tn_name_from_feature(feature: SeqFeature) -> Optional[str]:
    """Extract TN name from a GenBank feature if it's an insertion sequence.

    Returns TN name if feature is a TN element, None otherwise.
    """
    # Get text to search based on feature type
    qualifiers = feature.qualifiers
    if feature.type == "mobile_element":
        text = qualifiers.get("mobile_element_type", [""])[0]
    elif feature.type in ("misc_feature", "repeat_region"):
        text = qualifiers.get("note", [""])[0]
    else:
        return None

    # Check for "insertion sequence" and extract name
    if "insertion sequence" not in text.lower():
        return None

    # Try "insertion sequence:IS1" format
    match = re.search(r'insertion sequence[:\s]*(\S+)', text, re.IGNORECASE)
    if match:
        return match.group(1)

    # Try ISxxx pattern anywhere
    match = re.search(r'(IS\d+\w*)', text, re.IGNORECASE)
    if match:
        return match.group(1)

    return "unknown"

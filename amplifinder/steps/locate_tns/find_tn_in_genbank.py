"""GenBank TN element parsing."""
from __future__ import annotations

import re
from typing import List, Optional

from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord

from amplifinder.data_types.enums import Orientation
from amplifinder.data_types.record_types import RefTn, TnId


def find_tn_elements(gb_records: List[SeqRecord]) -> List[RefTn]:
    """Find TN elements from GenBank annotations.

    Looks for mobile_element features with 'insertion sequence' in the type.
    """
    results = []
    tn_id: TnId = 0

    for record in gb_records:
        scaf_name = record.name
        for feature in record.features:
            tn_name = _extract_tn_name_from_feature(feature)
            if tn_name is None:
                continue

            location = feature.location
            assert location.start >= 0 and location.end >= 0
            assert location.start < location.end
            tn_id += 1
            # BioPython locations are 0-based (start) and 0-based exclusive end
            # Convert to 1-based inclusive for consistency with BLAST coordinates
            results.append(RefTn(
                tn_id=tn_id,
                tn_name=tn_name,
                tn_scaf=scaf_name,
                loc_left=int(location.start) + 1,  # Convert 0-based to 1-based inclusive
                loc_right=int(location.end),       # BioPython end is exclusive, stored as 1-based inclusive
                orientation=Orientation.REVERSE if location.strand == -1 else Orientation.FORWARD,
                join=hasattr(location, 'parts') and len(location.parts) > 1,
            ))

    return results


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

"""GenBank parsing utilities using BioPython."""

import re
from pathlib import Path
from typing import List, Dict, Any, Optional

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature


def find_IS_elements(genbank_path: Path, ref_name: str) -> List[Dict[str, Any]]:
    """Find IS elements from GenBank annotations.
    
    Looks for mobile_element features with 'insertion sequence' in the type.
    
    Args:
        genbank_path: Path to GenBank file
        ref_name: Reference/scaffold name for output
    
    Returns:
        List of IS element dicts with keys: ID, IS_Name, IS_scaf, LocLeft, LocRight, Complement, Join
    """
    record = SeqIO.read(genbank_path, "genbank")
    
    records = []
    id_counter = 0
    
    for feature in record.features:
        IS_name = _extract_IS_name_from_feature(feature)
        if IS_name is None:
            continue
        
        # Get location info
        location = feature.location
        pos_left = int(location.start) + 1  # BioPython is 0-based
        pos_right = int(location.end)
        is_complement = location.strand == -1
        
        # Check for join (compound location)
        is_join = hasattr(location, 'parts') and len(location.parts) > 1
        
        id_counter += 1
        records.append({
            "ID": id_counter,
            "IS_Name": IS_name,
            "IS_scaf": ref_name,
            "LocLeft": pos_left,
            "LocRight": pos_right,
            "Complement": is_complement,
            "Join": is_join,
        })
    
    return records


def _extract_IS_name_from_feature(feature: SeqFeature) -> Optional[str]:
    """Extract IS name from a GenBank feature if it's an insertion sequence.
    
    Returns IS name if feature is an IS element, None otherwise.
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


"""JSON utilities."""

import re


def compact_short_lists(json_str: str) -> str:
    """Compact short lists to single lines."""
    # Pattern to match multi-line lists with short string items
    # Matches: key: [\n  "item1",\n  "item2",\n]
    def compact_list(match):
        prefix = match.group(1)  # Indentation and key
        items_str = match.group(2)  # All items
        items = re.findall(r'"([^"]+)"', items_str)
        # Only compact if <= 5 items and each item is short
        if len(items) <= 5 and all(len(item) < 30 for item in items):
            compact_items = ', '.join(f'"{item}"' for item in items)
            return f'{prefix}[{compact_items}]'
        return match.group(0)
    
    # Match: "key": [\n    "item1",\n    "item2",\n  ]
    pattern = r'(\s+"[^"]+":\s*)\[\s*\n((?:\s+"[^"]+",?\s*\n)*)\s+\]'
    return re.sub(pattern, compact_list, json_str, flags=re.MULTILINE)

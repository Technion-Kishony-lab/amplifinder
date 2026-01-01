# Plan: Refactor classify_structure.py

## 1. Add BaseRawEvent enum
```python
class BaseRawEvent(Enum):
    LJ = "locus-joining"      # Locus-joining (any pair that's not REF or XPOS)
    REF = "reference"          # Both junctions match reference TN at same location
    XPOS = "transposition"     # Short amplicon (< threshold), de novo
```

## 2. Compute base_raw_event for each tnjc2
For each tnjc2:
- **REF**: both jc_num_S == 0 and jc_num_E == 0 and same TN ID
- **XPOS**: not REF and amplicon_length < min_amplicon_length
- **LJ**: everything else (default)

Add as property to ClassifiedTnJc2

## 3. For LJ pairs: find matching single-locus pairs for each side

For each LJ tnjc2, for each side (S/E):
- Find other tnjc2s with base_raw_event = LJ or REF
- Match if they share the same junction ID (jc_num)
- Store as base_raw_event_S and base_raw_event_E (the matching tnjc2 records)

## 4. Assign final RawEvent classification

Based on base_raw_event, base_raw_event_S, base_raw_event_E:
- If base_raw_event = REF → RawEvent.REFERENCE
- If base_raw_event = XPOS → RawEvent.TRANSPOSITION
- If base_raw_event = LJ:
  - Both S and E have matches → FLANKED
  - Only S has match → HEMI_FLANKED_LEFT (account for origin spanning)
  - Only E has match → HEMI_FLANKED_RIGHT (account for origin spanning)
  - Neither has match → UNFLANKED
  - Multiple matches on either side → MULTIPLE_SINGLE_LOCUS

## 5. Update classification logic
Replace current classification algorithm (lines 46-187) with the new approach

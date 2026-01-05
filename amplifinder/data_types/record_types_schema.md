# Record Types Schema Documentation

## Inheritance Hierarchy

### Base Classes
```
Record (from amplifinder.data_types.records)
├── Pydantic BaseModel with schema support
└── Used by all Record subclasses
```

### Enum Types
```
ReversibleIntEnum(int, Enum)
├── Side
│   ├── LEFT = -1
│   └── RIGHT = 1
└── Orientation
    ├── FORWARD = 1
    ├── REVERSE = -1
    └── BOTH = 0

RawEvent(str, Enum)
├── REFERENCE
├── TRANSPOSITION
├── UNFLANKED
├── HEMI_FLANKED_LEFT
├── HEMI_FLANKED_RIGHT
├── FLANKED
├── MULTIPLE_SINGLE_LOCUS
└── UNRESOLVED

JunctionType(int, Enum)
├── LEFT_REF = 1
├── LEFT_IS_TRANS = 2
├── LEFT_MID_IS = 3
├── LOST_IS = 4
├── RIGHT_MID_IS = 5
├── RIGHT_IS_TRANS = 6
└── RIGHT_REF = 7

EventModifier(str, Enum)
├── ANCESTRAL
├── DE_NOVO
└── LOW_COVERAGE
```

### NamedTuple Types
```
Average(NamedTuple)
├── mean: float
├── median: float
└── mode: float

JunctionCoverage(NamedTuple)
├── spanning: int
├── left: int
└── right: int
```

### Record Inheritance Tree

#### TN Element Records
```
Record
├── RefTnSide
│   ├── tn_id: TnId (int)
│   ├── side: Side
│   └── distance: Optional[int]
│
└── SeqRefTnSide(RefTnSide)
    ├── [inherits: tn_id, side, distance]
    ├── offset: int
    └── seq_inward: str
```

#### TN Location Records
```
Record
└── RefTn
    ├── tn_id: TnId (int)
    ├── tn_name: str
    ├── tn_scaf: str
    ├── loc_left: int
    ├── loc_right: int
    ├── orientation: Orientation
    ├── join: bool
    └── [methods: length, get_sides(), get_junctions()]
```

#### Junction Records
```
Record
└── Junction
    ├── num: int
    ├── scaf1: str
    ├── pos1: int
    ├── dir1: Orientation
    ├── scaf2: str
    ├── pos2: int
    ├── dir2: Orientation
    ├── flanking1: int
    ├── flanking2: int
    └── [methods: switch_sides(), get_scaf_pos_dir_flank()]
    │
    ├── RefTnJunction(Junction)
    │   ├── [inherits all Junction fields]
    │   └── ref_tn_side: RefTnSide
    │
    └── TnJunction(Junction)
        ├── [inherits all Junction fields]
        ├── ref_tn_sides: List[RefTnSide]
        └── swapped: bool
```

#### TN Junction Pair (Amplicon) Records
```
Record
└── RawTnJc2
    ├── jc_num_L: int
    ├── jc_num_R: int
    ├── scaf: str
    ├── pos_scaf_L: int
    ├── pos_scaf_R: int
    ├── pos_tn_L: int
    ├── pos_tn_R: int
    ├── dir_scaf_L: Orientation
    ├── dir_scaf_R: Orientation
    ├── dir_tn_L: Orientation
    ├── dir_tn_R: Orientation
    ├── tn_ids: List[int]
    ├── tn_orientations: List[Orientation]
    ├── span_origin: bool
    ├── amplicon_length: int
    └── complementary_length: int
    │
    └── CoveredTnJc2(RawTnJc2)
        ├── [inherits all RawTnJc2 fields]
        ├── iso_amplicon_coverage: Average
        ├── iso_scaf_coverage: Average
        ├── anc_amplicon_coverage: Optional[Average]
        ├── anc_scaf_coverage: Optional[Average]
        ├── copy_number: float
        └── copy_number_vs_anc: Optional[float]
        │
        └── ClassifiedTnJc2(CoveredTnJc2)
            ├── [inherits all CoveredTnJc2 fields]
            ├── raw_event: RawEvent
            ├── shared_tn_ids: List[int]
            └── chosen_tn_id: Optional[int]
            │
            └── FilteredTnJc2(ClassifiedTnJc2)
                ├── [inherits all ClassifiedTnJc2 fields]
                └── analysis_dir: str
                │
                └── AnalyzedTnJc2(FilteredTnJc2)
                    ├── [inherits all FilteredTnJc2 fields]
                    ├── jc_cov_left: List[int]
                    ├── jc_cov_right: List[int]
                    ├── jc_cov_spanning: List[int]
                    ├── anc_jc_cov_left: Optional[List[int]]
                    ├── anc_jc_cov_right: Optional[List[int]]
                    ├── anc_jc_cov_spanning: Optional[List[int]]
                    ├── isolate_architecture: RawEvent
                    ├── ancestor_architecture: Optional[RawEvent]
                    ├── event: str
                    └── event_modifiers: List[EventModifier]
```

#### Other Records
```
Record
├── BlastHit
│   ├── query: str
│   ├── subject: str
│   ├── percent_identical: float
│   ├── length: int
│   ├── mismatch: int
│   ├── gapopen: int
│   ├── qstart: int
│   ├── qend: int
│   ├── sstart: int
│   ├── send: int
│   ├── evalue: float
│   └── bitscore: float
│
└── ExportedTnJc2
    ├── isolate: Optional[str]
    ├── Reference: Optional[str]
    ├── Positions_in_chromosome: Optional[str]
    ├── Direction_in_chromosome: Optional[str]
    ├── amplicon_length: Optional[int]
    ├── IS_element: Optional[str]
    ├── median_copy_number: Optional[float]
    ├── mode_copy_number: Optional[float]
    ├── Ancestor: Optional[str]
    ├── event: Optional[str]
    └── isolate_architecture: Optional[str]
```

## Field Usage Schema

### Type Dependencies (What types use which other types as fields)

#### RefTnSide
- Used in: `RefTnJunction.ref_tn_side`, `TnJunction.ref_tn_sides`, `RefTn.get_sides()`, `RefTn.get_junctions()`

#### SeqRefTnSide
- Extends: `RefTnSide`
- No direct field usage (specialized for sequence matching)

#### RefTn
- No direct field usage (used for generating RefTnSide and RefTnJunction)

#### Junction
- Base class for: `RefTnJunction`, `TnJunction`
- Fields use: `Orientation` (dir1, dir2)

#### RefTnJunction
- Extends: `Junction`
- Fields use: `RefTnSide` (ref_tn_side)
- Created by: `RefTn.get_junctions()`

#### TnJunction
- Extends: `Junction`
- Fields use: `List[RefTnSide]` (ref_tn_sides)

#### RawTnJc2
- Fields use: `Orientation` (dir_scaf_L, dir_scaf_R, dir_tn_L, dir_tn_R), `List[int]` (tn_ids), `List[Orientation]` (tn_orientations)

#### CoveredTnJc2
- Extends: `RawTnJc2`
- Fields use: `Average` (iso_amplicon_coverage, iso_scaf_coverage, anc_amplicon_coverage, anc_scaf_coverage)

#### ClassifiedTnJc2
- Extends: `CoveredTnJc2`
- Fields use: `RawEvent` (raw_event), `List[int]` (shared_tn_ids, chosen_tn_id)

#### FilteredTnJc2
- Extends: `ClassifiedTnJc2`
- No additional type dependencies

#### AnalyzedTnJc2
- Extends: `FilteredTnJc2`
- Fields use: `RawEvent` (isolate_architecture, ancestor_architecture), `List[EventModifier]` (event_modifiers)

#### ExportedTnJc2
- No type dependencies (all fields are primitives/strings)

#### Average
- Used in: `CoveredTnJc2` (iso_amplicon_coverage, iso_scaf_coverage, anc_amplicon_coverage, anc_scaf_coverage)

#### JunctionCoverage
- Not used as field in any Record type (standalone data structure)

#### Side
- Used in: `RefTnSide.side`, `RefTnJunction.ref_tn_side.side`

#### Orientation
- Used in: `Junction.dir1`, `Junction.dir2`, `RawTnJc2.dir_scaf_L`, `RawTnJc2.dir_scaf_R`, `RawTnJc2.dir_tn_L`, `RawTnJc2.dir_tn_R`, `RawTnJc2.tn_orientations`

#### RawEvent
- Used in: `ClassifiedTnJc2.raw_event`, `AnalyzedTnJc2.isolate_architecture`, `AnalyzedTnJc2.ancestor_architecture`

#### JunctionType
- Not used as field (used for indexing junction coverage lists)

#### EventModifier
- Used in: `AnalyzedTnJc2.event_modifiers`

## Pipeline Flow (Record Type Progression)

```
Step 1-6: RawTnJc2
    ↓
Step 7: CoveredTnJc2 (adds coverage)
    ↓
Step 8: ClassifiedTnJc2 (adds raw_event classification)
    ↓
Step 9: FilteredTnJc2 (adds analysis_dir)
    ↓
Step 12: AnalyzedTnJc2 (adds junction coverage analysis)
    ↓
Step 14: ExportedTnJc2 (export format)
```

## Type Aliases

- `TnId = int` - Used in: `RefTnSide.tn_id`, `RefTn.tn_id`, `RawTnJc2.tn_ids`, etc.
- `JunctionT = TypeVar("JunctionT", bound="Junction")` - Used for generic Junction methods


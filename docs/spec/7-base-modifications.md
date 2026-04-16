# Base Modifications (MM/ML Tags)

Base modifications encode epigenetic marks (e.g. 5mC methylation, 6mA, 5hmC) on aligned reads. The SAM 4.5 specification defines two auxiliary tags for this:

- **MM** (Z-type string): modification type and positions, e.g. `C+m,0,2,5;` means "cytosine with 5mC modification at the 1st, 4th, and 10th unmodified C in the query sequence"
- **ML** (B:C array): modification likelihoods (0–255 scale, 255 = highest confidence), one per modification call in MM

These tags are critical for bisulfite/TAPS methylation pipelines and nanopore base modification calling. seqair's TAPS pipeline already writes MM/ML tags during BAM rewriting (see `r[bam.owned_record.aux_array_encoding]`).

> **Sources:** MM/ML tag format is defined in [SAM1] §1.7 (Base modifications). htslib's `hts_base_mod_state` provides the reference parser. The modification type codes (m = 5mC, h = 5hmC, etc.) are defined in the SAM spec's base modification table and registered ChEBI ontology codes.

## Current support: raw pass-through

r[base_mod.passthrough]
MM and ML tags MUST be preserved as raw auxiliary bytes through all I/O paths (BAM read, BAM write, CRAM decode). No parsing or validation of the MM string or ML array contents is performed during pass-through. This ensures that files with base modification annotations are not silently corrupted.

r[base_mod.passthrough.bam]
BAM read → RecordStore → OwnedBamRecord → BAM write MUST preserve MM (Z-type) and ML (B:C array) tags byte-for-byte.

r[base_mod.passthrough.cram]
CRAM decode → RecordStore MUST preserve MM and ML tags. CRAM encodes these via the tag encoding map (`r[cram.compression.tag_encodings]`); the decoder MUST handle arbitrary-length Z and B:C values without truncation.

r[base_mod.passthrough.coexistence]
MM/ML tags MUST coexist with other auxiliary tags (NM, RG, MD, XS, etc.) without interference. Adding, removing, or updating other tags via `AuxData` methods MUST NOT corrupt MM/ML data.

## Planned: structured modification API

The following rules describe the planned structured API for base modifications. They are not yet implemented.

### Design choices

**Type-based API, not free functions.** A `BaseModState` struct holds the parsed MM + ML data for a single record and provides query methods. This mirrors seqair's "types are the primary abstraction" philosophy and parallels `CigarMapping` (pre-extracted CIGAR data) and `PileupAlignment` (pre-extracted alignment fields).

**Lazy parsing with caching.** Modification data is NOT parsed when a record enters the pileup active set. Instead, `BaseModState` is constructed on first access (e.g. when a caller asks for modification data at a column position). The parsed state is cached for subsequent queries on the same record. This avoids overhead when modifications aren't needed (the common case for non-methylation workflows).

**Dual coordinate access.** Modifications can be queried in query-space (matching the stored BAM sequence, zero-cost) or reference-space (mapped through the CIGAR). The pileup integration uses query-space internally since it already has `qpos`.

### BaseModState

r[base_mod.state]
`BaseModState` MUST be constructable from a record's MM string and ML array. It parses all modification entries in the MM string (not filtered by canonical base) and pairs them with the corresponding ML probabilities. The struct is immutable after construction.

r[base_mod.parse_mm]
The MM parser MUST handle the full SAM 4.5 modification string syntax:
- Multiple canonical bases and modification types (e.g. `C+m,0,2;C+h,1,3;A+a,0;`)
- Implicit vs explicit mode: a trailing `.` after the semicolon (e.g. `C+m.,0,2;`) means unlisted positions are definitively unmodified; without `.`, unlisted positions have unknown modification status. The API MUST preserve this distinction.
- ChEBI numeric codes (e.g. `C+27551,0,2;` for 5mC by ChEBI ID)
- Both `+` (top/forward strand) and `-` (bottom/reverse strand) modification strands

r[base_mod.resolve_positions]
`BaseModState` MUST resolve the relative skip-counts in MM into absolute query-sequence positions at construction time. For each modification entry:
1. Iterate through the query sequence (provided at construction)
2. Count only bases matching the canonical base for that entry
3. At each skip-count boundary, record the absolute query position and pair it with the corresponding ML probability

The resolved positions are stored internally as a sorted list for efficient lookup.

### Query API

r[base_mod.query_qpos]
`BaseModState` MUST provide `mod_at_qpos(qpos: usize) -> Option<&[Modification]>` to look up modifications at a given query-sequence position. Returns `None` if no modification is called at that position. Multiple modification types at the same position (e.g. 5mC and 5hmC on the same C) MUST all be returned.

r[base_mod.query_refpos]
`BaseModState` MUST provide `mod_at_ref_pos(ref_pos: i64, cigar: &CigarMapping) -> Option<&[Modification]>` to look up modifications at a reference position. This maps the reference position to a query position via the CIGAR, then delegates to `mod_at_qpos`. Returns `None` if the reference position falls in a deletion/ref-skip (no qpos) or has no modification.

r[base_mod.modification_struct]
Each `Modification` MUST carry:
- `mod_type`: the modification type (single-char code or ChEBI number)
- `probability`: the ML value (u8, 0–255)
- `canonical_base`: which base is modified (A, C, G, T)
- `strand`: `+` or `-`

r[base_mod.implicit_explicit]
`BaseModState` MUST track whether each canonical-base+strand combination uses implicit or explicit mode. A `is_unmodified(qpos: usize, canonical_base: Base) -> Option<bool>` accessor MUST return:
- `Some(true)` if the position is in an explicit-mode entry and not listed (definitively unmodified)
- `Some(false)` if the position has a modification call
- `None` if the position is in an implicit-mode entry and not listed (unknown status)

### Pileup integration

r[base_mod.pileup_integration]
The pileup engine MUST support exposing base modifications via a separate accessor on `PileupColumn`, NOT by extending `PileupOp` (which has a ≤16 byte size constraint). The accessor:

```
column.modification_at(alignment_idx: usize) -> Option<&[Modification]>
```

returns the modification(s) at the current column's reference position for the given alignment, or `None` if the alignment has no MM/ML tags or no modification at this position. Internally, this uses `mod_at_qpos` with the alignment's `qpos`.

`BaseModState` is constructed lazily on first call to `modification_at` for each alignment and cached for subsequent column positions within that alignment's active-set lifetime.

### Reverse complement

r[base_mod.reverse_complement]
For reverse-strand reads (flag 0x10), MM positions are encoded relative to the query sequence **as stored in BAM** (already reverse-complemented). `mod_at_qpos` operates in stored-sequence space and requires no strand adjustment — the MM positions directly index into the stored bases.

`mod_at_ref_pos` MUST also work correctly for reverse-strand reads: the CIGAR-to-qpos mapping already accounts for strand (the CIGAR describes the alignment of the stored sequence to the reference).

Callers who need to reason about the original biological strand (e.g. "this is a methylated C on the bottom strand") can combine the `strand` field from `Modification` with the read's reverse flag.

### Validation

r[base_mod.validation]
When parsing MM tags, the implementation MUST validate:
- The canonical base character is one of A, C, G, T (case-insensitive)
- The strand is `+` or `-`
- Position deltas are non-negative integers
- The total number of modification calls across all MM entries equals the ML array length

Validation failures MUST return typed errors, not panics.

# Indel review

## Summary

The indel implementation is well-designed. The type-safe `PileupOp` enum is a genuine improvement over htslib's flat-field model, and the code is clean, well-tested, and well-specified. However, there is one **blocking correctness issue** (a stale test that now asserts the wrong invariant), one **semantic design concern** worth discussing (orphan insertions), and several smaller items.

## Critical Issues (Blocking)

### 1. `pileup_depth_matches` test is stale and contradicts the new depth semantics

**`compare_pileup.rs:197-216`** — This test compares `rio.depth()` against `hts.alignments.len()`, which only counts non-deletion alignments (because `fetch_htslib_pileup` at line 49 filters with `if let Some(qpos) = a.qpos()`). The comment says:

```rust
// htslib depth includes deletions; alignments.len() only counts bases with qpos.
// seqair depth should match the non-deletion count.
```

But the new spec rule `r[pileup_indel.depth_includes_all]` says depth MUST include deletions and ref-skips. And `total_depth_with_deletions_matches_htslib` (line 275) tests exactly that — seqair depth should match `hts.depth`, not the filtered count.

This test passes today only because `chr19:6_105_700-6_105_800` happens to have no deletions. If you ever add test data with deletions in that region, or change the region, this test will either fail or assert the wrong thing. Right now it's testing a superseded invariant. **Fix: update to compare against `hts.depth`, or delete it since `total_depth_with_deletions_matches_htslib` already covers the correct behavior.**

## Required Changes

### 2. `deletion_ops_match_htslib` cross-validation is weak — it doesn't match by read identity

**`compare_pileup.rs:336-353`** — The cross-validation logic says: "when htslib reports `Indel::Del(len)` at pos P, check that seqair has *any* deletion with matching `del_len` at pos P+1." But it doesn't track *which read* the deletion belongs to. If two different reads have different deletion lengths at adjacent positions, this could pass spuriously. Matching by `record_idx`/`qname` or at minimum by flags would make this airtight.

### 3. `classify_op` insertion detection: `consumes_query` gate silently swallows soft-clips at lookup time

**`cigar.rs:269`** — `classify_op` enters the match/insertion branch for *any* `consumes_query` op type, which includes soft-clips (S). In `pos_info_linear`/`pos_info_bsearch`, soft-clips are skipped because they don't consume ref (`consumes_ref` check at line 298), so they never reach `classify_op`. But if the calling structure ever changes, `classify_op` would produce a `Match` for a soft-clip position. The function's contract is implicit — **add a `debug_assert!` that `op.op_type` is M/=/X when entering the match branch, not just "consumes query."**

### 4. `next_insertion_len` could overflow `u32` with pathological consecutive I ops

**`cigar.rs:249-258`** — `total += next.len` is a `u32` addition with no overflow check. A pathological CIGAR with many large consecutive I ops could silently wrap. Use `total = total.saturating_add(next.len)` or `checked_add` with a debug_assert. In practice impossible with real BAMs, but this is the kind of thing that bites during fuzzing.

## Suggestions

### 5. The "total del_len vs remaining del_len" design choice is correct but deserves a consumer warning

The spec says `del_len` is the same at every position within a deletion — the total D op length, not the remaining bases. This is the right call for variant calling (you want to know "this is a 5bp deletion" not "you're 3bp into a 5bp deletion"). But a consumer who naively iterates columns and sees `del_len=5` five times might think "five separate 5bp deletions" rather than "one 5bp deletion spanning 5 columns." The spec covers this, but `del_len()`'s doc comment could note: "all positions within the same D op report the same value."

### 6. The no-orphan-insertion rule (D-I-M) is defensible but lossy

**`cigar.rs:269-276`** — When a CIGAR has `...M D I M...`, the insertion after the deletion is invisible in the pileup. The spec (`r[pileup_indel.no_orphan_insertions]`) explicitly mandates this because there's "no anchor match position." This is a pragmatic choice — htslib does the same. But a downstream variant caller looking for complex indels (deletion + insertion = replacement) will have to re-parse the CIGAR for these reads. The spec documents this clearly, which is good.

Worth noting: the test `insertion_after_deletion_not_reported` (line 194) validates this, but there's no test for `M I D I M` — an insertion between two D ops that *does* have an anchor. Verify this case is handled correctly (the first I should be reported at the last M, the second I should be orphaned).

### 7. Dedup `to_remove` may contain duplicate entries from symmetric pairs

**`pileup.rs:621`** — `to_remove.dedup()` after `sort_unstable()` handles this, but the only way duplicates arise is if both reads in a pair try to remove each other (which the `seen_record_idx` gate prevents). The `dedup()` call is defensive. Not wrong, but if you're confident the invariant holds, the `dedup` is dead code. If you're not confident, add a `debug_assert!` that no duplicates exist before the dedup.

### 8. `PileupOp` size is 12 bytes — consider layout

`PileupOp::Match` carries `(u32, Base, u8)` = `(4 + 1 + 1)` + padding, `Insertion` adds another `u32` = 12 bytes with discriminant. `Deletion` is just `u32` + discriminant = 8 bytes. The enum is 12 bytes total (sized to largest variant). This is fine, but worth a `static_assert` so nobody accidentally adds a field that balloons it.

### 9. Proptest `arb_cigar_class` doesn't generate `D I M` patterns

**`helpers.rs:213-279`** — The "kitchen sink" pattern generates `H S M I M D M N M S H` but never `M D I M`. This means the no-orphan-insertion rule is only tested via a handwritten test, not via proptest. Add a `D I M` pattern class.

## Verdict

**Needs Discussion** — one stale test asserting a superseded invariant (#1) needs fixing before this is trustworthy. The cross-validation weakness (#2) is a test-quality issue, not a correctness bug. The design decisions (type-safe enum, total del_len, no-orphan insertions) are sound and well-documented.

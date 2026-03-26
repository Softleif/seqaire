//! Tests for performance-related spec rules.
//! after refactoring) rather than measuring wall-clock time.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::{cigar_bytes, cigar_op, make_record, make_record_with_cigar};
use proptest::prelude::*;
use seqair::bam::{RecordStore, cigar::CigarIndex, pileup::PileupEngine};

// ---- perf.reuse_alignment_vec ----
// Verified indirectly: if the engine reuses the vec internally, output must
// still be correct. This test ensures correctness is maintained after the
// optimization.

// r[verify perf.reuse_alignment_vec+2]
proptest! {
    #[test]
    fn reused_alignment_vec_produces_correct_depth(
        records in prop::collection::vec((0i32..200, 20u32..=80), 2..=15),
    ) {
        // Sort by position to match fetch_into's coordinate-sorted guarantee
        let mut sorted_records = records.clone();
        sorted_records.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        let mut intervals: Vec<(i64, i64)> = Vec::new();
        for &(offset, len) in &sorted_records {
            arena.push_raw(&make_record(0, offset, 99, 60, len)).unwrap();
            intervals.push((i64::from(offset), i64::from(offset) + i64::from(len) - 1));
        }

        let engine = PileupEngine::new(arena, 0, 300);
        for col in engine {
            let expected = intervals.iter()
                .filter(|&&(s, e)| col.pos() >= s && col.pos() <= e)
                .count();
            prop_assert_eq!(col.depth(), expected, "depth wrong at pos {}", col.pos());
        }
    }
}

// ---- perf.no_sorted_indices ----
// The engine must produce correct output when records are already sorted
// (which they always are from fetch_into).

// r[verify perf.no_sorted_indices]
#[test]
fn already_sorted_records_produce_correct_pileup() {
    let mut arena = RecordStore::new();
    // Push in sorted order (as fetch_into guarantees)
    arena.push_raw(&make_record(0, 10, 99, 60, 30)).unwrap();
    arena.push_raw(&make_record(0, 20, 99, 60, 30)).unwrap();
    arena.push_raw(&make_record(0, 30, 99, 60, 30)).unwrap();

    let engine = PileupEngine::new(arena, 0, 70);
    let columns: Vec<_> = engine.collect();

    // pos 10-19: 1 read
    assert_eq!(columns.iter().find(|c| c.pos() == 15).unwrap().depth(), 1);
    // pos 20-29: 2 reads
    assert_eq!(columns.iter().find(|c| c.pos() == 25).unwrap().depth(), 2);
    // pos 30-39: 3 reads
    assert_eq!(columns.iter().find(|c| c.pos() == 35).unwrap().depth(), 3);
    // pos 40-49: 2 reads (first ended)
    assert_eq!(columns.iter().find(|c| c.pos() == 45).unwrap().depth(), 2);
}

// ---- perf.avoid_redundant_arena_get ----
// Correctness preserved when arena.get() is called fewer times.

// r[verify perf.avoid_redundant_arena_get+2]
#[test]
fn single_arena_get_per_record_entry_still_correct() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 50)).unwrap();
    arena.push_raw(&make_record(0, 10, 99 | 0x100, 40, 50)).unwrap(); // secondary flag

    let mut engine = PileupEngine::new(arena, 0, 59);
    engine.set_filter(|flags, _aux| flags & 0x100 == 0);
    let columns: Vec<_> = engine.collect();

    // Only the mapq=60 read should pass the filter
    for col in &columns {
        if col.pos() < 10 {
            assert_eq!(col.depth(), 1, "only high-mapq read at pos {}", col.pos());
        }
    }
    // At pos 10+, still only 1 because second read is filtered
    let col = columns.iter().find(|c| c.pos() == 25).unwrap();
    assert_eq!(col.depth(), 1);
}

// ---- perf.cigar_no_to_vec ----
// After removing .to_vec(), CigarIndex must still produce correct qpos.

// r[verify perf.cigar_no_to_vec]
#[test]
fn cigar_index_from_arena_slab_correct() {
    let mut arena = RecordStore::new();
    // 30M 5D 20M — tests that cigar bytes from arena slab work without copying
    let raw = make_record_with_cigar(
        0,
        100,
        99,
        60,
        &[cigar_op(30, 0), cigar_op(5, 2), cigar_op(20, 0)],
        50,
    );
    arena.push_raw(&raw).unwrap();

    let engine = PileupEngine::new(arena, 100, 154);
    let columns: Vec<_> = engine.collect();

    // Before deletion: qpos = pos - 100
    let col = columns.iter().find(|c| c.pos() == 110).unwrap();
    assert_eq!(col.alignments().next().unwrap().qpos(), 10);

    // Inside deletion: no alignment
    assert!(
        columns.iter().find(|c| c.pos() == 132).is_none()
            || columns.iter().find(|c| c.pos() == 132).unwrap().depth() == 0
    );

    // After deletion: qpos = pos - 100 - 5 (deletion consumes 5 ref but 0 query)
    let col = columns.iter().find(|c| c.pos() == 140).unwrap();
    assert_eq!(col.alignments().next().unwrap().qpos(), 35);
}

// ---- perf.precompute_matches_indels ----

// r[verify perf.precompute_matches_indels]
#[test]
fn precomputed_matches_indels_accessible_from_record() {
    let mut arena = RecordStore::new();
    // 30M 5I 15M = 45 matches, 5 indels
    let raw = make_record_with_cigar(
        0,
        0,
        99,
        60,
        &[cigar_op(30, 0), cigar_op(5, 1), cigar_op(15, 0)],
        50,
    );
    arena.push_raw(&raw).unwrap();

    let r = arena.record(0);
    assert_eq!(r.matching_bases, 45);
    assert_eq!(r.indel_bases, 5);
}

// r[verify perf.precompute_matches_indels]
proptest! {
    #[test]
    fn precomputed_matches_indels_consistent_with_calc(
        n_match in 1u32..=200,
        n_ins in 0u32..=20,
        n_del in 0u32..=20,
    ) {
        let mut ops = vec![cigar_op(n_match, 0)];
        let mut seq_len = n_match;
        if n_ins > 0 {
            ops.push(cigar_op(n_ins, 1));
            seq_len += n_ins;
        }
        if n_del > 0 {
            ops.push(cigar_op(n_del, 2));
        }
        // Add trailing match so record is valid
        ops.push(cigar_op(10, 0));
        seq_len += 10;

        let raw = make_record_with_cigar(0, 0, 99, 60, &ops, seq_len);
        let mut arena = RecordStore::new();
        arena.push_raw(&raw).unwrap();
        let r = arena.record(0);

        prop_assert_eq!(r.matching_bases, n_match + 10);
        prop_assert_eq!(r.indel_bases, n_ins + n_del);
    }
}

// ---- perf.cigar_binary_search ----

// r[verify perf.cigar_binary_search]
#[test]
fn binary_search_correct_for_many_ops() {
    use seqair::bam::cigar::CigarIndex;

    // Simulate RNA-seq: 30M 5000N 30M 3000N 40M (5 ops, triggers binary search)
    let ops = cigar_bytes(&[
        cigar_op(30, 0),
        cigar_op(5000, 3),
        cigar_op(30, 0),
        cigar_op(3000, 3),
        cigar_op(40, 0),
    ]);
    let idx = CigarIndex::new(1000, &ops);

    // First M block: 1000-1029
    assert_eq!(idx.qpos_at(1000), Some(0));
    assert_eq!(idx.qpos_at(1029), Some(29));

    // N skip: 1030-6029 → None
    assert_eq!(idx.qpos_at(1030), None);
    assert_eq!(idx.qpos_at(5000), None);

    // Second M block: 6030-6059
    assert_eq!(idx.qpos_at(6030), Some(30));
    assert_eq!(idx.qpos_at(6059), Some(59));

    // Second N skip: 6060-9059 → None
    assert_eq!(idx.qpos_at(6060), None);

    // Third M block: 9060-9099
    assert_eq!(idx.qpos_at(9060), Some(60));
    assert_eq!(idx.qpos_at(9099), Some(99));
    assert_eq!(idx.qpos_at(9100), None);
}

// r[verify perf.cigar_binary_search]
proptest! {
    #[test]
    fn binary_search_matches_linear_for_many_ops(
        n_exons in 3usize..=8,
        exon_len in 20u32..=100,
        intron_len in 100u32..=5000,
    ) {
        // Build alternating M/N CIGAR (RNA-seq pattern)
        let mut ops = Vec::new();
        let mut total_query = 0u32;
        for i in 0..n_exons {
            ops.push(cigar_op(exon_len, 0)); // M
            total_query += exon_len;
            if i < n_exons - 1 {
                ops.push(cigar_op(intron_len, 3)); // N
            }
        }
        let bytes = cigar_bytes(&ops);
        let idx = CigarIndex::new(0, &bytes);

        // Compute total ref span
        let total_ref: i64 = ops.iter().map(|&op| {
            let len = i64::from(op >> 4);
            let op_type = (op & 0xF) as u8;
            if matches!(op_type, 0 | 2 | 3 | 7 | 8) { len } else { 0 }
        }).sum();

        // Check every position — result must be consistent
        let mut last_some: Option<usize> = None;
        for pos in 0..total_ref {
            let result = idx.qpos_at(pos);
            if let Some(qpos) = result {
                prop_assert!(qpos < total_query as usize);
                if let Some(prev) = last_some {
                    prop_assert!(qpos > prev, "qpos not monotonic at ref {pos}");
                }
                last_some = Some(qpos);
            }
        }
        // Out of range
        prop_assert_eq!(idx.qpos_at(total_ref), None);
    }
}

// ---- perf.arena_capacity_hint ----

// r[verify perf.arena_capacity_hint+2]
#[test]
fn arena_with_capacity_avoids_realloc() {
    // Pre-size based on estimated compressed bytes
    let mut arena = RecordStore::with_byte_hint(50_000);

    for i in 0..100 {
        arena.push_raw(&make_record(0, i * 10, 99, 60, 50)).unwrap();
    }
    assert_eq!(arena.len(), 100);

    // All records accessible
    for i in 0..100u32 {
        let r = arena.record(i);
        assert_eq!(r.pos, i64::from(i) * 10);
    }
}

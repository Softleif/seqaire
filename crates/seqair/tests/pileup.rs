//! Tests for PileupEngine: column iteration, filtering, max depth, edge cases.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
mod helpers;

use helpers::{cigar_op, make_record, make_record_with_cigar};
use proptest::prelude::*;
use seqair::bam::{RecordStore, pileup::PileupEngine};
use std::{cell::Cell, rc::Rc};

// ---- pileup.active_set + pileup.column_contents ----

// r[verify pileup.active_set]
// r[verify pileup.column_contents]
proptest! {
    #[test]
    fn depth_equals_covering_records(
        records in prop::collection::vec((0i32..500, 10u32..=100), 1..=20),
    ) {
        let mut sorted = records.clone();
        sorted.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        let mut intervals: Vec<(i64, i64)> = Vec::new();
        for &(offset, read_len) in &sorted {
            let pos = i64::from(offset);
            let end = pos + i64::from(read_len) - 1;
            arena.push_raw(&make_record(0, offset, 99, 60, read_len)).unwrap();
            intervals.push((pos, end));
        }

        let engine = PileupEngine::new(arena, 0, 600);
        for col in engine {
            let expected = intervals.iter()
                .filter(|&&(s, e)| col.pos() >= s && col.pos() <= e)
                .count();
            prop_assert_eq!(col.depth(), expected,
                "depth mismatch at pos {}", col.pos());
        }
    }
}

// ---- pileup.position_iteration ----

// r[verify pileup.position_iteration]
proptest! {
    #[test]
    fn every_covered_position_yields_a_column(
        records in prop::collection::vec((0i32..500, 10u32..=100), 1..=20),
    ) {
        let mut sorted = records.clone();
        sorted.sort_by_key(|&(offset, _)| offset);

        let mut arena = RecordStore::new();
        let mut intervals: Vec<(i64, i64)> = Vec::new();
        for &(offset, read_len) in &sorted {
            arena.push_raw(&make_record(0, offset, 99, 60, read_len)).unwrap();
            intervals.push((i64::from(offset), i64::from(offset) + i64::from(read_len) - 1));
        }

        let engine = PileupEngine::new(arena, 0, 600);
        let columns: Vec<_> = engine.collect();
        let col_positions: std::collections::HashSet<i64> = columns.iter().map(|c| c.pos()).collect();

        for pos in 0..=600i64 {
            let covered = intervals.iter().any(|&(s, e)| pos >= s && pos <= e);
            prop_assert_eq!(col_positions.contains(&pos), covered,
                "position {} coverage mismatch", pos);
        }
    }
}

// ---- pileup.read_filter ----

// r[verify pileup.read_filter]
#[test]
fn filter_evaluated_once_per_record() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 0, 99, 60, 20)).unwrap();

    let count = Rc::new(Cell::new(0usize));
    let count_clone = count.clone();
    let mut engine = PileupEngine::new(arena, 0, 19);
    engine.set_filter(move |_flags, _aux| {
        count_clone.set(count_clone.get() + 1);
        true
    });

    let columns: Vec<_> = engine.collect();
    assert_eq!(columns.len(), 20);
    assert_eq!(count.get(), 1, "filter called once per record, not per column");
}

// r[verify pileup.read_filter]
proptest! {
    #[test]
    fn filter_by_flags_excludes_correct_reads(
        pass_flags in prop::collection::vec(prop::bool::ANY, 1..=20),
    ) {
        let mut arena = RecordStore::new();
        for &pass in &pass_flags {
            let flags = if pass { 99 } else { 99 | 0x100 }; // 0x100 = secondary
            arena.push_raw(&make_record(0, 0, flags, 60, 10)).unwrap();
        }

        let expected = pass_flags.iter().filter(|&&p| p).count();
        let mut engine = PileupEngine::new(arena, 0, 0);
        engine.set_filter(move |flags, _aux| flags & 0x100 == 0);
        let columns: Vec<_> = engine.collect();

        if expected > 0 {
            prop_assert_eq!(columns[0].depth(), expected);
        } else {
            prop_assert!(columns.is_empty());
        }
    }
}

// ---- pileup.max_depth + pileup.max_depth_per_position ----

// r[verify pileup.max_depth]
proptest! {
    #[test]
    fn max_depth_never_exceeded(max in 1u32..=10, n in 1usize..=30, len in 10u32..=50) {
        let mut arena = RecordStore::new();
        for _ in 0..n {
            arena.push_raw(&make_record(0, 0, 99, 60, len)).unwrap();
        }

        let mut engine = PileupEngine::new(arena, 0, i64::from(len) - 1);
        engine.set_max_depth(max);
        for col in engine {
            prop_assert!(col.depth() <= max as usize);
        }
    }
}

// r[verify pileup.max_depth_per_position]
proptest! {
    #[test]
    fn reads_available_at_lower_coverage_despite_cap(extra in 3usize..=10) {
        let mut arena = RecordStore::new();
        for _ in 0..2 {
            arena.push_raw(&make_record(0, 0, 99, 60, 100)).unwrap();
        }
        for _ in 0..extra {
            arena.push_raw(&make_record(0, 0, 99, 60, 20)).unwrap();
        }

        let mut engine = PileupEngine::new(arena, 0, 99);
        engine.set_max_depth(3);
        let columns: Vec<_> = engine.collect();

        let col0 = columns.iter().find(|c| c.pos() == 0).unwrap();
        prop_assert_eq!(col0.depth(), 3); // capped

        let col50 = columns.iter().find(|c| c.pos() == 50).unwrap();
        prop_assert_eq!(col50.depth(), 2); // below cap
    }
}

// ---- pileup.qpos ----

// r[verify pileup.qpos]
proptest! {
    #[test]
    fn qpos_for_simple_cigar_is_offset_from_start(start in 0i32..100, len in 10u32..=50) {
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, start, 99, 60, len)).unwrap();

        let start_i64 = i64::from(start);
        let end_i64 = start_i64 + i64::from(len) - 1;
        let engine = PileupEngine::new(arena, start_i64, end_i64);

        for col in engine {
            for aln in col.alignments() {
                prop_assert_eq!(aln.qpos(), (col.pos() - start_i64) as usize);
            }
        }
    }
}

// ---- pileup.empty_positions_skipped ----

// r[verify pileup.empty_positions_skipped]
proptest! {
    #[test]
    fn no_empty_columns(pos1 in 0i32..100, len1 in 5u32..=20, gap in 50i32..=500, len2 in 5u32..=20) {
        let pos2 = pos1 + i32::from(len1 as u16) + gap;
        let mut arena = RecordStore::new();
        arena.push_raw(&make_record(0, pos1, 99, 60, len1)).unwrap();
        arena.push_raw(&make_record(0, pos2, 99, 60, len2)).unwrap();

        let end = i64::from(pos2) + i64::from(len2) + 10;
        let engine = PileupEngine::new(arena, i64::from(pos1), end);
        for col in engine {
            prop_assert!(col.depth() > 0, "empty column at {}", col.pos());
        }
    }
}

// ---- pileup.unmapped_excluded ----

// r[verify pileup.unmapped_excluded]
#[test]
fn unmapped_reads_excluded_from_pileup() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 100, 99, 60, 50)).unwrap();
    arena.push_raw(&make_record(0, 100, 0x4, 60, 50)).unwrap(); // unmapped
    arena.push_raw(&make_record(0, 100, 163, 60, 50)).unwrap();

    let engine = PileupEngine::new(arena, 100, 149);
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns[0].depth(), 2, "unmapped should be excluded");
}

// ---- pileup.zero_refspan_reads ----

// r[verify pileup.zero_refspan_reads]
#[test]
fn zero_refspan_read_handled_gracefully() {
    let raw = make_record_with_cigar(0, 50, 99, 60, &[cigar_op(10, 4)], 10); // 10S
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();
    assert_eq!(arena.record(0).end_pos, 50);

    let engine = PileupEngine::new(arena, 48, 52);
    let columns: Vec<_> = engine.collect();
    // Pure soft-clip has no ref-consuming ops → no qpos → no column
    assert!(columns.is_empty() || columns.iter().all(|c| c.pos() != 50 || c.depth() == 0));
}

// ---- pileup.soft_clip_at_position ----

// r[verify pileup.soft_clip_at_position]
#[test]
fn leading_softclip_does_not_extend_ref_range() {
    // 5S 20M at pos 100 → ref 100-119, query 5-24
    let raw = make_record_with_cigar(0, 100, 99, 60, &[cigar_op(5, 4), cigar_op(20, 0)], 25);
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();

    let engine = PileupEngine::new(arena, 95, 125);
    let columns: Vec<_> = engine.collect();
    let positions: Vec<i64> = columns.iter().map(|c| c.pos()).collect();

    assert_eq!(positions.first(), Some(&100), "should start at first ref-consuming pos");
    assert_eq!(positions.last(), Some(&119));
    assert_eq!(columns[0].alignments().next().unwrap().qpos(), 5);
}

// ---- bam.reader edge cases ----

// r[verify bam.reader.secondary_supplementary_included+2]
#[test]
fn secondary_and_supplementary_reads_in_arena() {
    let mut arena = RecordStore::new();
    arena.push_raw(&make_record(0, 100, 99, 60, 50)).unwrap();
    arena.push_raw(&make_record(0, 100, 0x100 | 99, 60, 50)).unwrap(); // secondary
    arena.push_raw(&make_record(0, 100, 0x800 | 99, 60, 50)).unwrap(); // supplementary

    assert_eq!(arena.len(), 3);

    let engine = PileupEngine::new(arena, 100, 149);
    let columns: Vec<_> = engine.collect();
    assert_eq!(columns[0].depth(), 3, "secondary+supplementary should appear by default");
}

// ---- bam.record.zero_refspan ----

// r[verify bam.record.zero_refspan]
#[test]
fn zero_refspan_endpos_equals_pos() {
    let raw = make_record_with_cigar(0, 100, 99, 60, &[cigar_op(10, 4)], 10); // 10S
    let mut arena = RecordStore::new();
    arena.push_raw(&raw).unwrap();
    let r = arena.record(0);
    assert_eq!(r.end_pos, 100, "zero-refspan should have end_pos == pos");
}

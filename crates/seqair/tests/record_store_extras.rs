//! Tests for `RecordStore<U>` per-record extras and `PileupEngine` columns_with_store.
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    reason = "test code"
)]
#![allow(
    clippy::arithmetic_side_effects,
    clippy::cast_possible_truncation,
    reason = "test code with known small values"
)]

mod helpers;

use helpers::make_record;
use seqair::bam::Pos0;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::RecordStore;

/// Helper: push N synthetic records at positions 100, 101, ... with 10M CIGAR.
fn store_with_n_records(n: u32) -> RecordStore {
    let mut store = RecordStore::new();
    for i in 0..n {
        let pos = 100 + i as i32;
        store.push_raw(&make_record(0, pos, 99, 60, 10)).unwrap();
    }
    store
}

// ---- RecordStore extras ----

// r[verify record_store.extras.with_extras]
// r[verify record_store.extras.access]
#[test]
fn with_extras_computes_per_record_data() {
    let store = store_with_n_records(5);

    // Compute extras: store each record's position as the extra.
    let store = store.with_extras(|idx, s| s.record(idx).pos.as_i32());

    assert_eq!(store.len(), 5);
    for i in 0..5u32 {
        assert_eq!(*store.extra(i), 100 + i as i32);
    }
}

// r[verify record_store.extras.with_extras]
#[test]
fn with_extras_preserves_slab_data() {
    let store = store_with_n_records(3);

    // Capture original data before transformation.
    let original_positions: Vec<i32> = (0..3u32).map(|i| store.record(i).pos.as_i32()).collect();
    let original_mapqs: Vec<u8> = (0..3u32).map(|i| store.record(i).mapq).collect();

    let store = store.with_extras(|idx, s| s.record(idx).mapq as u16 * 10);

    // Slab data must be moved, not copied — verify it's still correct.
    for i in 0..3u32 {
        assert_eq!(store.record(i).pos.as_i32(), original_positions[i as usize]);
        assert_eq!(store.record(i).mapq, original_mapqs[i as usize]);
    }
    // Extras computed correctly.
    assert_eq!(*store.extra(0), 600);
    assert_eq!(*store.extra(1), 600);
    assert_eq!(*store.extra(2), 600);
}

// r[verify record_store.extras.access]
#[test]
fn extra_mut_allows_modification() {
    let store = store_with_n_records(3);
    let mut store = store.with_extras(|_, _| 0u32);

    *store.extra_mut(1) = 42;

    assert_eq!(*store.extra(0), 0);
    assert_eq!(*store.extra(1), 42);
    assert_eq!(*store.extra(2), 0);
}

// r[verify record_store.extras.with_extras]
#[test]
fn with_extras_closure_can_read_aux() {
    // Build a record with aux data and verify the closure can access it.
    let aux = b"NMCd"; // NM:C:100 (tag NM, type C=u8, value 100)
    let raw = helpers_make_test_record_with_aux(0, 100, 99, 60, 10, aux);

    let mut store = RecordStore::new();
    store.push_raw(&raw).unwrap();

    let store = store.with_extras(|idx, s| {
        let aux_bytes = s.aux(idx);
        aux_bytes.len()
    });

    assert_eq!(*store.extra(0), aux.len());
}

// r[verify record_store.extras.clear]
#[test]
fn clear_clears_extras_and_retains_capacity() {
    let store = store_with_n_records(50);
    let mut store = store.with_extras(|idx, s| s.record(idx).pos.as_i32());

    assert_eq!(store.len(), 50);
    let records_cap = store.records_capacity();
    let extras_cap = store.extras_capacity();
    assert!(extras_cap >= 50);

    store.clear();

    assert_eq!(store.len(), 0);
    assert!(store.records_capacity() >= records_cap);
    assert!(store.extras_capacity() >= extras_cap);
}

// r[verify record_store.extras.strip]
#[test]
fn strip_extras_preserves_slab_capacity() {
    let store = store_with_n_records(50);
    let store = store.with_extras(|idx, s| s.record(idx).pos.as_i32());

    let records_cap = store.records_capacity();
    let names_cap = store.names_capacity();
    let qual_cap = store.qual_capacity();

    let store = store.strip_extras();

    // Slab capacities preserved.
    assert!(store.records_capacity() >= records_cap);
    assert!(store.names_capacity() >= names_cap);
    assert!(store.qual_capacity() >= qual_cap);

    // Records still accessible (they were moved, not cleared).
    assert_eq!(store.len(), 50);
}

// r[verify record_store.extras.generic_param]
#[test]
fn unit_extras_is_zero_cost() {
    // Vec<()> elements are ZST — no heap allocation for element storage.
    assert_eq!(std::mem::size_of::<()>(), 0);

    // A RecordStore<()> (the default) should work identically to the old RecordStore.
    let store = store_with_n_records(3);
    assert_eq!(store.len(), 3);
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
}

// r[verify record_store.extras.sort_dedup_unit_only]
#[test]
fn sort_and_dedup_work_on_unit_store() {
    // These methods are only on RecordStore<()> — this is a compile-time check.
    let mut store = RecordStore::new();
    // Push records out of order.
    store.push_raw(&make_record(0, 200, 99, 60, 10)).unwrap();
    store.push_raw(&make_record(0, 100, 99, 60, 10)).unwrap();

    store.sort_by_pos();
    assert_eq!(store.record(0).pos, Pos0::new(100).unwrap());
    assert_eq!(store.record(1).pos, Pos0::new(200).unwrap());
}

// ---- PileupEngine extras ----

// r[verify pileup.extras.with_extras]
#[test]
fn engine_with_extras_preserves_settings() {
    let store = store_with_n_records(5);
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    engine.set_max_depth(10);
    engine.set_filter(|flags, _aux| !flags.is_duplicate());

    let engine = engine.with_extras(|idx, s| s.record(idx).mapq);

    // Verify the engine still works — iterate and check columns are produced.
    let columns: Vec<_> = engine.collect();
    assert!(!columns.is_empty());
}

// r[verify pileup.extras.columns_with_store]
#[test]
fn columns_with_store_yields_same_columns_as_iterator() {
    let store = store_with_n_records(5);

    // Collect columns via Iterator.
    let store_copy = store_with_n_records(5);
    let engine = PileupEngine::new(store_copy, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    let iter_columns: Vec<_> = engine.map(|c| (c.pos(), c.depth())).collect();

    // Collect columns via columns_with_store.
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());
    let mut cws_columns = Vec::new();
    let mut cols = engine.columns_with_store();
    while let Some((col, _store)) = cols.next_column() {
        cws_columns.push((col.pos(), col.depth()));
    }

    assert_eq!(iter_columns, cws_columns);
}

// r[verify pileup.extras.columns_with_store]
#[test]
fn columns_with_store_provides_store_access() {
    let store = store_with_n_records(3);
    let store = store.with_extras(|idx, s| s.record(idx).pos.as_i32());
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());

    let mut saw_extras = false;
    let mut cols = engine.columns_with_store();
    while let Some((col, store)) = cols.next_column() {
        for aln in col.alignments() {
            let extra = store.extra(aln.record_idx());
            // Extra should be the record's position.
            assert!(*extra >= 100 && *extra <= 102);
            saw_extras = true;
        }
    }
    assert!(saw_extras, "should have seen at least one alignment with extras");
}

// r[verify pileup.extras.recover_store]
#[test]
fn recover_store_works_with_extras_engine() {
    // Build a store manually, transform to extras, feed into engine, then strip.
    let store = store_with_n_records(10);
    let store = store.with_extras(|idx, s| s.record(idx).mapq as u32);
    let mut engine = PileupEngine::new(store, Pos0::new(100).unwrap(), Pos0::new(109).unwrap());

    // Consume all columns.
    let mut cols = engine.columns_with_store();
    while let Some((_col, _store)) = cols.next_column() {}
    drop(cols);

    // take_store + strip_extras (the same operation recover_store does internally).
    let store = engine.take_store().expect("store should be available");
    let store: RecordStore = store.strip_extras();

    // The recovered store should be empty (cleared by take_store) but have capacity.
    assert_eq!(store.len(), 0);
    assert!(store.records_capacity() > 0);
}

// ---- helpers ----

/// Build a BAM record with aux data for testing aux access in with_extras.
fn helpers_make_test_record_with_aux(
    tid: i32,
    pos: i32,
    flags: u16,
    mapq: u8,
    seq_len: u32,
    aux: &[u8],
) -> Vec<u8> {
    let name = b"read\0";
    let name_len = name.len();
    let cigar_op = seq_len << 4; // N x M
    let seq_bytes = (seq_len as usize).div_ceil(2);

    let total = 32 + name_len + 4 + seq_bytes + seq_len as usize + aux.len();
    let mut raw = vec![0u8; total];

    raw[0..4].copy_from_slice(&tid.to_le_bytes());
    raw[4..8].copy_from_slice(&pos.to_le_bytes());
    raw[8] = name_len as u8;
    raw[9] = mapq;
    raw[12..14].copy_from_slice(&1u16.to_le_bytes()); // 1 cigar op
    raw[14..16].copy_from_slice(&flags.to_le_bytes());
    raw[16..20].copy_from_slice(&seq_len.to_le_bytes());
    raw[20..32].fill(0);

    raw[32..32 + name_len].copy_from_slice(name);
    let cigar_start = 32 + name_len;
    raw[cigar_start..cigar_start + 4].copy_from_slice(&cigar_op.to_le_bytes());

    let seq_start = cigar_start + 4;
    for i in 0..seq_bytes {
        raw[seq_start + i] = 0x11; // A,A
    }

    let qual_start = seq_start + seq_bytes;
    for i in 0..seq_len as usize {
        raw[qual_start + i] = 30;
    }

    let aux_start = qual_start + seq_len as usize;
    raw[aux_start..aux_start + aux.len()].copy_from_slice(aux);

    raw
}

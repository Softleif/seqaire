//! Full-stack IndexedBamReader fuzz: BAM bytes → tmpfile → open → fetch_into → pileup.
//!
//! Uses tmpfiles because RegionBuf requires seekable File handles.
//! Seeded with real test.bam + test.bam.bai via symlinks for high coverage.
//! Input: raw bytes split into BAM content + BAI content.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::reader::IndexedBamReader;
use seqair::bam::record_store::RecordStore;
use seqair_types::{Offset, Pos, Zero};
use std::io::Write;

fuzz_target!(|data: &[u8]| {
    // Need at least 8 bytes: 4 for split point, 4 minimum content
    if data.len() < 8 || data.len() > 256 * 1024 {
        return;
    }

    // First 4 bytes encode the split point between BAM and BAI data
    let split = u32::from_le_bytes([data[0], data[1], data[2], data[3]]) as usize;
    let split = split % (data.len() - 4); // Ensure split is within bounds
    let bam_data = &data[4..4 + split];
    let bai_data = &data[4 + split..];

    // Need minimum sizes for both
    if bam_data.len() < 4 || bai_data.len() < 8 {
        return;
    }

    let dir = match tempfile::tempdir() {
        Ok(d) => d,
        Err(_) => return,
    };

    let bam_path = dir.path().join("test.bam");
    let bai_path = dir.path().join("test.bam.bai");

    if std::fs::File::create(&bam_path).and_then(|mut f| f.write_all(bam_data)).is_err() {
        return;
    }
    if std::fs::File::create(&bai_path).and_then(|mut f| f.write_all(bai_data)).is_err() {
        return;
    }

    // Open with the real reader (format detection, header parsing, index loading)
    let mut reader = match IndexedBamReader::open(&bam_path) {
        Ok(r) => r,
        Err(_) => return,
    };

    let header = reader.header();
    let target_count = header.target_count();

    // Try fetching from first contig
    if target_count == 0 {
        return;
    }

    let start = Pos::<Zero>::new(0).unwrap();
    let end = match start.checked_add_offset(Offset::new(10_000)) {
        Some(p) => p,
        None => return,
    };

    let mut store = RecordStore::new();
    let _ = reader.fetch_into(0, start, end, &mut store);

    if store.len() == 0 {
        return;
    }

    // Run pileup over fetched records
    let mut engine = PileupEngine::new(store, start, end);
    engine.set_max_depth(50);
    for col in engine.by_ref().take(1000) {
        let _depth = col.depth();
        for aln in col.alignments() {
            let _op = aln.op();
            let _base = aln.base();
        }
    }
});

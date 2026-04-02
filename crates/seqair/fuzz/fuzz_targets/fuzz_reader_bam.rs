#![no_main]

use arbitrary::Arbitrary;
use libfuzzer_sys::fuzz_target;
use seqair::bam::bgzf::BgzfReader;
use seqair::bam::header::BamHeader;
use seqair::bam::index::BamIndex;
use seqair::bam::pileup::PileupEngine;
use seqair::bam::record_store::RecordStore;
use seqair_types::{Offset, Pos, Zero};

/// Full-stack BAM reader fuzz: parse BGZF stream → BAM header → records → pileup.
/// Uses cursor-backed BgzfReader (no file I/O).
#[derive(Arbitrary, Debug)]
struct BamReaderInput {
    /// Raw BAM file content (BGZF-compressed stream)
    bam_data: Vec<u8>,
    /// Raw BAI index content
    bai_data: Vec<u8>,
    /// Region to query
    region_start: u16,
    region_len: u16,
}

fuzz_target!(|input: BamReaderInput| {
    // Limit input sizes
    if input.bam_data.len() > 128 * 1024 || input.bai_data.len() > 16 * 1024 {
        return;
    }

    // Parse BAI index
    let _index = match BamIndex::from_bytes(&input.bai_data) {
        Ok(idx) => idx,
        Err(_) => return,
    };

    // Open BGZF cursor and parse BAM header
    let mut bgzf = BgzfReader::from_cursor(input.bam_data);
    let header = match BamHeader::parse(&mut bgzf) {
        Ok(h) => h,
        Err(_) => return,
    };

    let _target_count = header.target_count();

    // Read records through the BGZF layer into RecordStore
    let mut store = RecordStore::new();
    let max_records = 64;
    let mut buf = [0u8; 4];
    for _ in 0..max_records {
        // Read 4-byte block_size
        if bgzf.read_exact_into(&mut buf).is_err() {
            break;
        }
        let block_size = i32::from_le_bytes(buf) as usize;
        if block_size < 32 || block_size > 100_000 {
            break;
        }

        // Read record body
        let mut record_buf = vec![0u8; block_size];
        if bgzf.read_exact_into(&mut record_buf).is_err() {
            break;
        }

        let _ = store.push_raw(&record_buf);
    }

    if store.len() == 0 {
        return;
    }

    // Run pileup over a small region
    let region_start = match Pos::<Zero>::new(input.region_start as u32) {
        Some(p) => p,
        None => return,
    };
    let region_len = (input.region_len as u32).min(500);
    let region_end = match region_start.checked_add_offset(Offset::new(region_len as i64)) {
        Some(p) => p,
        None => return,
    };

    let mut engine = PileupEngine::new(store, region_start, region_end);
    engine.set_max_depth(50);
    for col in engine.by_ref().take(500) {
        let _depth = col.depth();
        for aln in col.alignments() {
            let _op = aln.op();
            let _base = aln.base();
        }
    }
});

//! Full-stack fuzz through IndexedReader<Cursor>: exercises the unified reader
//! enum with cursor-backed I/O for BAM and CRAM paths.
#![no_main]

use libfuzzer_sys::fuzz_target;
use seqair::{
    bam::{pileup::PileupEngine, record_store::RecordStore},
    reader::FuzzReaders,
};
use seqair_fuzz::indexed_reader::{FastaInput, Input, ReadersInput};
use seqair_types::{Offset, Pos, Zero};

fuzz_target!(|data: Input| {
    run(data);
});

fn run(data: Input) {
    let FastaInput { fasta_gz, fai, gzi } = data.alignment;
    let reader = match data.readers {
        ReadersInput::Bam { bam, bai } => {
            FuzzReaders::from_bam_bytes(bam, &bai, fasta_gz, &fai, &gzi)
        }
        ReadersInput::Sam { sam: _, sai: _ } => return,
        ReadersInput::Cram { cram, crai } => {
            FuzzReaders::from_cram_bytes(cram, &crai, fasta_gz, &fai, &gzi)
        }
    };
    let Ok(mut reader) = reader else {
        return;
    };

    let header = reader.header();
    let target_count = header.target_count();
    if target_count == 0 {
        return;
    }
    for i in 0..target_count.min(50) {
        let _ = header.target_name(i as u32);
        let _ = header.target_len(i as u32);
    }

    let start = Pos::<Zero>::new(0).unwrap();
    let Some(end) = start.checked_add_offset(Offset::new(10_000)) else {
        return;
    };

    let mut store = RecordStore::new();
    let _ = reader.fetch_into(0, start, end, &mut store);

    if store.len() == 0 {
        return;
    }

    let mut engine = PileupEngine::new(store, start, end);
    engine.set_max_depth(50);
    for col in engine.by_ref().take(1000) {
        let _depth = col.depth();
        let _mdepth = col.match_depth();
        for aln in col.alignments() {
            let _op = aln.op();
            let _base = aln.base();
            let _qual = aln.qual();
        }
    }
}

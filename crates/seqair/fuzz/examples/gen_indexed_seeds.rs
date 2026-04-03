//! Generate seeds for `fuzz_reader_indexed` from real test data.
//!
//! Uses `Input::encode` for exact round-trip fidelity — every byte of the seed
//! is used by the fuzz target, no waste.
//!
//! Run: cd crates/seqair/fuzz && cargo run --example gen_indexed_seeds

use seqair_fuzz::indexed_reader::Input;
use std::fs;
use std::io::Read;
use std::path::Path;

fn main() {
    let data_dir = Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/data");
    let seed_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("seeds/fuzz_reader_indexed");
    fs::create_dir_all(&seed_dir).unwrap();

    let fai = fs::read(data_dir.join("test.fasta.gz.fai")).unwrap();
    let gzi = fs::read(data_dir.join("test.fasta.gz.gzi")).unwrap();

    // === BAM seed ===
    let bam = fs::read(data_dir.join("test.bam")).unwrap();
    let bai = fs::read(data_dir.join("test.bam.bai")).unwrap();
    let small_bam = truncate_bgzf(&bam, 5);

    let bam_seed = Input::encode(0, &small_bam, &bai, &[], &fai, &gzi);
    // Verify round-trip
    let parsed = Input::parse(&bam_seed).expect("BAM seed must round-trip");
    assert_eq!(parsed.data1.len(), small_bam.len());
    assert_eq!(parsed.data2.len(), bai.len());
    fs::write(seed_dir.join("bam_seed"), &bam_seed).unwrap();
    println!(
        "bam_seed: {} bytes (bam={}, bai={}, fai={}, gzi={})",
        bam_seed.len(),
        small_bam.len(),
        bai.len(),
        fai.len(),
        gzi.len(),
    );

    // === CRAM seed ===
    let cram = fs::read(data_dir.join("test.cram")).unwrap();
    let crai_gz = fs::read(data_dir.join("test.cram.crai")).unwrap();
    let crai_text = decompress_gz(&crai_gz);
    let small_cram = &cram[..cram.len().min(50_000)];

    let cram_seed = Input::encode(2, small_cram, &crai_text, &[], &fai, &gzi);
    let parsed = Input::parse(&cram_seed).expect("CRAM seed must round-trip");
    assert_eq!(parsed.data1.len(), small_cram.len());
    assert_eq!(parsed.data2.len(), crai_text.len());
    fs::write(seed_dir.join("cram_seed"), &cram_seed).unwrap();
    println!(
        "cram_seed: {} bytes (cram={}, crai={}, fai={}, gzi={})",
        cram_seed.len(),
        small_cram.len(),
        crai_text.len(),
        fai.len(),
        gzi.len(),
    );
}

fn truncate_bgzf(data: &[u8], max_blocks: usize) -> Vec<u8> {
    let mut pos = 0;
    let mut blocks = 0;
    while pos < data.len() && blocks < max_blocks {
        if data.get(pos..pos + 2) != Some(&[0x1f, 0x8b]) {
            break;
        }
        let xlen = u16::from_le_bytes([data[pos + 10], data[pos + 11]]) as usize;
        let extra = &data[pos + 12..pos + 12 + xlen];
        let mut ep = 0;
        let mut bsize = None;
        while ep + 4 <= extra.len() {
            let (si1, si2) = (extra[ep], extra[ep + 1]);
            let slen = u16::from_le_bytes([extra[ep + 2], extra[ep + 3]]) as usize;
            if si1 == b'B' && si2 == b'C' && slen == 2 {
                bsize = Some(u16::from_le_bytes([extra[ep + 4], extra[ep + 5]]));
                break;
            }
            ep += 4 + slen;
        }
        let Some(bsize) = bsize else { break };
        pos += bsize as usize + 1;
        blocks += 1;
    }
    data[..pos].to_vec()
}

fn decompress_gz(data: &[u8]) -> Vec<u8> {
    let mut decoder = flate2::read::GzDecoder::new(data);
    let mut out = Vec::new();
    decoder.read_to_end(&mut out).unwrap_or(0);
    out
}

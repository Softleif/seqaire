//! Generate seeds for `fuzz_reader_indexed` from real test data.
//!
//! Concatenates test files as raw bytes — Arbitrary will split them however it
//! can, and the fuzzer will mutate to find valid splits. Also verifies that the
//! seed can be parsed by Arbitrary (even if the extracted sizes aren't perfect,
//! the fuzzer will fix that through mutation).
//!
//! Run: cd crates/seqair/fuzz && cargo run --example gen_indexed_seeds

use arbitrary::Unstructured;
use seqair_fuzz::indexed_reader::{Input, ReadersInput};
use std::fs;
use std::io::Read;
use std::path::Path;

fn main() {
    // CARGO_MANIFEST_DIR = crates/seqair/fuzz
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

    // Shared FASTA data
    let fai = fs::read(data_dir.join("test.fasta.gz.fai")).unwrap();
    let gzi = fs::read(data_dir.join("test.fasta.gz.gzi")).unwrap();

    // === BAM seed ===
    let bam = fs::read(data_dir.join("test.bam")).unwrap();
    let bai = fs::read(data_dir.join("test.bam.bai")).unwrap();
    let small_bam = truncate_bgzf(&bam, 5);

    // Variant 0 = Bam. Place discriminant byte at end so Arbitrary picks it up.
    let bam_seed = build_seed(0, &small_bam, &bai, &fai, &gzi);
    report_parse(&bam_seed, "BAM");
    fs::write(seed_dir.join("bam_seed"), &bam_seed).unwrap();
    println!("  bam_seed: {} bytes written", bam_seed.len());

    // === CRAM seed ===
    let cram = fs::read(data_dir.join("test.cram")).unwrap();
    let crai_gz = fs::read(data_dir.join("test.cram.crai")).unwrap();
    let crai_text = decompress_gz(&crai_gz);
    let small_cram = &cram[..cram.len().min(50_000)];

    let cram_seed = build_seed(2, small_cram, &crai_text, &fai, &gzi);
    report_parse(&cram_seed, "CRAM");
    fs::write(seed_dir.join("cram_seed"), &cram_seed).unwrap();
    println!("  cram_seed: {} bytes written", cram_seed.len());
}

/// Build a seed by concatenating data with the variant discriminant at the end.
///
/// Arbitrary reads lengths from the end of the buffer. By putting all the valid
/// file data in the blob and the variant discriminant as the last byte, Arbitrary
/// will:
/// 1. Read the variant byte (last byte) → picks BAM (0) or CRAM (2)
/// 2. Read lengths for Vec fields from the end → gets portions of the data
/// 3. Read byte content from the front → gets the actual file bytes
///
/// Even if the exact split isn't perfect, the fuzzer will mutate the length
/// bytes to find valid combinations. The key value is having all the right
/// byte patterns in the seed.
fn build_seed(variant: u8, data1: &[u8], data2: &[u8], fai: &[u8], gzi: &[u8]) -> Vec<u8> {
    let mut buf = Vec::new();
    // Front: all the file data the fuzzer needs to discover
    buf.extend_from_slice(data1);
    buf.extend_from_slice(data2);
    buf.extend_from_slice(fai);
    buf.extend_from_slice(gzi);
    // End: variant discriminant (Arbitrary reads this last)
    buf.push(variant % 3);
    buf
}

/// Try to parse a seed and report what Arbitrary extracted.
fn report_parse(seed: &[u8], label: &str) {
    let mut u = Unstructured::new(seed);
    match u.arbitrary::<Input>() {
        Ok(input) => {
            let variant = match &input.readers {
                ReadersInput::Bam { bam, bai } => {
                    format!("Bam(bam={}, bai={})", bam.len(), bai.len())
                }
                ReadersInput::Sam { sam, sai } => {
                    format!("Sam(sam={}, sai={})", sam.len(), sai.len())
                }
                ReadersInput::Cram { cram, crai } => {
                    format!("Cram(cram={}, crai={})", cram.len(), crai.len())
                }
            };
            println!(
                "  {label}: {variant}, fasta_gz={}, fai={}, gzi={}",
                input.alignment.fasta_gz.len(),
                input.alignment.fai.len(),
                input.alignment.gzi.len(),
            );
        }
        Err(e) => {
            println!("  {label}: parse incomplete ({e}) — fuzzer will refine via mutation");
        }
    }
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

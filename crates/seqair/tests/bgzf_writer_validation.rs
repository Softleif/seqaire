//! BGZF writer integration tests: validate seqair's BGZF output is
//! readable by htslib tools (bgzip, samtools).
#![allow(
    clippy::unwrap_used,
    clippy::expect_used,
    clippy::panic,
    clippy::indexing_slicing,
    clippy::arithmetic_side_effects,
    reason = "test code"
)]

use seqair::io::BgzfWriter;
use std::process::Command;

/// Write BGZF data, decompress with bgzip -d, verify content matches.
#[test]
fn bgzip_decompresses_seqair_output() {
    let dir = tempfile::tempdir().unwrap();
    let bgzf_path = dir.path().join("test.gz");
    let expected = b"Hello, BGZF world! This is a test of BGZF compression.\n";

    {
        let file = std::fs::File::create(&bgzf_path).unwrap();
        let mut writer = BgzfWriter::new(file);
        writer.write_all(expected).unwrap();
        writer.finish().unwrap();
    }

    // Decompress with bgzip
    let output =
        Command::new("bgzip").args(["-d", "-c"]).arg(&bgzf_path).output().expect("bgzip not found");
    assert!(output.status.success(), "bgzip -d failed");
    assert_eq!(output.stdout, expected);
}

/// Large BGZF output spanning many blocks, decompressed by bgzip.
#[test]
fn bgzip_decompresses_large_output() {
    let dir = tempfile::tempdir().unwrap();
    let bgzf_path = dir.path().join("large.gz");

    // 500KB of patterned data — forces ~8 BGZF blocks
    let data: Vec<u8> = (0..500_000).map(|i| (i % 251) as u8).collect();

    {
        let file = std::fs::File::create(&bgzf_path).unwrap();
        let mut writer = BgzfWriter::new(file);
        writer.write_all(&data).unwrap();
        writer.finish().unwrap();
    }

    let output =
        Command::new("bgzip").args(["-d", "-c"]).arg(&bgzf_path).output().expect("bgzip not found");
    assert!(output.status.success(), "bgzip -d failed");
    assert_eq!(output.stdout.len(), data.len());
    assert_eq!(output.stdout, data);
}

/// BGZF written in many small chunks (simulating record-by-record writes).
#[test]
fn bgzip_decompresses_chunked_writes() {
    let dir = tempfile::tempdir().unwrap();
    let bgzf_path = dir.path().join("chunked.gz");

    let mut expected = Vec::new();
    {
        let file = std::fs::File::create(&bgzf_path).unwrap();
        let mut writer = BgzfWriter::new(file);

        // Write 1000 small chunks
        for i in 0..1000u32 {
            let chunk = format!("record_{i:04}\t{}\tACGTACGT\n", i * 100);
            writer.write_all(chunk.as_bytes()).unwrap();
            expected.extend_from_slice(chunk.as_bytes());
        }
        writer.finish().unwrap();
    }

    let output =
        Command::new("bgzip").args(["-d", "-c"]).arg(&bgzf_path).output().expect("bgzip not found");
    assert!(output.status.success(), "bgzip -d failed");
    assert_eq!(output.stdout, expected);
}

/// Empty BGZF file (just EOF block) is valid.
#[test]
fn bgzip_accepts_empty_output() {
    let dir = tempfile::tempdir().unwrap();
    let bgzf_path = dir.path().join("empty.gz");

    {
        let file = std::fs::File::create(&bgzf_path).unwrap();
        let writer = BgzfWriter::new(file);
        writer.finish().unwrap();
    }

    let output =
        Command::new("bgzip").args(["-d", "-c"]).arg(&bgzf_path).output().expect("bgzip not found");
    assert!(output.status.success(), "bgzip -d should accept empty BGZF");
    assert!(output.stdout.is_empty(), "decompressed empty BGZF should be empty");
}

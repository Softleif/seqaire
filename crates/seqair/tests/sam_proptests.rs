//! Property-based tests for SAM parsing internals.
#![allow(clippy::unwrap_used, clippy::expect_used, clippy::panic, clippy::indexing_slicing)]
use proptest::prelude::*;
use seqair_types::Base;

// We can't directly call the private parse functions, but we can test the
// public API with generated inputs. For CIGAR, we generate valid CIGAR strings,
// create SAM records, push through the reader, and verify round-trip properties.

/// Generate a valid CIGAR string with 1-10 operations.
fn arb_cigar() -> impl Strategy<Value = String> {
    prop::collection::vec(
        (1..200u32, prop::sample::select(vec!['M', 'I', 'D', 'S', 'N', '=', 'X'])),
        1..10,
    )
    .prop_map(|ops| ops.iter().map(|(len, op)| format!("{len}{op}")).collect::<String>())
}

/// Compute query-consuming length from a CIGAR string.
fn cigar_query_len(cigar: &str) -> usize {
    let mut len = 0;
    let mut num_start = 0;
    for (i, c) in cigar.char_indices() {
        if c.is_ascii_digit() {
            continue;
        }
        let n: usize = cigar[num_start..i].parse().unwrap_or(0);
        match c {
            'M' | 'I' | 'S' | '=' | 'X' => len += n,
            _ => {}
        }
        num_start = i + c.len_utf8();
    }
    len
}

/// Compute reference-consuming length from a CIGAR string.
fn cigar_ref_len(cigar: &str) -> usize {
    let mut len = 0;
    let mut num_start = 0;
    for (i, c) in cigar.char_indices() {
        if c.is_ascii_digit() {
            continue;
        }
        let n: usize = cigar[num_start..i].parse().unwrap_or(0);
        match c {
            'M' | 'D' | 'N' | '=' | 'X' => len += n,
            _ => {}
        }
        num_start = i + c.len_utf8();
    }
    len
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(200))]

    // r[verify sam.record.cigar_parse]
    // r[verify sam.edge.long_cigar]
    #[test]
    fn cigar_roundtrip_preserves_end_pos(
        cigar in arb_cigar(),
        pos in 100i64..10000,
    ) {
        let query_len = cigar_query_len(&cigar);
        if query_len == 0 {
            return Ok(());
        }

        let ref_len = cigar_ref_len(&cigar);
        let expected_end_pos = if ref_len > 0 {
            pos + ref_len as i64 - 1
        } else {
            pos
        };

        let seq: String = (0..query_len).map(|_| 'A').collect();
        let qual: String = (0..query_len).map(|_| 'I').collect();
        let sam_text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
             read1\t99\tchr1\t{pos_1}\t60\t{cigar}\t=\t200\t150\t{seq}\t{qual}\n",
            pos_1 = pos + 1, // SAM is 1-based
        );

        let dir = tempfile::tempdir()?;
        let sam_gz = create_indexed_sam_proptest(dir.path(), &sam_text);
        let mut reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz)?;
        let tid = reader.header().tid("chr1").unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, 1, 100_000_000, &mut store)?;

        prop_assert!(!store.is_empty(), "should fetch at least 1 record");
        let rec = store.record(0);
        prop_assert_eq!(rec.pos, pos);
        prop_assert_eq!(rec.end_pos, expected_end_pos);
        prop_assert_eq!(rec.seq_len as usize, query_len);
    }

    // r[verify sam.record.seq_decode]
    #[test]
    fn seq_decode_roundtrip(
        bases in prop::collection::vec(prop::sample::select(vec![b'A', b'C', b'G', b'T']), 1..200),
    ) {
        let seq_str: String = bases.iter().map(|&b| b as char).collect();
        let qual_str: String = bases.iter().map(|_| 'I').collect();
        let cigar = format!("{}M", bases.len());
        let sam_text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
             read1\t99\tchr1\t100\t60\t{cigar}\t=\t200\t150\t{seq_str}\t{qual_str}\n"
        );

        let dir = tempfile::tempdir()?;
        let sam_gz = create_indexed_sam_proptest(dir.path(), &sam_text);
        let mut reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz)?;
        let tid = reader.header().tid("chr1").unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, 1, 100_000_000, &mut store)?;

        prop_assert_eq!(store.len(), 1);
        let seq = store.seq(0);
        prop_assert_eq!(seq.len(), bases.len());
        for (i, (&expected, &actual)) in bases.iter().zip(seq.iter()).enumerate() {
            let expected_base = match expected {
                b'A' => Base::A,
                b'C' => Base::C,
                b'G' => Base::G,
                b'T' => Base::T,
                _ => Base::Unknown,
            };
            prop_assert_eq!(actual, expected_base, "pos {}: expected {:?}", i, expected_base);
        }
    }

    // r[verify sam.record.qual_decode]
    #[test]
    fn qual_decode_roundtrip(
        // Exclude quality 9: Phred+33 encodes it as ASCII '*' (42), which SAM
        // treats as "quality unavailable" when it is the entire QUAL field.
        quals in prop::collection::vec(
            (0u8..=40).prop_filter("not SAM unknown-qual sentinel", |&q| q + 33 != b'*'),
            1..200,
        ),
    ) {
        let seq_str: String = quals.iter().map(|_| 'A').collect();
        let qual_str: String = quals.iter().map(|&q| (q + 33) as char).collect();
        let cigar = format!("{}M", quals.len());
        let sam_text = format!(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:100000000\n\
             read1\t99\tchr1\t100\t60\t{cigar}\t=\t200\t150\t{seq_str}\t{qual_str}\n"
        );

        let dir = tempfile::tempdir()?;
        let sam_gz = create_indexed_sam_proptest(dir.path(), &sam_text);
        let mut reader = seqair::sam::reader::IndexedSamReader::open(&sam_gz)?;
        let tid = reader.header().tid("chr1").unwrap();
        let mut store = seqair::bam::RecordStore::new();
        reader.fetch_into(tid, 1, 100_000_000, &mut store)?;

        prop_assert_eq!(store.len(), 1);
        let stored_qual = store.qual(0);
        prop_assert_eq!(stored_qual.len(), quals.len());
        for (i, (&expected, &actual)) in quals.iter().zip(stored_qual.iter()).enumerate() {
            prop_assert_eq!(actual, expected, "pos {}", i);
        }
    }
}

fn create_indexed_sam_proptest(dir: &std::path::Path, sam_text: &str) -> std::path::PathBuf {
    use std::process::Command;

    let sam_path = dir.join("test.sam");
    let sam_gz_path = dir.join("test.sam.gz");

    std::fs::write(&sam_path, sam_text).expect("write SAM");

    let status = Command::new("bgzip")
        .arg("-c")
        .arg(&sam_path)
        .stdout(std::fs::File::create(&sam_gz_path).expect("create .sam.gz"))
        .status()
        .expect("bgzip not found");
    assert!(status.success(), "bgzip failed");

    let status = Command::new("tabix")
        .args(["-p", "sam"])
        .arg(&sam_gz_path)
        .status()
        .expect("tabix not found");
    assert!(status.success(), "tabix failed");

    sam_gz_path
}
